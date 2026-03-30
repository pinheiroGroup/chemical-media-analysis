"""
Downstream ML: link BW25113 medium chemistry to KinBiont aHPM growth parameters.

Workflow:
  1. Load batch fit results (results/batch_fit_results.csv)
  2. Load GrowthDataEvaluation.xlsx → map curve labels to Condition IDs + reference K/r
  3. Load Medium composition.xlsx   → 44 compound concentrations per condition
  4. Join and aggregate: mean growth params per condition
  5. Correlation analysis: Spearman rank correlations (compounds → params)
  6. Decision-tree regression: feature importance for gr, N_max, exit_lag_rate
  7. Figures: validation scatter, correlation heatmap, feature importance bars
  8. Save results to results/ml_results/

Input:
  results/batch_fit_results.csv
  ../chemical-media-dataset/xlsx_raw/BW25113_GrowthDataEvaluation.xlsx
  ../chemical-media-dataset/xlsx_raw/BW25113_Medium composition.xlsx

Output (CSVs):
  results/ml_results/condition_means.csv
  results/ml_results/correlations.csv
  results/ml_results/feature_importance_{gr,Nmax,lag}.csv
  results/ml_results/validation_kinbiont_vs_reference.csv
Output (figures, 300 dpi PDF + PNG):
  results/ml_results/fig_validation.pdf
  results/ml_results/fig_correlation_heatmap.pdf
  results/ml_results/fig_feature_importance.pdf
"""

using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames
using XLSX
using Statistics: mean, cor
using StatsBase: corspearman
using DecisionTree
using CairoMakie

const XLSX_DIR    = "/media/aivuk/64fce268-2613-4033-b39f-537ae2d28805/roms/pinheiroTech/chemical-media-dataset/xlsx_raw"
const RESULTS_DIR = joinpath(@__DIR__, "results")
const ML_DIR      = joinpath(RESULTS_DIR, "ml_results")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function read_xlsx_as_df(path::String, sheet::Int=1)::DataFrame
    xf   = XLSX.readxlsx(path)
    sh   = xf[XLSX.sheetnames(xf)[sheet]]
    data = XLSX.getdata(sh)
    # First row is header
    header = [string(x) for x in data[1, :]]
    rows   = [data[i, :] for i in 2:size(data, 1)]
    df     = DataFrame([col => [r[j] for r in rows] for (j, col) in enumerate(header)])
    return df
end

function safe_float(x)::Float64
    x === nothing || ismissing(x) && return NaN
    try
        Float64(x)
    catch
        NaN
    end
end

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

function load_data()
    println("Loading batch fit results …")
    fit_df = CSV.read(joinpath(RESULTS_DIR, "batch_fit_results.csv"), DataFrame)
    filter!(r -> r.converged, fit_df)
    println("  $(nrow(fit_df)) converged fits")

    println("Loading GrowthDataEvaluation …")
    eval_df = read_xlsx_as_df(joinpath(XLSX_DIR, "BW25113_GrowthDataEvaluation.xlsx"))
    rename!(eval_df, strip.(names(eval_df)))
    println("  $(nrow(eval_df)) rows, columns: $(names(eval_df))")

    println("Loading Medium composition …")
    comp_df = read_xlsx_as_df(joinpath(XLSX_DIR, "BW25113_Medium composition.xlsx"))
    rename!(comp_df, strip.(names(comp_df)))
    println("  $(nrow(comp_df)) rows, $(ncol(comp_df)) columns")

    return fit_df, eval_df, comp_df
end

# ---------------------------------------------------------------------------
# Join: fit results → condition ID → compound concentrations
# ---------------------------------------------------------------------------

function build_joined(fit_df, eval_df, comp_df)::DataFrame
    # eval_df has columns: Label, Assay ID, Condition ID, K, r, K_info, r_info
    # Standardise column names (strip spaces, handle variations)
    label_col     = findfirst(n -> occursin("label", lowercase(n)), names(eval_df))
    cond_col      = findfirst(n -> occursin("condition", lowercase(n)), names(eval_df))
    k_col         = findfirst(n -> n == "K" || lowercase(n) == "k", names(eval_df))
    r_col         = findfirst(n -> n == "r" || lowercase(n) == "r", names(eval_df))

    eval_df[!, :Label_str]    = string.(eval_df[!, names(eval_df)[label_col]])
    eval_df[!, :ConditionID]  = string.(eval_df[!, names(eval_df)[cond_col]])
    eval_df[!, :K_ref]        = safe_float.(eval_df[!, names(eval_df)[k_col]])
    eval_df[!, :r_ref]        = safe_float.(eval_df[!, names(eval_df)[r_col]])

    # Merge fit results with condition IDs
    joined = leftjoin(fit_df, select(eval_df, :Label_str => :label, :ConditionID, :K_ref, :r_ref);
                      on = :label)
    dropmissing!(joined, :ConditionID)
    println("  Joined: $(nrow(joined)) curves with condition IDs")

    # comp_df first column is condition identifier
    cond_col_comp = findfirst(n -> occursin("condition", lowercase(n)), names(comp_df))
    comp_df[!, :ConditionID] = string.(comp_df[!, names(comp_df)[cond_col_comp]])

    # Compound columns: all after the first 3 metadata columns
    meta_cols = names(comp_df)[1:3]
    compound_cols = setdiff(names(comp_df), meta_cols)
    # Convert compound columns to Float64
    for col in compound_cols
        comp_df[!, col] = safe_float.(comp_df[!, col])
    end

    final = leftjoin(joined, select(comp_df, :ConditionID, compound_cols...); on = :ConditionID)
    dropmissing!(final, compound_cols[1])   # drop rows with no composition info
    println("  Final joined: $(nrow(final)) curves with composition data")
    return final, compound_cols
end

# ---------------------------------------------------------------------------
# Aggregate per condition
# ---------------------------------------------------------------------------

function aggregate_per_condition(joined::DataFrame, compound_cols::Vector{String})::DataFrame
    param_cols = [:gr, :exit_lag_rate, :N_max, :shape]
    group_cols = vcat([:ConditionID], Symbol.(compound_cols))

    agg = combine(groupby(joined, :ConditionID)) do g
        row = Dict{Symbol, Any}(:ConditionID => g.ConditionID[1])
        for col in param_cols
            row[col] = mean(filter(!isnan, g[!, col]))
        end
        row[:K_ref] = mean(filter(!isnan, g[!, :K_ref]))
        row[:r_ref] = mean(filter(!isnan, g[!, :r_ref]))
        # compound concentrations are constant within a condition
        for col in compound_cols
            row[Symbol(col)] = g[1, col]
        end
        DataFrame(row)
    end
    println("  Aggregated: $(nrow(agg)) unique conditions")
    return agg
end

# ---------------------------------------------------------------------------
# Spearman correlation: compounds → growth params
# ---------------------------------------------------------------------------

function correlation_analysis(cond_df::DataFrame, compound_cols::Vector{String})::DataFrame
    params = [:gr, :exit_lag_rate, :N_max, :shape]
    X = Matrix{Float64}(cond_df[!, compound_cols])
    rows = map(compound_cols) do comp
        row = Dict{Symbol,Any}(:compound => comp)
        for p in params
            y = Float64.(cond_df[!, p])
            mask = .!isnan.(y) .& .!isnan.(X[:, findfirst(==(comp), compound_cols)])
            rho = length(mask) >= 3 ? corspearman(X[mask, findfirst(==(comp), compound_cols)], y[mask]) : NaN
            row[p] = rho
        end
        row
    end
    return DataFrame(rows)
end

# ---------------------------------------------------------------------------
# Decision tree feature importance
# ---------------------------------------------------------------------------

function tree_importance(cond_df::DataFrame, compound_cols::Vector{String},
                         target::Symbol; max_depth=5)::DataFrame
    X = Matrix{Float64}(cond_df[!, compound_cols])
    y = Float64.(cond_df[!, target])
    mask = .!isnan.(y) .& all(.!isnan.(X), dims=2)[:]
    Xm, ym = X[mask, :], y[mask]
    nrow(DataFrame(Xm, :auto)) < 10 && return DataFrame(compound=compound_cols, importance=fill(NaN, length(compound_cols)))

    model = build_forest(ym, Xm; n_trees=100, max_depth=max_depth, rng=42)
    imp   = impurity_importance(model)
    return sort(DataFrame(compound=compound_cols, importance=imp), :importance; rev=true)
end

# ---------------------------------------------------------------------------
# Validation: KinBiont gr vs reference r, N_max vs K
# ---------------------------------------------------------------------------

function validation_comparison(cond_df::DataFrame)::DataFrame
    select(cond_df, :ConditionID, :gr, :r_ref, :N_max, :K_ref)
end

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

function fig_validation(val_df::DataFrame, rho_gr::Float64, rho_Nmax::Float64)
    fig = Figure(size=(900, 420), fontsize=13)

    for (col, (x_col, y_col, xlabel, ylabel, rho)) in enumerate([
        (:gr,   :r_ref,  "Reference r (h⁻¹)",          "Fitted gr (h⁻¹)",    rho_gr),
        (:N_max, :K_ref, "Reference K (OD)",             "Fitted N_max (OD)", rho_Nmax),
    ])
        ax = Axis(fig[1, col];
            xlabel = xlabel,
            ylabel = ylabel,
            title  = "ρ = $(round(rho; digits=3))",
        )
        x = val_df[!, x_col]; y = val_df[!, y_col]
        ok = .!isnan.(x) .& .!isnan.(y)
        scatter!(ax, x[ok], y[ok]; color=(:steelblue, 0.5), markersize=5)
        # y = x reference line
        lo = min(minimum(x[ok]), minimum(y[ok]))
        hi = max(maximum(x[ok]), maximum(y[ok]))
        lines!(ax, [lo, hi], [lo, hi]; color=:black, linestyle=:dash, linewidth=1)
    end

    return fig
end

function fig_correlation_heatmap(corr_df::DataFrame)
    params      = [:gr, :exit_lag_rate, :N_max, :shape]
    param_labels = ["gr", "exit_lag_rate", "N_max", "shape"]
    compounds   = corr_df.compound
    Z = [corr_df[i, p] for i in 1:nrow(corr_df), p in params]

    # Sort compounds by |ρ| with gr (most informative first)
    order = sortperm(abs.(Z[:, 1]); rev=true)
    Z_sorted = Z[order, :]
    comp_sorted = compounds[order]

    n_comp = length(comp_sorted)
    fig = Figure(size=(520, max(300, 18 * n_comp + 80)), fontsize=11)
    ax  = Axis(fig[1, 1];
        xticks = (1:length(params), param_labels),
        yticks = (1:n_comp, comp_sorted),
        xticklabelrotation = 0.0,
    )
    hm = heatmap!(ax, 1:length(params), 1:n_comp, Z_sorted';
        colormap = :RdBu, colorrange = (-1, 1))
    Colorbar(fig[1, 2], hm; label = "Spearman ρ", width = 15)

    return fig
end

function fig_feature_importance(imp_dfs::Vector{DataFrame}, targets::Vector{String};
                                 top_n::Int = 15)
    n = length(imp_dfs)
    fig = Figure(size=(360 * n, 420), fontsize=12)

    for (col, (imp_df, target)) in enumerate(zip(imp_dfs, targets))
        df  = imp_df[1:min(top_n, nrow(imp_df)), :]
        ax  = Axis(fig[1, col];
            xlabel = "Importance",
            title  = target,
            yticks = (1:nrow(df), reverse(df.compound)),
        )
        barplot!(ax, 1:nrow(df), reverse(df.importance);
            direction = :x, color = :steelblue)
    end

    return fig
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    mkpath(ML_DIR)

    fit_df, eval_df, comp_df = load_data()

    println("\nJoining datasets …")
    joined, compound_cols = build_joined(fit_df, eval_df, comp_df)
    CSV.write(joinpath(ML_DIR, "joined_curves.csv"), joined)

    println("\nAggregating per condition …")
    cond_df = aggregate_per_condition(joined, compound_cols)
    CSV.write(joinpath(ML_DIR, "condition_means.csv"), cond_df)

    println("\nValidation: KinBiont vs reference K/r …")
    val_df = validation_comparison(cond_df)
    CSV.write(joinpath(ML_DIR, "validation_kinbiont_vs_reference.csv"), val_df)
    gr_ok   = .!isnan.(val_df.gr)   .& .!isnan.(val_df.r_ref)
    Nmax_ok = .!isnan.(val_df.N_max) .& .!isnan.(val_df.K_ref)
    rho_gr   = corspearman(val_df.gr[gr_ok],    val_df.r_ref[gr_ok])
    rho_Nmax = corspearman(val_df.N_max[Nmax_ok], val_df.K_ref[Nmax_ok])
    println("  ρ(gr, r_ref)    = $(round(rho_gr;   digits=3))")
    println("  ρ(N_max, K_ref) = $(round(rho_Nmax; digits=3))")

    println("\nSpearman correlations (compounds → growth params) …")
    corr_df = correlation_analysis(cond_df, compound_cols)
    CSV.write(joinpath(ML_DIR, "correlations.csv"), corr_df)
    println("  Top 5 correlated compounds with gr:")
    top5 = sort(corr_df, :gr; rev=true)[1:min(5, nrow(corr_df)), :]
    for r in eachrow(top5)
        println("    $(rpad(r.compound, 30))  ρ=$(round(r.gr; digits=3))")
    end

    println("\nDecision-tree feature importance …")
    imp_dfs  = DataFrame[]
    imp_names = [("gr", "gr"), ("Nmax", "N_max"), ("lag", "exit_lag_rate")]
    for (fname, param_str) in imp_names
        param  = Symbol(param_str)
        imp_df = tree_importance(cond_df, compound_cols, param)
        CSV.write(joinpath(ML_DIR, "feature_importance_$(fname).csv"), imp_df)
        push!(imp_dfs, imp_df)
        println("  Top-3 predictors of $(param):")
        for r in eachrow(imp_df[1:min(3, nrow(imp_df)), :])
            println("    $(rpad(r.compound, 30))  importance=$(round(r.importance; digits=4))")
        end
    end

    println("\nGenerating figures …")

    f1 = fig_validation(val_df, rho_gr, rho_Nmax)
    save(joinpath(ML_DIR, "fig_validation.pdf"), f1; pt_per_unit=1)
    save(joinpath(ML_DIR, "fig_validation.png"), f1; px_per_unit=2)
    println("  fig_validation saved")

    f2 = fig_correlation_heatmap(corr_df)
    save(joinpath(ML_DIR, "fig_correlation_heatmap.pdf"), f2; pt_per_unit=1)
    save(joinpath(ML_DIR, "fig_correlation_heatmap.png"), f2; px_per_unit=2)
    println("  fig_correlation_heatmap saved")

    f3 = fig_feature_importance(imp_dfs, ["gr", "N_max", "exit_lag_rate"])
    save(joinpath(ML_DIR, "fig_feature_importance.pdf"), f3; pt_per_unit=1)
    save(joinpath(ML_DIR, "fig_feature_importance.png"), f3; px_per_unit=2)
    println("  fig_feature_importance saved")

    println("\nAll results saved to $(ML_DIR)")
end

main()
