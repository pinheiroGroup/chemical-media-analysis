"""
Downstream ML: link BW25113 medium chemistry to KinBiont aHPM growth parameters.

Workflow:
  1. Load batch fit results (results/batch_fit_results.csv)
  2. Load GrowthDataEvaluation.xlsx → map curve labels to Condition IDs + reference K/r
  3. Load Medium composition.xlsx   → 44 compound concentrations per condition
  4. Join and aggregate: mean growth params per condition
  5. Correlation analysis: Spearman rank correlations (compounds → params)
  6. Decision-tree regression: feature importance for gr, N_max, exit_lag_rate
  7. Save results to results/ml_results/

Input:
  results/batch_fit_results.csv
  ../chemical-media-dataset/xlsx_raw/BW25113_GrowthDataEvaluation.xlsx
  ../chemical-media-dataset/xlsx_raw/BW25113_Medium composition.xlsx

Output:
  results/ml_results/condition_means.csv       (per-condition mean params)
  results/ml_results/correlations.csv          (Spearman ρ, compounds × params)
  results/ml_results/feature_importance_gr.csv
  results/ml_results/feature_importance_Nmax.csv
  results/ml_results/feature_importance_lag.csv
  results/ml_results/validation_kinbiont_vs_reference.csv (gr vs r, N_max vs K)
"""

using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames
using XLSX
using Statistics: mean, cor
using StatsBase: corspearman
using DecisionTree

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
    # Quick Spearman correlations between fitted and reference
    gr_ok  = .!isnan.(val_df.gr) .& .!isnan.(val_df.r_ref)
    Nmax_ok = .!isnan.(val_df.N_max) .& .!isnan.(val_df.K_ref)
    println("  ρ(gr, r_ref)   = $(round(corspearman(val_df.gr[gr_ok],   val_df.r_ref[gr_ok]);   digits=3))")
    println("  ρ(N_max, K_ref) = $(round(corspearman(val_df.N_max[Nmax_ok], val_df.K_ref[Nmax_ok]); digits=3))")

    println("\nSpearman correlations (compounds → growth params) …")
    corr_df = correlation_analysis(cond_df, compound_cols)
    CSV.write(joinpath(ML_DIR, "correlations.csv"), corr_df)
    println("  Top 5 correlated compounds with gr:")
    top5 = sort(corr_df, :gr; rev=true)[1:min(5, nrow(corr_df)), :]
    for r in eachrow(top5)
        println("    $(rpad(r.compound, 30))  ρ=$(round(r.gr; digits=3))")
    end

    println("\nDecision-tree feature importance …")
    for (param, fname) in [(:gr, "gr"), (:N_max, "Nmax"), (:exit_lag_rate, "lag")]
        imp_df = tree_importance(cond_df, compound_cols, param)
        CSV.write(joinpath(ML_DIR, "feature_importance_$(fname).csv"), imp_df)
        println("  Top-3 predictors of $(param):")
        for r in eachrow(imp_df[1:min(3, nrow(imp_df)), :])
            println("    $(rpad(r.compound, 30))  importance=$(round(r.importance; digits=4))")
        end
    end

    println("\nAll results saved to $(ML_DIR)")
end

main()
