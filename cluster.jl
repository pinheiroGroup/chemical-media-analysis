"""
Cluster all 13,608 BW25113 chemical-media growth curves with KinBiont.

Workflow:
  1. Load all 7 rounds from data/Round0{1-7}.csv
  2. Build a single data matrix (curves × timepoints), handling missing values
  3. Run WCSS elbow sweep (k=2..10) using KinBiont clustering pipeline
     (z-scored k-means, trend test, constant pre-screening)
  4. Cluster with optimal k
  5. Compare cluster assignments with medium composition metadata

Run:
  julia --threads auto --project=. cluster.jl

Input:  data/Round0{1-7}.csv
        ../chemical-media-dataset/xlsx_raw/BW25113_Medium\ composition.xlsx (optional)

Output:
  results/clustering/cluster_assignments.csv  — label, round, cluster
  results/clustering/wcss_sweep.csv           — k, wcss
  results/clustering/centroids.csv            — cluster, timepoint, mean_od
  results/clustering/fig_elbow.pdf/png
  results/clustering/fig_centroids.pdf/png
"""

using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames
using Statistics: mean, std
using Kinbiont
using CairoMakie

const DATA_DIR    = joinpath(@__DIR__, "data")
const RESULTS_DIR = joinpath(@__DIR__, "results", "clustering")

# ---------------------------------------------------------------------------
# Load all rounds into a unified time × curve matrix
# ---------------------------------------------------------------------------

function load_all_rounds()::Tuple{Vector{Float64}, DataFrame}
    println("Loading rounds 1–7…")
    all_dfs = DataFrame[]
    for r in 1:7
        path = joinpath(DATA_DIR, "Round$(lpad(r, 2, '0')).csv")
        isfile(path) || (println("  Round $r not found, skipping"); continue)
        df = CSV.read(path, DataFrame; missingstring="")
        println("  Round $r: $(ncol(df)-1) curves, $(nrow(df)) time points")
        push!(all_dfs, df)
    end
    isempty(all_dfs) && error("No data files found in $DATA_DIR")

    # All rounds share the same time column; use the first
    times = Float64.(all_dfs[1][!, :Time_h])

    # Trim all rounds to common length
    n_tp = minimum(nrow(df) for df in all_dfs)
    times = times[1:n_tp]

    # Build label → OD vector mapping
    label_od = Dict{String, Vector{Float64}}()
    for df in all_dfs
        for col in names(df)[2:end]
            od_raw = df[1:n_tp, col]
            # Replace missing with NaN, then NaN-fill with column mean
            od = [ismissing(x) ? NaN : Float64(x) for x in od_raw]
            label_od[col] = od
        end
    end

    labels   = sort(collect(keys(label_od)))
    println("Total: $(length(labels)) curves, $n_tp time points")
    return times, DataFrame(:label => labels,
        [Symbol("t$i") => [label_od[l][i] for l in labels] for i in 1:n_tp]...)
end

# ---------------------------------------------------------------------------
# Build clustering matrix: replace NaN with column mean
# ---------------------------------------------------------------------------

function build_matrix(df::DataFrame, n_tp::Int)::Matrix{Float64}
    mat = Matrix{Float64}(df[!, 2:end])   # n_curves × n_tp
    for j in 1:n_tp
        col = mat[:, j]
        valid = filter(!isnan, col)
        fill_val = isempty(valid) ? 0.0 : mean(valid)
        mat[isnan.(col), j] .= fill_val
    end
    return mat
end

function find_elbow(ks, wcss_vals)
    length(ks) < 3 && return ks[end]
    d2 = diff(diff(wcss_vals))
    return ks[argmax(d2) + 1]
end

# ---------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------

function run_clustering(mat::Matrix{Float64}, times::Vector{Float64},
                        labels::Vector{String}; k_range=2:10)
    gd = GrowthData(mat, times, labels)

    println("\nWCSS sweep k=$(first(k_range))..$(last(k_range))…")
    ks        = collect(k_range)
    wcss_vals = Float64[]
    for k in ks
        opts = FitOptions(
            cluster                    = true,
            n_clusters                 = k,
            cluster_prescreen_constant = true,
            cluster_tol_const          = 1.5,
        )
        proc = preprocess(gd, opts)
        push!(wcss_vals, something(proc.wcss, 0.0))
        print("  k=$k  WCSS=$(round(wcss_vals[end]; digits=1))\n")
    end

    opt_k = find_elbow(ks, wcss_vals)
    println("Optimal k = $opt_k")

    println("\nFinal clustering with k=$opt_k…")
    opts_final = FitOptions(
        cluster                    = true,
        n_clusters                 = opt_k,
        cluster_prescreen_constant = true,
        cluster_tol_const          = 1.5,
    )
    proc_final = preprocess(gd, opts_final)
    assignments = something(proc_final.clusters, ones(Int, size(mat, 1)))

    return (ks=ks, wcss=wcss_vals, optimal_k=opt_k,
            assignments=assignments, centroids=proc_final.centroids)
end

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

function fig_elbow(ks, wcss_vals, opt_k)
    fig = Figure(size=(520, 360), fontsize=13)
    ax  = Axis(fig[1,1];
        xlabel = "Number of clusters (k)",
        ylabel = "WCSS",
        xticks = ks,
        title  = "KinBiont k-means sweep — BW25113 chemical media",
    )
    lines!(ax,  ks, wcss_vals; color=:steelblue, linewidth=2)
    scatter!(ax, ks, wcss_vals; color=:steelblue, markersize=7)
    idx = findfirst(==(opt_k), ks)
    scatter!(ax, [opt_k], [wcss_vals[idx]];
             color=:steelblue, markersize=16, marker=:star5,
             label="Optimal k=$opt_k")
    axislegend(ax; position=:rt, framevisible=false)
    return fig
end

function fig_centroids(centroids, times, opt_k, assignments)
    counts = [sum(assignments .== k) for k in 1:opt_k]
    palette = cgrad(:tab10, opt_k; categorical=true)

    fig = Figure(size=(700, 420), fontsize=12)
    ax  = Axis(fig[1,1];
        xlabel = "Time (h)",
        ylabel = "OD (z-scored)",
        title  = "Cluster centroids (z-scored growth shapes)",
    )
    for k in 1:opt_k
        c = centroids[k, :]
        lines!(ax, times, c; color=palette[k], linewidth=2,
               label="Cluster $k (n=$(counts[k]))")
    end
    axislegend(ax; position=:rb, framevisible=false, labelsize=10)
    return fig
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    mkpath(RESULTS_DIR)

    times, curve_df = load_all_rounds()
    n_tp     = ncol(curve_df) - 1
    labels   = curve_df.label
    mat      = build_matrix(curve_df, n_tp)

    cl = run_clustering(mat, times, labels)

    # Cluster size summary
    for k in 1:cl.optimal_k
        n = sum(cl.assignments .== k)
        println("  Cluster $k: $n curves ($(round(100n/length(labels);digits=1))%)")
    end

    # Save WCSS sweep
    CSV.write(joinpath(RESULTS_DIR, "wcss_sweep.csv"),
              DataFrame(k=cl.ks, wcss=cl.wcss))
    println("\nSaved wcss_sweep.csv")

    # Save assignments
    rounds = [split(l, r"(?<=Curve)")[1] for l in labels]   # best-effort round label
    asgn_df = DataFrame(label=labels, cluster=cl.assignments)
    CSV.write(joinpath(RESULTS_DIR, "cluster_assignments.csv"), asgn_df)
    println("Saved cluster_assignments.csv")

    # Save centroids (long format)
    cent_rows = []
    for k in 1:cl.optimal_k, (i, t) in enumerate(times)
        push!(cent_rows, (cluster=k, time_h=t, mean_od_z=cl.centroids[k, i]))
    end
    CSV.write(joinpath(RESULTS_DIR, "centroids.csv"), DataFrame(cent_rows))
    println("Saved centroids.csv")

    # Figures
    println("\nGenerating figures…")
    f1 = fig_elbow(cl.ks, cl.wcss, cl.optimal_k)
    save(joinpath(RESULTS_DIR, "fig_elbow.pdf"), f1; pt_per_unit=1)
    save(joinpath(RESULTS_DIR, "fig_elbow.png"), f1; px_per_unit=2)
    println("  fig_elbow saved")

    f2 = fig_centroids(cl.centroids, times, cl.optimal_k, cl.assignments)
    save(joinpath(RESULTS_DIR, "fig_centroids.pdf"), f2; pt_per_unit=1)
    save(joinpath(RESULTS_DIR, "fig_centroids.png"), f2; px_per_unit=2)
    println("  fig_centroids saved")

    println("\nAll results in $RESULTS_DIR")
end

main()
