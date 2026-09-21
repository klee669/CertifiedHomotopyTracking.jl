include(joinpath(@__DIR__, "..", "common", "higher_order_predictor.jl"))
include(joinpath(@__DIR__, "..", "common", "benchmark_systems.jl"))
include(joinpath(@__DIR__, "..", "common", "csv_io.jl"))

using Statistics

const INITIAL_RADIUS = 0.1

function run_benchmarks(; full=false, max_paths=nothing, gamma_mode=:diagonal_start)
    dimensions = full ? (3:6) : (3:3)
    rows = NamedTuple[]
    for spec in benchmark_specs(; degree=2, dimensions)
        H, roots = make_benchmark(spec; gamma_mode)
        count = isnothing(max_paths) ? (full ? length(roots) : 1) : min(max_paths, length(roots))
        path_stats = NamedTuple[]
        println("$(system_name(spec)): tracking $count/$(length(roots)) paths")
        for path_index in 1:count
            _, stats = tracking_constant_apriori(
                H, copy(roots[path_index]), INITIAL_RADIUS;
                iterations_count=true, show_display=false, final_refine=false,
            )
            push!(path_stats, stats)
            println("  path $path_index: $(stats.iterations) iterations")
        end
        push!(rows, (
            system=system_name(spec),
            paths=count,
            avg_iterations=mean(s.iterations for s in path_stats),
            avg_min_dt=mean(s.min_dt for s in path_stats),
            avg_median_dt=mean(s.median_dt for s in path_stats),
            avg_min_radius=mean(s.min_radius for s in path_stats),
            avg_eta=mean(s.max_eta for s in path_stats),
            seed=spec.seed,
            gamma_mode=String(gamma_mode),
        ))
    end
    output = joinpath(@__DIR__, "results", "generated", "benchmark_summary.csv")
    write_csv(output, rows)
    println("wrote $output")
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    full = get(ENV, "CHT_FULL", "0") == "1"
    max_paths_text = get(ENV, "CHT_MAX_PATHS", "")
    max_paths = isempty(max_paths_text) ? nothing : parse(Int, max_paths_text)
    gamma_mode = Symbol(get(ENV, "CHT_GAMMA_MODE", "diagonal_start"))
    run_benchmarks(; full, max_paths, gamma_mode)
end
