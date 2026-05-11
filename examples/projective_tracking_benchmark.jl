# Run with:
#   julia --project=. examples/projective_tracking_benchmark.jl
#
# Optional arguments:
#   --precision=256
#   --max-roots=8
#   --seed=20260511
#   --trials=1
#   --mode=both        # affine, projective, both
#   --benchmark=both   # highlevel, diagnostic, both
#   --vertices=4
#   --h-init=0.1
#   --max-paths=0      # diagnostic only; 0 means no cap
#   --root-match=heuristic  # heuristic, certified, certified_or_heuristic

using CertifiedHomotopyTracking

const DEFAULT_PRECISION_BITS = 256
const DEFAULT_MAX_ROOTS = 8
const DEFAULT_RANDOM_SEED = 20260511
const DEFAULT_TRIALS = 1
const DEFAULT_MODE = "both"
const DEFAULT_BENCHMARK = "both"
const DEFAULT_NUM_VERTICES = 4
const DEFAULT_H_INIT = 0.1
const DEFAULT_MAX_PATHS = 0
const DEFAULT_ROOT_MATCH = "heuristic"

function parse_options(args)
    options = Dict(
        "precision" => string(DEFAULT_PRECISION_BITS),
        "max-roots" => string(DEFAULT_MAX_ROOTS),
        "seed" => string(DEFAULT_RANDOM_SEED),
        "trials" => string(DEFAULT_TRIALS),
        "mode" => DEFAULT_MODE,
        "benchmark" => DEFAULT_BENCHMARK,
        "vertices" => string(DEFAULT_NUM_VERTICES),
        "h-init" => string(DEFAULT_H_INIT),
        "max-paths" => string(DEFAULT_MAX_PATHS),
        "root-match" => DEFAULT_ROOT_MATCH,
    )

    for arg in args
        if startswith(arg, "--") && occursin("=", arg)
            key, value = split(arg[3:end], "=", limit=2)
            options[key] = value
        else
            error("Unsupported argument: $arg")
        end
    end

    mode = lowercase(options["mode"])
    mode in ("affine", "projective", "both") || error("--mode must be affine, projective, or both.")
    benchmark = lowercase(options["benchmark"])
    benchmark in ("highlevel", "diagnostic", "both") || error("--benchmark must be highlevel, diagnostic, or both.")
    root_match = lowercase(options["root-match"])
    root_match in ("heuristic", "certified", "certified_or_heuristic") || error("--root-match must be heuristic, certified, or certified_or_heuristic.")

    return (
        precision_bits = parse(Int, options["precision"]),
        max_roots = parse(Int, options["max-roots"]),
        seed = parse(UInt64, options["seed"]),
        trials = parse(Int, options["trials"]),
        mode = mode,
        benchmark = benchmark,
        num_vertices = parse(Int, options["vertices"]),
        h_init = parse(Float64, options["h-init"]),
        max_paths = parse(Int, options["max-paths"]),
        root_match = Symbol(root_match),
    )
end

mutable struct LCG
    state::UInt64
end

function next_float!(rng::LCG)
    rng.state = 6364136223846793005 * rng.state + 1442695040888963407
    return Float64(rng.state >> 11) / Float64(UInt64(1) << 53)
end

function unit_complex!(CC, rng::LCG)
    theta = 2pi * next_float!(rng)
    return CC(cos(theta), sin(theta))
end

function benchmark_system()
    @variables x y z λ
    @variables u1 u2 u3

    F = [
        2 * (x - u1) - 6 * λ * (x^2 + y^2)^2 * x,
        2 * (y - u2) - 6 * λ * (x^2 + y^2)^2 * y,
        2 * (z - u3) + 4 * λ * z^3,
        0 * u1 + z^4 - (x^2 + y^2)^3,
    ]

    vars = [x, y, z, λ]
    pars = [u1, u2, u3]
    return F, vars, pars
end

function initial_vertices(CC, pars, num_vertices::Int, seed::UInt64)
    num_vertices >= 1 || error("num_vertices must be positive.")
    rng = LCG(seed)

    bp = [CC(0.09868, 0.675389), CC(0.423238, 0.713082), CC(0.592351, 0.144969)]
    x0 = [
        CC(1.23836, -0.422501),
        CC(1.19574, -1.0474),
        CC(2.08916, 1.85256),
        CC(-0.0126777, 0.0505892),
    ]

    vertices = [vertex(bp, [x0])]
    for _ in 2:num_vertices
        rand_u = [unit_complex!(CC, rng) for _ in 1:length(pars)]
        push!(vertices, vertex(rand_u))
    end
    return vertices
end

function complete_graph!(vertices)
    edges = Edge[]
    for i in 1:length(vertices)-1
        for j in i+1:length(vertices)
            e = edge(vertices[i], vertices[j])
            push!(edges, e)
            push!(vertices[i].Edges, e)
            push!(vertices[j].Edges, e)
        end
    end
    return edges
end

function acb_norm(z)
    return Float64(abs(z))
end

function vector_norm(v)
    isempty(v) && return 0.0
    return maximum(acb_norm.(v))
end

mutable struct BenchmarkMetrics
    mode::String
    trial::Int
    total_time::Float64
    edge_calls::Int
    path_attempts::Int
    successes::Int
    failures::Int
    max_coord_norm::Float64
    iterations::Int
    rejected_steps::Int
    accepted_steps::Int
    correspondences::Int
    path_time::Float64
    max_path_time::Float64
end

function BenchmarkMetrics(mode::String, trial::Int)
    return BenchmarkMetrics(mode, trial, 0.0, 0, 0, 0, 0, 0.0, 0, 0, 0, 0, 0.0, 0.0)
end

function patch_template_for(nvars::Int)
    raw = ComplexF64[cis(0.37 * k) for k in 0:nvars-1]
    scale = sqrt(sum(abs2, raw))
    return ComplexF64[a / scale for a in raw]
end

function patch_vector_for(CC, nvars::Int)
    return [CC(real(a), imag(a)) for a in patch_template_for(nvars)]
end

mutable struct HighLevelMetrics
    mode::String
    trial::Int
    total_time::Float64
    edges::Int
    correspondences::Int
    complete_edges::Int
end

function compile_benchmark_system(F, vars, pars; projective::Bool)
    patch_vector = projective ? patch_template_for(length(vars) + 1) : nothing
    return redirect_stdout(devnull) do
        compile_edge_homotopy(F, vars, pars; projective=projective, patch_vector=patch_vector)
    end
end

function run_highlevel_trial(mode::String, trial::Int, options)
    F, vars, pars = benchmark_system()
    CC = AcbField(options.precision_bits)
    vertices = initial_vertices(CC, pars, options.num_vertices, options.seed + UInt64(trial - 1))
    projective = mode == "projective"
    compiled = compile_benchmark_system(F, vars, pars; projective=projective)

    elapsed = @elapsed edges = solve_monodromy(compiled, vertices; max_roots=options.max_roots, root_match=options.root_match)
    correspondences = sum(length(e.correspondence12) for e in edges)
    complete_edges = count(e -> length(e.correspondence12) >= options.max_roots, edges)
    return HighLevelMetrics(mode, trial, elapsed, length(edges), correspondences, complete_edges)
end

function track_edge_diagnostic!(
    compiled_sys,
    e::Edge,
    from1to2::Bool,
    metrics::BenchmarkMetrics;
    patch_vector=AcbFieldElem[],
    h_init=0.1,
    max_paths=0,
    root_match::Symbol=:heuristic,
)
    max_paths > 0 && metrics.path_attempts >= max_paths && return

    if from1to2
        source_v, target_v = e.node1, e.node2
        c_forward, c_backward = e.correspondence12, e.correspondence21
    else
        source_v, target_v = e.node2, e.node1
        c_forward, c_backward = e.correspondence21, e.correspondence12
    end

    tracked_indices = Set(x[1] for x in c_forward)
    untracked_indices = [i for i in 1:length(source_v.sols) if !(i in tracked_indices)]
    isempty(untracked_indices) && return

    metrics.edge_calls += 1
    sys_edge = make_edge_system(compiled_sys, source_v.base_point, target_v.base_point; patch_vector=patch_vector)

    for src_idx in untracked_indices
        max_paths > 0 && metrics.path_attempts >= max_paths && return

        metrics.path_attempts += 1
        elapsed = @elapsed result = track_path(sys_edge, source_v.sols[src_idx]; t_end=1.0, h_init=h_init)
        metrics.path_time += elapsed
        metrics.max_path_time = max(metrics.max_path_time, elapsed)
        metrics.iterations += result.iterations
        metrics.accepted_steps += result.accepted_steps
        metrics.rejected_steps += result.rejected_steps

        if succeeded(result)
            metrics.successes += 1
            tracked_root = certified_region(result)
            if !isempty(projective_solution(result)) && metrics.mode == "projective"
                metrics.max_coord_norm = max(metrics.max_coord_norm, vector_norm(projective_solution(result)))
            else
                metrics.max_coord_norm = max(metrics.max_coord_norm, vector_norm(tracked_root))
            end

            dest_idx = if root_match == :heuristic
                search_point(tracked_root, target_v.sols)
            elseif root_match == :certified
                search_point_certified(sys_edge, tracked_root, target_v.sols)
            elseif root_match == :certified_or_heuristic
                idx = search_point_certified(sys_edge, tracked_root, target_v.sols)
                idx === nothing ? search_point(tracked_root, target_v.sols) : idx
            else
                error("Unknown root_match=$root_match")
            end
            if dest_idx === nothing
                push!(target_v.sols, tracked_root)
                dest_idx = length(target_v.sols)
                push!(c_forward, (src_idx, dest_idx))
                push!(c_backward, (dest_idx, src_idx))
            elseif !any(p -> p[2] == dest_idx, c_forward)
                push!(c_forward, (src_idx, dest_idx))
                push!(c_backward, (dest_idx, src_idx))
            end
        else
            metrics.failures += 1
            if !isempty(projective_solution(result)) && metrics.mode == "projective"
                metrics.max_coord_norm = max(metrics.max_coord_norm, vector_norm(projective_solution(result)))
            else
                metrics.max_coord_norm = max(metrics.max_coord_norm, vector_norm(certified_region(result)))
            end
        end
    end

    sort!(c_forward)
    sort!(c_backward)
end

function solve_monodromy_diagnostic!(
    compiled_sys,
    vertices,
    edges,
    metrics::BenchmarkMetrics;
    max_roots=8,
    patch_vector=AcbFieldElem[],
    h_init=0.1,
    max_paths=0,
    root_match::Symbol=:heuristic,
)
    iter_stagnant = 0
    total_correspondences = 0

    while true
        current_correspondences = sum(length(e.correspondence12) for e in edges)
        if current_correspondences == total_correspondences
            iter_stagnant += 1
        else
            iter_stagnant = 0
            total_correspondences = current_correspondences
        end

        if all(e -> length(e.correspondence12) == max_roots, edges)
            break
        end
        max_paths > 0 && metrics.path_attempts >= max_paths && break
        iter_stagnant > 10 && break

        for e in edges
            track_edge_diagnostic!(compiled_sys, e, true, metrics; patch_vector=patch_vector, h_init=h_init, max_paths=max_paths, root_match=root_match)
            max_paths > 0 && metrics.path_attempts >= max_paths && break
            track_edge_diagnostic!(compiled_sys, e, false, metrics; patch_vector=patch_vector, h_init=h_init, max_paths=max_paths, root_match=root_match)
            max_paths > 0 && metrics.path_attempts >= max_paths && break
        end
    end

    metrics.correspondences = sum(length(e.correspondence12) for e in edges)
    return edges
end

function run_trial(mode::String, trial::Int, options)
    F, vars, pars = benchmark_system()
    CC = AcbField(options.precision_bits)
    vertices = initial_vertices(CC, pars, options.num_vertices, options.seed + UInt64(trial - 1))
    edges = complete_graph!(vertices)
    metrics = BenchmarkMetrics(mode, trial)

    projective = mode == "projective"
    compiled = compile_benchmark_system(F, vars, pars; projective=projective)
    patch_vector = projective ? patch_vector_for(CC, length(vars) + 1) : AcbFieldElem[]

    metrics.total_time = @elapsed solve_monodromy_diagnostic!(
        compiled,
        vertices,
        edges,
        metrics;
        max_roots=options.max_roots,
        patch_vector=patch_vector,
        h_init=options.h_init,
        max_paths=options.max_paths,
        root_match=options.root_match,
    )
    return metrics
end

function print_highlevel_table(rows)
    header = (
        "mode", "trial", "total_time", "edges", "correspondences", "complete_edges",
        "failures", "max_coord_norm", "iters", "rejects"
    )
    println(join(header, "  "))
    for r in rows
        values = (
            r.mode,
            string(r.trial),
            string(round(r.total_time; digits=3)),
            string(r.edges),
            string(r.correspondences),
            string(r.complete_edges),
            "unavailable",
            "unavailable",
            "unavailable",
            "unavailable",
        )
        println(join(values, "  "))
    end
end

function print_diagnostic_table(rows)
    header = (
        "mode", "trial", "total_time", "edge_calls", "paths",
        "success", "failure", "max_coord_norm", "iters", "rejects",
        "avg_path_time", "max_path_time", "corr"
    )
    println(join(header, "  "))
    for r in rows
        values = (
            r.mode,
            string(r.trial),
            string(round(r.total_time; digits=3)),
            string(r.edge_calls),
            string(r.path_attempts),
            string(r.successes),
            string(r.failures),
            string(round(r.max_coord_norm; sigdigits=5)),
            string(r.iterations),
            string(r.rejected_steps),
            string(round(r.path_attempts == 0 ? 0.0 : r.path_time / r.path_attempts; digits=3)),
            string(round(r.max_path_time; digits=3)),
            string(r.correspondences),
        )
        println(join(values, "  "))
    end
end

function main(args=ARGS)
    options = parse_options(args)
    modes = options.mode == "both" ? ["affine", "projective"] : [options.mode]
    run_highlevel = options.benchmark in ("highlevel", "both")
    run_diagnostic = options.benchmark in ("diagnostic", "both")

    println("Projective tracking benchmark")
    println("precision=$(options.precision_bits), max_roots=$(options.max_roots), seed=$(options.seed), trials=$(options.trials), vertices=$(options.num_vertices), h_init=$(options.h_init), max_paths=$(options.max_paths), benchmark=$(options.benchmark), root_match=$(options.root_match)")
    println("Projective mode is enabled when compiling the homotopy: compile_edge_homotopy(...; projective=true).")
    println()

    if run_highlevel
        println("High-level solve_monodromy benchmark")
        println("This measures the user-facing call: edges = solve_monodromy(compiled_homotopy, vertices; max_roots=max_roots).")
        println("Per-path failures, coordinate norms, iterations, and rejections are unavailable from current API in this mode.")
        highlevel_rows = HighLevelMetrics[]
        for trial in 1:options.trials
            for mode in modes
                push!(highlevel_rows, run_highlevel_trial(mode, trial, options))
            end
        end
        print_highlevel_table(highlevel_rows)
        println()
    end

    if run_diagnostic
        println("Instrumented diagnostic benchmark")
        println("This uses a script-local diagnostic loop around track_path and does not claim to be identical to solve_monodromy.")
        println("max_coord_norm is measured on returned roots; path-internal coordinate maxima are unavailable from current API.")
        diagnostic_rows = BenchmarkMetrics[]
        for trial in 1:options.trials
            for mode in modes
                push!(diagnostic_rows, run_trial(mode, trial, options))
            end
        end
        print_diagnostic_table(diagnostic_rows)
    end

    return nothing
end

main()
