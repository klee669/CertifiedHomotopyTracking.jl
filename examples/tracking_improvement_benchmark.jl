using CertifiedHomotopyTracking
using Statistics

# Compare the current fixed-precision tracker against the opt-in adaptive
# precision wrapper on deterministic examples.
#
# Run:
#   julia --project=. examples/tracking_improvement_benchmark.jl
#
# Optional environment variables:
#   CHT_BENCH_CASES=small,surface
#   CHT_BENCH_REPEATS=1
#   CHT_BENCH_PRECISIONS=96,128,192,256
#   CHT_SURFACE_PROFILE=nearby   # nearby, original
#   CHT_SURFACE_SEGMENTS=1

function _env_int(name, default)
    return parse(Int, get(ENV, name, string(default)))
end

function _env_precisions(name, default)
    raw = get(ENV, name, join(default, ","))
    return Tuple(parse(Int, strip(x)) for x in split(raw, ",") if !isempty(strip(x)))
end

function _env_cases()
    raw = lowercase(get(ENV, "CHT_BENCH_CASES", "small,surface"))
    return Set(strip(x) for x in split(raw, ",") if !isempty(strip(x)))
end

const REPEATS = _env_int("CHT_BENCH_REPEATS", 1)
const PRECISIONS = _env_precisions("CHT_BENCH_PRECISIONS", (96, 128, 192, 256))
const CASES = _env_cases()
const SURFACE_PROFILE = lowercase(get(ENV, "CHT_SURFACE_PROFILE", "nearby"))
const SURFACE_SEGMENTS = _env_int("CHT_SURFACE_SEGMENTS", 1)

function _convert_root(sys, root)
    return AcbFieldElem[sys.CC(z) for z in root]
end

function _final_residual(sys, res)
    root = _convert_root(sys, certified_region(res))
    return norm_inf(evaluate_H(sys, root, sys.CC(1)))
end

function _precision_of(res)
    root = certified_region(res)
    isempty(root) && return 0
    return precision(parent(root[1]))
end

function _path_stats(sys, result, elapsed, bytes)
    residual = try
        _final_residual(sys, result)
    catch
        Inf
    end
    return (
        paths = 1,
        success = succeeded(result) ? 1 : 0,
        accepted = result.accepted_steps,
        rejected = result.rejected_steps,
        iterations = result.iterations,
        time = elapsed,
        bytes = bytes,
        max_residual = residual,
        max_krawczyk = result.final_krawczyk_norm,
        precision = _precision_of(result),
        statuses = string(result.status),
    )
end

function _combine_stats(stats)
    isempty(stats) && return (
        paths = 0,
        success = 0,
        accepted = 0,
        rejected = 0,
        iterations = 0,
        time = 0.0,
        bytes = 0,
        max_residual = Inf,
        max_krawczyk = Inf,
        precision = 0,
        statuses = "",
    )
    return (
        paths = sum(s.paths for s in stats),
        success = sum(s.success for s in stats),
        accepted = sum(s.accepted for s in stats),
        rejected = sum(s.rejected for s in stats),
        iterations = sum(s.iterations for s in stats),
        time = sum(s.time for s in stats),
        bytes = sum(s.bytes for s in stats),
        max_residual = maximum(s.max_residual for s in stats),
        max_krawczyk = maximum(s.max_krawczyk for s in stats),
        precision = maximum(s.precision for s in stats),
        statuses = join(unique(s.statuses for s in stats), ","),
    )
end

function _print_header()
    println(join((
        "case",
        "mode",
        "paths",
        "success",
        "accepted",
        "rejected",
        "iters",
        "time_s",
        "alloc_MB",
        "max_residual",
        "max_krawczyk",
        "max_precision",
        "statuses",
    ), "\t"))
end

function _print_row(case_name, mode, stats)
    println(join((
        case_name,
        mode,
        string(stats.paths),
        string(stats.success),
        string(stats.accepted),
        string(stats.rejected),
        string(stats.iterations),
        string(round(stats.time; sigdigits=6)),
        string(round(stats.bytes / 2.0^20; sigdigits=6)),
        string(round(stats.max_residual; sigdigits=6)),
        string(round(stats.max_krawczyk; sigdigits=6)),
        string(stats.precision),
        stats.statuses,
    ), "\t"))
end

function _run_once(mode, sys, start; kwargs...)
    runner = mode == :fixed256 ?
        () -> track_path(sys, start; kwargs...) :
        () -> track_path_adaptive_precision(sys, start; precisions=PRECISIONS, kwargs...)
    GC.gc()
    timed = @timed runner()
    return timed.value, timed.time, timed.bytes
end

function _warmup_trackers()
    @variables x
    CC = AcbField(128)
    H = straight_line_homotopy([x - 2], [x - 1], [x]; CCRing=CC)
    start = [CC(1)]
    track_path(H, start; h_init=0.05)
    track_path_adaptive_precision(H, start; precisions=(96, 128), h_init=0.05)
    return nothing
end

function _median_run(case_name, mode, path_builder; kwargs...)
    path_builder(mode; kwargs...)
    runs = []
    for _ in 1:REPEATS
        stats = path_builder(mode; kwargs...)
        push!(runs, stats)
    end
    sort!(runs, by = s -> s.time)
    _print_row(case_name, string(mode), runs[cld(length(runs), 2)])
end

function small_quadratic_path(mode; kwargs...)
    @variables x y
    CC = AcbField(256)
    F = [x^2 + 3*y - 4, y^2 + 3]
    G = [x^2 - 1, y^2 - 1]
    H = straight_line_homotopy(F, G, [x, y]; CCRing=CC)
    start = [CC(1), CC(-1)]
    res, elapsed, bytes = _run_once(mode, H, start; h_init=0.05)
    return _path_stats(H, res, elapsed, bytes)
end

function surface_setup()
    @variables x y z λ
    @variables u1 u2 u3
    CC = AcbField(256)

    F = [
        2*(x-u1)-6*λ*(x^2+y^2)^2*x,
        2*(y-u2)-6*λ*(x^2+y^2)^2*y,
        2*(z-u3)+4*λ*z^3,
        0*u1+z^4-(x^2+y^2)^3,
    ]
    vars = [x, y, z, λ]
    pars = [u1, u2, u3]

    base = [CC(.09868,.675389), CC(.423238,.713082), CC(.592351,.144969)]
    x0 = [CC(1.23836,-.422501), CC(1.19574,-1.0474), CC(2.08916,1.85256), CC(-.0126777,.0505892)]
    vertices = if SURFACE_PROFILE == "nearby"
        [
            base,
            base .+ [CC(0.02,-0.01), CC(-0.015,0.012), CC(0.01,0.018)],
            base .+ [CC(0.04,-0.02), CC(-0.03,0.024), CC(0.02,0.036)],
        ]
    elseif SURFACE_PROFILE == "original"
        [
            base,
            [CC(cis(0.31)), CC(cis(1.27)), CC(cis(2.49))],
            [CC(cis(2.11)), CC(cis(0.73)), CC(cis(4.02))],
            [CC(cis(5.62)), CC(cis(3.18)), CC(cis(1.45))],
        ]
    else
        error("Unknown CHT_SURFACE_PROFILE=$SURFACE_PROFILE. Use nearby or original.")
    end
    compiled = compile_edge_homotopy(F, vars, pars)
    return compiled, vertices, x0
end

function surface_paths(mode; kwargs...)
    compiled, vertices, x = surface_setup()
    stats = []
    max_segments = min(SURFACE_SEGMENTS, length(vertices) - 1)

    for i in 1:max_segments
        sys = make_edge_system(compiled, vertices[i], vertices[i + 1])
        res, elapsed, bytes = _run_once(mode, sys, x; t_end=1.0, h_init=0.1, show_progress=false)
        push!(stats, _path_stats(sys, res, elapsed, bytes))
        succeeded(res) || break
        x = _convert_root(sys, certified_region(res))
    end

    return _combine_stats(stats)
end

println("Tracking improvement benchmark")
println("repeats=$REPEATS, precisions=$PRECISIONS, surface_profile=$SURFACE_PROFILE, surface_segments=$SURFACE_SEGMENTS")
println("fixed256 = track_path at the system precision; adaptive = track_path_adaptive_precision")
println()
_warmup_trackers()
_print_header()

if "small" in CASES
    _median_run("small_quadratic", :fixed256, small_quadratic_path)
    _median_run("small_quadratic", :adaptive, small_quadratic_path)
end

if "surface" in CASES
    _median_run("surface_optimization", :fixed256, surface_paths)
    _median_run("surface_optimization", :adaptive, surface_paths)
end
