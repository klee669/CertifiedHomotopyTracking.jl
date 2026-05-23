using CertifiedHomotopyTracking

# Compare fixed 256-bit tracking against the adaptive precision wrapper on the
# 27-lines symmetric surface loops.
#
# Run:
#   julia --project=. examples/27_lines_on_symmetric_surface/27_lines_benchmark.jl
#
# Optional environment variables:
#   CHT_27_MAX_ROOTS=1
#   CHT_27_LOOP=red        # red, green, both
#   CHT_27_PRECISIONS=96,128,192,256
#   CHT_27_SHOW_PROGRESS=false

function _env_int(name, default)
    return parse(Int, get(ENV, name, string(default)))
end

function _env_bool(name, default)
    raw = lowercase(get(ENV, name, default ? "true" : "false"))
    return raw in ("1", "true", "yes", "on")
end

function _env_precisions(name, default)
    raw = get(ENV, name, join(default, ","))
    return Tuple(parse(Int, strip(x)) for x in split(raw, ",") if !isempty(strip(x)))
end

const MAX_ROOTS = _env_int("CHT_27_MAX_ROOTS", 27)
const LOOP_MODE = lowercase(get(ENV, "CHT_27_LOOP", "red"))
const PRECISIONS = _env_precisions("CHT_27_PRECISIONS", (96, 128, 192, 256))
const SHOW_PROGRESS = _env_bool("CHT_27_SHOW_PROGRESS", true)

function _warmup_trackers()
    @variables x
    CC = AcbField(128)
    H = straight_line_homotopy([x - 2], [x - 1], [x]; CCRing=CC)
    start = [CC(1)]
    track_path(H, start; h_init=0.05)
    track_path_adaptive_precision(H, start; precisions=(96, 128), h_init=0.05)
    return nothing
end

_warmup_trackers()

ENV["CHT_27_LINES_SETUP_ONLY"] = "1"
include("27_lines_symmetric.jl")
delete!(ENV, "CHT_27_LINES_SETUP_ONLY")

function _convert_root(sys, root)
    return AcbFieldElem[sys.CC(z) for z in root]
end

function _precision_of(res)
    root = certified_region(res)
    isempty(root) && return 0
    return precision(parent(root[1]))
end

function _final_residual(sys, res)
    root = _convert_root(sys, certified_region(res))
    return norm_inf(evaluate_H(sys, root, sys.CC(1)))
end

function _run_path(mode, sys, start)
    runner = mode == :fixed256 ?
        () -> track_path(sys, start; t_end=1.0, h_init=0.1, show_progress=SHOW_PROGRESS) :
        () -> track_path_adaptive_precision(sys, start; precisions=PRECISIONS, t_end=1.0, h_init=0.1, show_progress=SHOW_PROGRESS)

    GC.gc()
    timed = @timed runner()
    res = timed.value
    return (
        result = res,
        time = timed.time,
        bytes = timed.bytes,
        residual = try
            _final_residual(sys, res)
        catch
            Inf
        end,
    )
end

function _combine(path_stats)
    return (
        paths = length(path_stats),
        success = count(s -> succeeded(s.result), path_stats),
        accepted = sum(s.result.accepted_steps for s in path_stats),
        rejected = sum(s.result.rejected_steps for s in path_stats),
        iterations = sum(s.result.iterations for s in path_stats),
        time = sum(s.time for s in path_stats),
        bytes = sum(s.bytes for s in path_stats),
        max_residual = maximum(s.residual for s in path_stats),
        max_krawczyk = maximum(s.result.final_krawczyk_norm for s in path_stats),
        max_precision = maximum(_precision_of(s.result) for s in path_stats),
        statuses = join(unique(string(s.result.status) for s in path_stats), ","),
    )
end

function _print_header()
    println(join((
        "loop",
        "root",
        "mode",
        "mapped",
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

function _print_row(loop_name, root_idx, mode, mapped, stats)
    println(join((
        loop_name,
        string(root_idx),
        string(mode),
        string(mapped),
        string(stats.paths),
        string(stats.success),
        string(stats.accepted),
        string(stats.rejected),
        string(stats.iterations),
        string(round(stats.time; sigdigits=6)),
        string(round(stats.bytes / 2.0^20; sigdigits=6)),
        string(round(stats.max_residual; sigdigits=6)),
        string(round(stats.max_krawczyk; sigdigits=6)),
        string(stats.max_precision),
        stats.statuses,
    ), "\t"))
end

function run_triangle_loop(mode, bp, a, b, x0)
    path_stats = []

    F1 = make_edge_system(compiled_homotopy, bp, a)
    s1 = _run_path(mode, F1, x0)
    push!(path_stats, s1)
    succeeded(s1.result) || return path_stats, nothing

    F2 = make_edge_system(compiled_homotopy, a, b)
    s2 = _run_path(mode, F2, _convert_root(F1, certified_region(s1.result)))
    push!(path_stats, s2)
    succeeded(s2.result) || return path_stats, nothing

    F3 = make_edge_system(compiled_homotopy, b, bp)
    s3 = _run_path(mode, F3, _convert_root(F2, certified_region(s2.result)))
    push!(path_stats, s3)
    succeeded(s3.result) || return path_stats, nothing

    x3 = _convert_root(F3, certified_region(s3.result))
    mapped = search_point_certified(F3, x3, p_list)
    return path_stats, mapped
end

function run_loop(loop_name, bp, a, b)
    n_roots = min(MAX_ROOTS, length(p_list))
    for root_idx in 1:n_roots
        for mode in (:fixed256, :adaptive)
            stats_list, mapped = run_triangle_loop(mode, bp, a, b, p_list[root_idx])
            stats = _combine(stats_list)
            _print_row(loop_name, root_idx, mode, mapped, stats)
        end
    end
end

println("27-lines tracking benchmark")
println("max_roots=$MAX_ROOTS, loop=$LOOP_MODE, precisions=$PRECISIONS, show_progress=$SHOW_PROGRESS")
println("fixed256 = track_path at 256 bits; adaptive = track_path_adaptive_precision")
println()
_print_header()

if LOOP_MODE in ("red", "both")
    run_loop("red", red1, red2, red3)
end

if LOOP_MODE in ("green", "both")
    run_loop("green", green1, green2, green3)
end
