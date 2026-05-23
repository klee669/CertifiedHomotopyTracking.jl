using Pkg

const PACKAGE_ROOT = normpath(joinpath(@__DIR__, ".."))
PACKAGE_ROOT in LOAD_PATH || push!(LOAD_PATH, PACKAGE_ROOT)

using BenchmarkTools
using CertifiedHomotopyTracking
using Printf

const DEFAULT_CASES = ("readme_straight_line", "small_dense", "monodromy_edge", "surface_optimization")
const CASE_FILTER = Set(
    strip(name) for name in split(get(ENV, "CHT_BENCH_CASES", join(DEFAULT_CASES, ",")), ",")
    if !isempty(strip(name))
)
const SAMPLES = parse(Int, get(ENV, "CHT_BENCH_SAMPLES", "3"))
const SECONDS = parse(Float64, get(ENV, "CHT_BENCH_SECONDS", "10"))
const PREC_BITS = parse(Int, get(ENV, "CHT_BENCH_PREC_BITS", "256"))

BenchmarkTools.DEFAULT_PARAMETERS.samples = SAMPLES
BenchmarkTools.DEFAULT_PARAMETERS.seconds = SECONDS
BenchmarkTools.DEFAULT_PARAMETERS.evals = 1

struct TrackPathBenchmark
    name::String
    description::String
    system::HCSystem
    start::Vector{AcbFieldElem}
    kwargs::NamedTuple
end

function run_track_path(case::TrackPathBenchmark)
    return track_path(case.system, copy(case.start); case.kwargs...)
end

function readme_straight_line_case()
    @variables x y
    CC = AcbField(PREC_BITS)

    F = [x^2 + 3*y - 4, y^2 + 3]
    G = [x^2 - 1, y^2 - 1]
    H = straight_line_homotopy(F, G, [x, y]; CCRing=CC, gamma=1 + im)
    start = [CC(1), CC(-1)]

    return TrackPathBenchmark(
        "readme_straight_line",
        "README straight-line homotopy example",
        H,
        start,
        (; h_init=0.05, show_progress=false),
    )
end

function small_dense_case()
    @variables x y z
    CC = AcbField(PREC_BITS)

    F = [
        x^2 - 1 + 0.11*x*y - 0.07*x*z + 0.05*y*z + 0.03*x - 0.02*y,
        y^2 - 1 - 0.09*x*y + 0.08*x*z + 0.04*y*z + 0.02*z,
        z^2 - 1 + 0.06*x*y - 0.05*x*z + 0.07*y*z - 0.03*x,
    ]
    G = [x^2 - 1, y^2 - 1, z^2 - 1]
    H = straight_line_homotopy(F, G, [x, y, z]; CCRing=CC, gamma=1 + im)
    start = [CC(1), CC(-1), CC(1)]

    return TrackPathBenchmark(
        "small_dense",
        "Small dense quadratic polynomial system",
        H,
        start,
        (; h_init=0.05, show_progress=false),
    )
end

function monodromy_edge_case()
    @variables x y
    @variables p q
    CC = AcbField(PREC_BITS)

    F = [
        p*x^2 + 3*y - 4,
        y^2 + q,
    ]
    compiled = compile_edge_homotopy(F, [x, y], [p, q]; homogeneous=false)
    p_start = [CC(1), CC(-1)]
    p_end = [CC(cis(0.31)), CC(cis(1.27))]
    system = make_edge_system(compiled, p_start, p_end)
    start = [CC(1), CC(1)]

    return TrackPathBenchmark(
        "monodromy_edge",
        "One deterministic monodromy edge from the README-style parameter system",
        system,
        start,
        (; t_end=1.0, h_init=0.1, show_progress=false),
    )
end

function surface_optimization_case()
    @variables x y z λ
    @variables u1 u2 u3
    CC = AcbField(PREC_BITS)

    F = [
        2*(x - u1) - 6*λ*(x^2 + y^2)^2*x,
        2*(y - u2) - 6*λ*(x^2 + y^2)^2*y,
        2*(z - u3) + 4*λ*z^3,
        0*u1 + z^4 - (x^2 + y^2)^3,
    ]
    vars = [x, y, z, λ]
    pars = [u1, u2, u3]
    compiled = compile_edge_homotopy(F, vars, pars; homogeneous=false)

    base = [
        CC(0.09868, 0.675389),
        CC(0.423238, 0.713082),
        CC(0.592351, 0.144969),
    ]
    target = base .+ [
        CC(0.02, -0.01),
        CC(-0.015, 0.012),
        CC(0.01, 0.018),
    ]
    start = [
        CC(1.23836, -0.422501),
        CC(1.19574, -1.0474),
        CC(2.08916, 1.85256),
        CC(-0.0126777, 0.0505892),
    ]
    system = make_edge_system(compiled, base, target)

    return TrackPathBenchmark(
        "surface_optimization",
        "One deterministic edge from examples/irreducible_surface_optimization/surface_optimization.jl",
        system,
        start,
        (; t_end=1.0, h_init=0.1, show_progress=false),
    )
end

function all_cases()
    return [
        readme_straight_line_case(),
        small_dense_case(),
        monodromy_edge_case(),
        surface_optimization_case(),
    ]
end

function selected_cases()
    builders = Dict(
        "readme_straight_line" => readme_straight_line_case,
        "small_dense" => small_dense_case,
        "monodromy_edge" => monodromy_edge_case,
        "surface_optimization" => surface_optimization_case,
    )
    unknown = setdiff(CASE_FILTER, Set(keys(builders)))
    isempty(unknown) || error("Unknown benchmark case(s): $(join(sort!(collect(unknown)), ","))")

    selected = [name for name in DEFAULT_CASES if name in CASE_FILTER]
    isempty(selected) && error("No benchmark cases selected. Set CHT_BENCH_CASES to one or more of $(join(DEFAULT_CASES, ","))")
    return [builders[name]() for name in selected]
end

function result_stats(result)
    return (
        success=succeeded(result),
        status=string(result.status),
        iterations=result.iterations,
        accepted_steps=result.accepted_steps,
        rejected_steps=result.rejected_steps,
        validation_attempts=result.accepted_steps + result.rejected_steps,
        refinement_attempts="not_accessible",
        final_radius=result.final_radius,
        final_krawczyk_norm=result.final_krawczyk_norm,
    )
end

function trial_stats(trial)
    med = median(trial)
    best = minimum(trial)
    return (
        median_time_ns=med.time,
        min_time_ns=best.time,
        memory_bytes=med.memory,
        allocations=med.allocs,
    )
end

function print_header()
    println(join((
        "case",
        "success",
        "status",
        "median_time_ns",
        "min_time_ns",
        "memory_bytes",
        "allocations",
        "iterations",
        "accepted_steps",
        "rejected_steps",
        "validation_attempts",
        "refinement_attempts",
        "final_radius",
        "final_krawczyk_norm",
    ), "\t"))
end

function print_row(case::TrackPathBenchmark, rstats, bstats)
    println(join((
        case.name,
        string(rstats.success),
        rstats.status,
        string(round(Int, bstats.median_time_ns)),
        string(round(Int, bstats.min_time_ns)),
        string(bstats.memory_bytes),
        string(bstats.allocations),
        string(rstats.iterations),
        string(rstats.accepted_steps),
        string(rstats.rejected_steps),
        string(rstats.validation_attempts),
        string(rstats.refinement_attempts),
        @sprintf("%.6e", rstats.final_radius),
        @sprintf("%.6e", rstats.final_krawczyk_norm),
    ), "\t"))
end

function run_benchmarks()
    println("# CertifiedHomotopyTracking track_path benchmarks")
    println("# samples=$SAMPLES seconds=$SECONDS evals=1 precision_bits=$PREC_BITS")
    println("# cases=$(join(sort!(collect(CASE_FILTER)), ","))")
    print_header()

    for case in selected_cases()
        result = run_track_path(case)
        rstats = result_stats(result)
        trial = @benchmark run_track_path($case) evals=1
        bstats = trial_stats(trial)
        print_row(case, rstats, bstats)
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_benchmarks()
end
