include("../complexity_higher_order_predictor.jl")

using Random
using Statistics
using Printf

const CCi = AcbField()
const ORDER = 3
const INITIAL_RADIUS = 0.1
const RESULT_DIR = joinpath(@__DIR__, "results_order3")
const DEFAULT_GAMMA_MODE = Symbol(get(ENV, "CHT_GAMMA_MODE", "diagonal_start"))
const DEFAULT_MAX_ITERATIONS = parse(Int, get(ENV, "CHT_MAX_ITERATIONS", string(typemax(Int))))
const DEFAULT_SEED_OFFSET = get(ENV, "CHT_FRESH_SEED", "0") == "1" ?
    Int(mod(time_ns(), 1_000_000_000)) :
    parse(Int, get(ENV, "CHT_SEED_OFFSET", "0"))

struct BenchmarkSpec
    family::Symbol
    n::Int
    d::Int
    seed::Int
end

function system_name(spec::BenchmarkSpec)
    if spec.family == :katsura
        return "katsura$(spec.n)"
    elseif spec.family == :random_dense
        return "random(2^$(spec.n))"
    end
    error("Unknown benchmark family $(spec.family)")
end

actual_seed(spec::BenchmarkSpec) = spec.seed + DEFAULT_SEED_OFFSET

function create_polynomial_ring(CC, n)
    var_names = ["x_$i" for i in 1:n]
    push!(var_names, "η")
    return CC[var_names...]
end

function partitions_of_degree(n::Int, k::Int)
    if k == 1
        return [[n]]
    elseif n == 0
        return [zeros(Int, k)]
    end

    result = Vector{Vector{Int}}()
    for i in 0:n
        for p in partitions_of_degree(n - i, k - 1)
            push!(result, vcat(i, p))
        end
    end
    return result
end

function generate_random_dense_system(x_vars, n::Int, d::Int)
    R = parent(x_vars[1])
    system = elem_type(typeof(R))[]
    for _ in 1:n
        poly = zero(R)
        for total_deg in 0:d
            for part in partitions_of_degree(total_deg, n)
                term = R(CCi(abs(randn())))
                for (var_idx, power) in enumerate(part)
                    power == 0 && continue
                    term *= x_vars[var_idx]^power
                end
                poly += term
            end
        end
        push!(system, poly - 1)
    end
    return system
end

function bezout_roots(CC, n, d)
    single_roots = [
        CC(cos(2 * k * pi / d) + im * sin(2 * k * pi / d))
        for k in 0:d - 1
    ]
    return [[x...] for x in Iterators.product(fill(single_roots, n)...)]
end

function katsura_start_system(x_vars, n::Int)
    start_polys = elem_type(typeof(parent(x_vars[1])))[x_vars[1] - 1]
    for i in 2:n
        push!(start_polys, x_vars[i]^2 - 1)
    end
    return start_polys
end

function katsura_start_roots(CC, n)
    solutions = Vector{Vector{AcbFieldElem}}()
    for bits in 0:2^(n - 1) - 1
        sol = Vector{AcbFieldElem}(undef, n)
        sol[1] = CC(1)
        for i in 2:n
            sol[i] = ((bits >> (i - 2)) & 1 == 1) ? CC(1) : CC(-1)
        end
        push!(solutions, sol)
    end
    return solutions
end

function katsura_system(x_vars, n::Int)
    R = parent(x_vars[1])
    m = n - 1
    u(i) = abs(i) <= m ? x_vars[abs(i) + 1] : zero(R)

    F = elem_type(typeof(R))[
        -one(R) + sum(u(i) for i in -m:m; init=zero(R)),
    ]

    for i in 0:m - 1
        push!(
            F,
            -u(i) + sum(u(j) * u(i - j) for j in -m:m; init=zero(R)),
        )
    end
    return F
end

function random_gamma()
    return CCi(rand(ComplexF64))
end

function linear_homotopy_with_gamma(G, F, t; gamma_mode=DEFAULT_GAMMA_MODE)
    if gamma_mode == :none
        return [(1 - t) * G[i] + t * F[i] for i in eachindex(G)]
    elseif gamma_mode == :scalar_start
        gamma = random_gamma()
        return [(1 - t) * gamma * G[i] + t * F[i] for i in eachindex(G)]
    elseif gamma_mode == :scalar_target
        gamma = random_gamma()
        return [(1 - t) * G[i] + t * gamma * F[i] for i in eachindex(G)]
    elseif gamma_mode == :diagonal_start
        return [
            (1 - t) * random_gamma() * G[i] + t * F[i]
            for i in eachindex(G)
        ]
    end
    error("Unknown gamma mode $gamma_mode. Use none, scalar_start, scalar_target, or diagonal_start.")
end

function make_benchmark(spec::BenchmarkSpec; gamma_mode=DEFAULT_GAMMA_MODE)
    Random.seed!(actual_seed(spec))
    R, vars = create_polynomial_ring(CCi, spec.n)
    x_vars = vars[1:spec.n]
    HR, (t) = R["t"]

    if spec.family == :katsura
        G = katsura_start_system(x_vars, spec.n)
        F = katsura_system(x_vars, spec.n)
        roots = katsura_start_roots(CCi, spec.n)
    elseif spec.family == :random_dense
        G = [x_vars[i]^spec.d - 1 for i in 1:spec.n]
        F = generate_random_dense_system(x_vars, spec.n, spec.d)
        roots = vec(bezout_roots(CCi, spec.n, spec.d))
    else
        error("Unknown benchmark family $(spec.family)")
    end

    H = linear_homotopy_with_gamma(G, F, t; gamma_mode)
    return H, roots
end

function csv_escape(x)
    s = string(x)
    if occursin(',', s) || occursin('"', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

function write_csv(path, header, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join((csv_escape(getfield(row, Symbol(h))) for h in header), ","))
        end
    end
end

function append_csv_row(path, header, row)
    mkpath(dirname(path))
    needs_header = !isfile(path) || filesize(path) == 0
    open(path, "a") do io
        needs_header && println(io, join(header, ","))
        println(io, join((csv_escape(getfield(row, Symbol(h))) for h in header), ","))
        flush(io)
    end
end

function table_5_1_path_records(specs; max_paths=nothing, start_path=1, gamma_mode=DEFAULT_GAMMA_MODE)
    rows = NamedTuple[]
    for spec in specs
        name = system_name(spec)
        H, roots = make_benchmark(spec; gamma_mode)
        path_count = isnothing(max_paths) ? length(roots) : min(max_paths, length(roots))
        println("5.1 order-$ORDER: $name, paths=$start_path:$path_count/$(length(roots)), gamma_mode=$gamma_mode, seed=$(actual_seed(spec))")
        flush(stdout)

        for path_index in start_path:path_count
            println("  starting path $path_index")
            flush(stdout)
            stats_ref = Ref{Any}()
            elapsed_sec = @elapsed silence_stdout() do
                tracking_higher_order_apriori(
                    H,
                    copy(roots[path_index]),
                    INITIAL_RADIUS;
                    order=ORDER,
                    iterations_count=true,
                    show_display=false,
                    final_refine=false,
                    max_iterations=DEFAULT_MAX_ITERATIONS,
                    compare_constant=false,
                ) |> result -> (stats_ref[] = result[2])
            end
            stats = stats_ref[]
            row = (
                system=name,
                family=String(spec.family),
                n=spec.n,
                seed=actual_seed(spec),
                paths_total=length(roots),
                path_index=path_index,
                status=String(stats.status),
                final_t=stats.final_t,
                elapsed_sec=elapsed_sec,
                iterations=stats.iterations,
                min_dt=stats.min_dt,
                median_dt=stats.median_dt,
                max_dt=stats.max_dt,
                mean_dt=stats.mean_dt,
                min_radius=stats.min_radius,
                median_radius=stats.median_radius,
                mean_radius_theory_gap=stats.mean_radius_theory_gap,
                max_eta=stats.max_eta,
                min_r_theory=stats.min_r_theory,
            )
            push!(rows, row)
            append_csv_row(
                joinpath(RESULT_DIR, "table_5_1_order3_paths.partial.csv"),
                collect(keys(row)),
                row,
            )
            min_dt_str = @sprintf("%.3e", stats.min_dt)
            median_dt_str = @sprintf("%.3e", stats.median_dt)
            elapsed_str = @sprintf("%.2f", elapsed_sec)
            status_str = stats.status == :success ? "" : ", status=$(stats.status), final_t=$(stats.final_t)"
            println("  path $path_index: iters=$(stats.iterations), min_dt=$min_dt_str, median_dt=$median_dt_str, elapsed=$(elapsed_str)s$status_str")
            flush(stdout)
        end
    end
    return rows
end

function summarize_5_1(path_rows)
    systems = unique(row.system for row in path_rows)
    rows = NamedTuple[]
    for system in systems
        rs = filter(row -> row.system == system, path_rows)
        push!(rows, (
            system=system,
            paths=length(rs),
            total_elapsed_sec=sum(row.elapsed_sec for row in rs),
            avg_iterations=mean(row.iterations for row in rs),
            avg_min_dt=mean(row.min_dt for row in rs),
            avg_median_dt=mean(row.median_dt for row in rs),
            avg_min_radius=mean(row.min_radius for row in rs),
            avg_radius_theory_gap=mean(row.mean_radius_theory_gap for row in rs),
            max_eta=maximum(row.max_eta for row in rs),
        ))
    end
    return rows
end

function silence_stdout(f)
    return redirect_stdout(devnull) do
        f()
    end
end

function track_iteration_count(H, point, method)
    if method == :adaptive_constant
        _, iters = silence_stdout() do
            CertifiedHomotopyTracking.tracking_without_predictor(
                H,
                copy(point);
                r=INITIAL_RADIUS,
                iterations_count=true,
            )
        end
        return (iterations=iters, status=:success, final_t=1.0)
    elseif method == :apriori_constant
        _, stats = silence_stdout() do
            tracking_constant_apriori(
                H,
                copy(point),
                INITIAL_RADIUS;
                iterations_count=true,
                show_display=false,
                final_refine=false,
            )
        end
        return (iterations=stats.iterations, status=:success, final_t=1.0)
    elseif method == :truncated_hermite
        _, iters = silence_stdout() do
            CertifiedHomotopyTracking.track(
                H,
                copy(point);
                r=INITIAL_RADIUS,
                tracking="truncate",
                iterations_count=true,
                show_display=false,
            )
        end
        return (iterations=iters, status=:success, final_t=1.0)
    elseif method == :nontruncated_hermite
        _, iters = silence_stdout() do
            CertifiedHomotopyTracking.track(
                H,
                copy(point);
                r=INITIAL_RADIUS,
                tracking="non-truncate",
                iterations_count=true,
                show_display=false,
            )
        end
        return (iterations=iters, status=:success, final_t=1.0)
    elseif method == :order3_apriori
        _, stats = silence_stdout() do
            tracking_higher_order_apriori(
                H,
                copy(point),
                INITIAL_RADIUS;
                order=ORDER,
                iterations_count=true,
                show_display=false,
                final_refine=false,
                max_iterations=DEFAULT_MAX_ITERATIONS,
                compare_constant=false,
            )
        end
        return (iterations=stats.iterations, status=stats.status, final_t=stats.final_t)
    end
    error("Unknown method $method")
end

const ALL_5_3_METHODS = [
    :adaptive_constant,
    :apriori_constant,
    :truncated_hermite,
    :nontruncated_hermite,
    :order3_apriori,
]

function table_5_3_records(specs; max_paths=nothing, start_path=1, methods=ALL_5_3_METHODS, gamma_mode=DEFAULT_GAMMA_MODE)
    rows = NamedTuple[]

    for spec in specs
        name = system_name(spec)
        H, roots = make_benchmark(spec; gamma_mode)
        path_count = isnothing(max_paths) ? length(roots) : min(max_paths, length(roots))
        println("5.3 comparison: $name, paths=$start_path:$path_count/$(length(roots)), gamma_mode=$gamma_mode, seed=$(actual_seed(spec))")
        flush(stdout)

        for path_index in start_path:path_count
            for method in methods
                result_ref = Ref{Any}()
                elapsed_sec = @elapsed begin
                    result_ref[] = track_iteration_count(H, roots[path_index], method)
                end
                result = result_ref[]
                row = (
                    system=name,
                    family=String(spec.family),
                    n=spec.n,
                    seed=actual_seed(spec),
                    paths_total=length(roots),
                    path_index=path_index,
                    method=String(method),
                    status=String(result.status),
                    final_t=result.final_t,
                    elapsed_sec=elapsed_sec,
                    iterations=result.iterations,
                )
                push!(rows, row)
                append_csv_row(
                    joinpath(RESULT_DIR, "table_5_3_order3_paths.partial.csv"),
                    collect(keys(row)),
                    row,
                )
                elapsed_str = @sprintf("%.2f", elapsed_sec)
                status_str = result.status == :success ? "" : ", status=$(result.status), final_t=$(result.final_t)"
                println("  path $path_index $(method): $(result.iterations), elapsed=$(elapsed_str)s$status_str")
                flush(stdout)
            end
        end
    end

    return rows
end

function summarize_5_3(iter_rows)
    systems = unique(row.system for row in iter_rows)
    methods = unique(row.method for row in iter_rows)
    rows = NamedTuple[]
    for system in systems
        for method in methods
            rs = filter(row -> row.system == system && row.method == method, iter_rows)
            isempty(rs) && continue
            push!(rows, (
                system=system,
                method=method,
                paths=length(rs),
                total_elapsed_sec=sum(row.elapsed_sec for row in rs),
                avg_iterations=mean(row.iterations for row in rs),
            ))
        end
    end
    return rows
end

function print_summary(title, rows)
    println()
    println(title)
    for row in rows
        println(row)
    end
end

function benchmark_specs(; full=false, only_family=nothing, only_n=nothing)
    ns = full ? (3:6) : (3:3)
    specs = BenchmarkSpec[]
    for n in ns
        push!(specs, BenchmarkSpec(:katsura, n, 2, 10_000 + n))
    end
    for n in ns
        push!(specs, BenchmarkSpec(:random_dense, n, 2, 20_000 + n))
    end
    if !isnothing(only_family)
        fam = Symbol(only_family)
        specs = filter(spec -> spec.family == fam, specs)
    end
    if !isnothing(only_n)
        specs = filter(spec -> spec.n == only_n, specs)
    end
    return specs
end

function run_all_tables_3rd_order(;
    full=false,
    max_paths=nothing,
    start_path=1,
    methods=ALL_5_3_METHODS,
    only_family=nothing,
    only_n=nothing,
    gamma_mode=DEFAULT_GAMMA_MODE,
    skip_5_1=false,
    skip_5_3=false,
)
    specs = benchmark_specs(; full, only_family, only_n)
    println("Seed offset: $DEFAULT_SEED_OFFSET")
    path_rows = NamedTuple[]
    summary_51 = NamedTuple[]
    iter_rows = NamedTuple[]
    summary_53 = NamedTuple[]

    if !skip_5_1
        path_rows = table_5_1_path_records(specs; max_paths, start_path, gamma_mode)
        summary_51 = summarize_5_1(path_rows)
        if !isempty(path_rows)
            write_csv(joinpath(RESULT_DIR, "table_5_1_order3_paths.csv"), collect(keys(first(path_rows))), path_rows)
            write_csv(joinpath(RESULT_DIR, "table_5_1_order3_summary.csv"), collect(keys(first(summary_51))), summary_51)
        end
        print_summary("5.1 order-$ORDER summary", summary_51)
    end

    if !skip_5_3
        iter_rows = table_5_3_records(specs; max_paths, start_path, methods, gamma_mode)
        summary_53 = summarize_5_3(iter_rows)
        if !isempty(iter_rows)
            write_csv(joinpath(RESULT_DIR, "table_5_3_order3_paths.csv"), collect(keys(first(iter_rows))), iter_rows)
            write_csv(joinpath(RESULT_DIR, "table_5_3_order3_summary.csv"), collect(keys(first(summary_53))), summary_53)
        end
        print_summary("5.3 order-$ORDER comparison summary", summary_53)
    end
    println()
    println("Wrote CSV files to $RESULT_DIR")

    return (; path_rows, summary_51, iter_rows, summary_53)
end

if abspath(PROGRAM_FILE) == @__FILE__
    full = get(ENV, "CHT_FULL_TABLES", "0") == "1"
    max_paths_env = get(ENV, "CHT_MAX_PATHS", "")
    max_paths = isempty(max_paths_env) ? nothing : parse(Int, max_paths_env)
    start_path = parse(Int, get(ENV, "CHT_START_PATH", "1"))
    all_methods = get(ENV, "CHT_ALL_METHODS", "0") == "1"
    methods = (full && !all_methods) ? [:order3_apriori] : ALL_5_3_METHODS
    only_family_env = get(ENV, "CHT_ONLY_FAMILY", "")
    only_family = isempty(only_family_env) ? nothing : only_family_env
    only_n_env = get(ENV, "CHT_ONLY_N", "")
    only_n = isempty(only_n_env) ? nothing : parse(Int, only_n_env)
    skip_5_1 = get(ENV, "CHT_SKIP_5_1", "0") == "1"
    skip_5_3 = get(ENV, "CHT_SKIP_5_3", "0") == "1"
    run_all_tables_3rd_order(; full, max_paths, start_path, methods, only_family, only_n, skip_5_1, skip_5_3)
end
