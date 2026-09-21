include(joinpath(@__DIR__, "..", "common", "higher_order_predictor.jl"))
include(joinpath(@__DIR__, "..", "common", "benchmark_systems.jl"))
include(joinpath(@__DIR__, "..", "common", "csv_io.jl"))

function fixed_cubic_problem(n::Int; seed=20240626)
    spec = BenchmarkSpec(:random_dense, n, 3, seed + n)
    H, roots = make_benchmark(spec; gamma_mode=:diagonal_start)
    return H, first(roots)
end

function run_tradeoff(; dimensions=(3,), max_order=6, max_iterations=50_000)
    rows = NamedTuple[]
    for n in dimensions
        H, point = fixed_cubic_problem(n)
        for order in 0:max_order
            result = Ref{Any}()
            elapsed = @elapsed begin
                if order == 0
                    _, stats = tracking_constant_apriori(
                        H, copy(point), 0.1;
                        iterations_count=true, show_display=false, final_refine=false,
                        max_iterations,
                    )
                    result[] = (status=stats.status, iterations=stats.iterations)
                else
                    _, stats = tracking_higher_order_apriori(
                        H, copy(point), 0.1;
                        order, iterations_count=true, show_display=false,
                        final_refine=false, max_iterations,
                    )
                    result[] = (status=stats.status, iterations=stats.iterations)
                end
            end
            push!(rows, (
                variables=n,
                bezout_degree=3^n,
                predictor_order=order,
                status=String(result[].status),
                iterations=result[].iterations,
                elapsed_sec=elapsed,
                seed=20240626 + n,
            ))
            println("3^$n order $order: $(result[].iterations) iterations, $(round(elapsed; digits=2)) s")
        end
    end
    output = joinpath(@__DIR__, "results", "generated", "order_tradeoff.csv")
    write_csv(output, rows)
    println("wrote $output")
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    full = get(ENV, "CHT_FULL", "0") == "1"
    dimensions = full ? (3, 4, 5, 6) : (3,)
    run_tradeoff(; dimensions)
end
