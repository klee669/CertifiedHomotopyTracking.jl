include("complexity_higher_order_random_katsura.jl")

using Random

function fixed_random_dense_problem(; n=5, d=3, seed=20240626)
    Random.seed!(seed)

    R, vars = create_polynomial_ring(CCi, n)
    x_vars = vars[1:n]
    HR, (t) = R["t"]

    G = transpose(hcat([x_vars[i]^d - 1 for i in 1:n]))
    F = generate_random_dense_system(x_vars, n, d)
    H = random_diagonal_homotopy(G, F, t)
    point = first(vec(generate_all_roots(CCi, n, d)))

    return H, point
end

function order_tradeoff_random(;
    n=5,
    d=3,
    seed=20240626,
    max_order=6,
    initial_radius=INITIAL_RADIUS,
)
    H, point = fixed_random_dense_problem(; n, d, seed)
    rows = []

    constant_elapsed = @elapsed begin
        _, constant_stats = tracking_constant_apriori(
            H,
            copy(point),
            initial_radius;
            iterations_count=true,
            show_display=true,
        )
    end
    push!(rows, (
        method="constant",
        order=0,
        iterations=constant_stats.iterations,
        mean_dt=constant_stats.mean_dt,
        median_dt=constant_stats.median_dt,
        elapsed=constant_elapsed,
        iter_per_second=constant_stats.iterations / constant_elapsed,
    ))

    track_elapsed = @elapsed begin
        _, track_iterations = track(
            H,
            copy(point);
            r=initial_radius,
            iterations_count=true,
            show_display=true,
        )
    end
    push!(rows, (
        method="track-default",
        order=-1,
        iterations=track_iterations,
        mean_dt=1 / track_iterations,
        median_dt=NaN,
        elapsed=track_elapsed,
        iter_per_second=track_iterations / track_elapsed,
    ))

    for order in 1:max_order
        elapsed = @elapsed begin
            _, stats = tracking_higher_order_apriori(
                H,
                copy(point),
                initial_radius;
                order,
                iterations_count=true,
                show_display=true,
            )
        end
        push!(rows, (
            method="apriori",
            order=order,
            iterations=stats.iterations,
            mean_dt=stats.mean_dt,
            median_dt=stats.median_dt,
            elapsed=elapsed,
            iter_per_second=stats.iterations / elapsed,
        ))
    end

    println("Fixed random dense system: n=$n, d=$d, seed=$seed")
    println("method, order, iterations, iter_saved_from_previous_order, mean_dt, median_dt, elapsed_sec, iter_per_sec")
    previous_apriori_iterations = nothing
    for row in rows
        iter_saved = ""
        if row.method == "apriori"
            iter_saved = isnothing(previous_apriori_iterations) ? "" : string(previous_apriori_iterations - row.iterations)
            previous_apriori_iterations = row.iterations
        end
        println(join((
            row.method,
            row.order,
            row.iterations,
            iter_saved,
            row.mean_dt,
            row.median_dt,
            round(row.elapsed; digits=4),
            round(row.iter_per_second; digits=4),
        ), ", "))
    end

    apriori_rows = filter(row -> row.method == "apriori", rows)
    best_iteration_row = argmin(row -> row.iterations, apriori_rows)
    best_elapsed_row = argmin(row -> row.elapsed, apriori_rows)
    println("best apriori iterations: order=$(best_iteration_row.order), iterations=$(best_iteration_row.iterations)")
    println("best apriori elapsed: order=$(best_elapsed_row.order), elapsed=$(round(best_elapsed_row.elapsed; digits=4)) sec")

    return rows
end

rows = order_tradeoff_random(max_order=8)
