include("../complexity_higher_order_predictor.jl")

using Random
using Statistics

const RUNS = 10          # Change this to 100 for the full experiment.
const ORDER = 3         # Higher predictor order q.
const INITIAL_RADIUS = 0.1

CCi = AcbField()

function create_polynomial_ring(CCi, n)
    var_names = ["x_$i" for i in 1:n]
    push!(var_names, "η")
    return CCi[var_names...]
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

function generate_random_dense_system(vars, n::Int, d::Int)
    R = parent(vars[1])
    system = elem_type(typeof(R))[]

    for _ in 1:n
        poly = zero(R)
        for total_deg in 0:d
            for part in partitions_of_degree(total_deg, n)
                term = R(CCi(abs(randn())))
                for (var_idx, power) in enumerate(part)
                    power == 0 && continue
                    term *= vars[var_idx]^power
                end
                poly += term
            end
        end
        push!(system, poly - 1)
    end

    return transpose(hcat(system))
end

function generate_all_roots(CCi, n, d)
    single_roots = [
        CCi(cos(2 * k * pi / d) + im * sin(2 * k * pi / d))
        for k in 0:d - 1
    ]
    return [[x...] for x in Iterators.product(fill(single_roots, n)...)]
end

function katsura_start_system(x_vars, n::Int)
    start_polys = elem_type(typeof(parent(x_vars[1])))[x_vars[1] - 1]
    for i in 2:n
        push!(start_polys, x_vars[i]^2 - 1)
    end
    return transpose(hcat(start_polys))
end

function katsura_start_roots(CCi, n)
    solutions = Vector{Vector{AcbFieldElem}}()
    for bits in 0:2^(n - 1) - 1
        sol = Vector{AcbFieldElem}(undef, n)
        sol[1] = CCi(1)
        for i in 2:n
            sol[i] = ((bits >> (i - 2)) & 1 == 1) ? CCi(1) : CCi(-1)
        end
        push!(solutions, sol)
    end
    return solutions
end

function katsura_system(x_vars, n::Int)
    R = parent(x_vars[1])
    F = elem_type(typeof(R))[]
    for i in 0:n - 1
        poly = zero(R)
        for j in -n:n
            k = i - j
            (0 <= abs(j) <= n - 1 && 0 <= abs(k) <= n - 1) || continue
            xj = x_vars[abs(j) + 1]
            xk = x_vars[abs(k) + 1]
            poly += xj * xk
        end
        push!(F, poly - x_vars[i + 1])
    end
    F[1] += x_vars[1] - 1
    return transpose(hcat(F))
end

function random_diagonal_homotopy(G, F, t)
    n = length(G)
    return [
        (1 - t) * CCi(rand(ComplexF64)) * G[i] + t * F[i]
        for i in 1:n
    ]
end

function run_random_dense_smoke(; n=4, d=2, runs=RUNS, order=ORDER)
    R, vars = create_polynomial_ring(CCi, n)
    x_vars = vars[1:n]
    HR, (t) = R["t"]

    G = transpose(hcat([x_vars[i]^d - 1 for i in 1:n]))
    roots = vec(generate_all_roots(CCi, n, d))
    stats_list = []

    for i in 1:runs
        random_system = generate_random_dense_system(x_vars, n, d)
        H = random_diagonal_homotopy(G, random_system, t)
        point = roots[i]

#        _, constant_stats = tracking_constant_apriori(
#            H,
#            copy(point),
#            INITIAL_RADIUS;
#            iterations_count=true,
#            show_display=false,
#        )
        _, stats = tracking_higher_order_apriori(
            H,
            copy(point),
            INITIAL_RADIUS;
            order,
            iterations_count=true,
            show_display=false,
        )
        _, hermite_iterations = track(
            H,
            copy(point);
            r=INITIAL_RADIUS,
            iterations_count=true,
            show_display=false,
        )
        push!(stats_list, stats)
        # println("random run $i constant: ", constant_stats)
        println("random run $i order-$order: ", stats)
        println("random run $i track default: iterations = ", hermite_iterations)
        # println("random run $i iteration ratio constant/order-$order = ", constant_stats.iterations / stats.iterations)
        println("random run $i iteration ratio track/order-$order = ", hermite_iterations / stats.iterations)
    end

    return stats_list
end

function run_katsura_smoke(; n=4, runs=RUNS, order=ORDER)
    R, vars = create_polynomial_ring(CCi, n)
    x_vars = vars[1:n]
    HR, (t) = R["t"]

    G = katsura_start_system(x_vars, n)
    F = katsura_system(x_vars, n)
    roots = katsura_start_roots(CCi, n)
    stats_list = []

    for i in 1:runs
        H = random_diagonal_homotopy(G, F, t)
        point = roots[i]

#        _, constant_stats = tracking_constant_apriori(
#            H,
#            copy(point),
#            INITIAL_RADIUS;
#            iterations_count=true,
#            show_display=false,
#        )
        _, stats = tracking_higher_order_apriori(
            H,
            copy(point),
            INITIAL_RADIUS;
            order,
            iterations_count=true,
            show_display=false,
        )
        _, hermite_iterations = track(
            H,
            copy(point);
            r=INITIAL_RADIUS,
            iterations_count=true,
            show_display=false,
        )
        push!(stats_list, stats)
#        println("katsura run $i constant: ", constant_stats)
        println("katsura run $i order-$order: ", stats)
        println("katsura run $i track default: iterations = ", hermite_iterations)
#        println("katsura run $i iteration ratio constant/order-$order = ", constant_stats.iterations / stats.iterations)
        println("katsura run $i iteration ratio track/order-$order = ", hermite_iterations / stats.iterations)
    end

    return stats_list
end

if abspath(PROGRAM_FILE) == @__FILE__
    random_stats = run_random_dense_smoke()
    katsura_stats = run_katsura_smoke()
end
