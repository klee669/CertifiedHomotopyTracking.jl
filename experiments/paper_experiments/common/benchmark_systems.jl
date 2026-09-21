using CertifiedHomotopyTracking
using Random

const CC = AcbField()

struct BenchmarkSpec
    family::Symbol
    n::Int
    degree::Int
    seed::Int
end

function system_name(spec::BenchmarkSpec)
    spec.family == :katsura && return "katsura$(spec.n)"
    spec.family == :random_dense && return "random($(spec.degree)^$(spec.n))"
    error("unknown benchmark family: $(spec.family)")
end

function polynomial_ring(n::Int)
    names = ["x_$i" for i in 1:n]
    push!(names, "eta")
    return CC[names...]
end

function degree_partitions(total::Int, variables::Int)
    variables == 1 && return [[total]]
    total == 0 && return [zeros(Int, variables)]
    result = Vector{Vector{Int}}()
    for exponent in 0:total
        for tail in degree_partitions(total - exponent, variables - 1)
            push!(result, vcat(exponent, tail))
        end
    end
    return result
end

function random_dense_system(x, degree::Int)
    ring = parent(first(x))
    n = length(x)
    equations = elem_type(typeof(ring))[]
    for _ in 1:n
        polynomial = zero(ring)
        for total_degree in 0:degree
            for exponents in degree_partitions(total_degree, n)
                term = ring(CC(abs(randn())))
                for (variable, exponent) in zip(x, exponents)
                    exponent == 0 || (term *= variable^exponent)
                end
                polynomial += term
            end
        end
        push!(equations, polynomial - 1)
    end
    return equations
end

bezout_start_system(x, degree::Int) = [variable^degree - 1 for variable in x]

function bezout_roots(n::Int, degree::Int)
    roots = [CC(cis(2 * pi * k / degree)) for k in 0:degree - 1]
    return [[coordinates...] for coordinates in Iterators.product(fill(roots, n)...)]
end

function katsura_start_system(x)
    return [x[1] - 1; [x[i]^2 - 1 for i in 2:length(x)]]
end

function katsura_start_roots(n::Int)
    return [
        [CC(1); [((bits >> (i - 2)) & 1 == 1) ? CC(1) : CC(-1) for i in 2:n]]
        for bits in 0:2^(n - 1) - 1
    ]
end

# Macaulay2 ExampleSystems convention: u(-i)=u(i), u(i)=0 for i>n-1.
function katsura_system(x)
    ring = parent(first(x))
    m = length(x) - 1
    u(i) = abs(i) <= m ? x[abs(i) + 1] : zero(ring)
    equations = elem_type(typeof(ring))[
        -one(ring) + sum(u(i) for i in -m:m; init=zero(ring)),
    ]
    for i in 0:m - 1
        push!(equations, -u(i) + sum(u(j) * u(i - j) for j in -m:m; init=zero(ring)))
    end
    return equations
end

random_gamma() = CC(rand(ComplexF64))

function gamma_homotopy(G, F, t; mode::Symbol=:diagonal_start)
    mode == :none && return [(1 - t) * G[i] + t * F[i] for i in eachindex(G)]
    if mode == :scalar_start
        gamma = random_gamma()
        return [(1 - t) * gamma * G[i] + t * F[i] for i in eachindex(G)]
    end
    mode == :diagonal_start && return [
        (1 - t) * random_gamma() * G[i] + t * F[i] for i in eachindex(G)
    ]
    error("gamma mode must be none, scalar_start, or diagonal_start")
end

function make_benchmark(spec::BenchmarkSpec; gamma_mode::Symbol=:diagonal_start)
    Random.seed!(spec.seed)
    ring, variables = polynomial_ring(spec.n)
    x = variables[1:spec.n]
    _, t = ring["t"]
    if spec.family == :katsura
        G, F = katsura_start_system(x), katsura_system(x)
        roots = katsura_start_roots(spec.n)
    elseif spec.family == :random_dense
        G, F = bezout_start_system(x, spec.degree), random_dense_system(x, spec.degree)
        roots = bezout_roots(spec.n, spec.degree)
    else
        error("unknown benchmark family: $(spec.family)")
    end
    return gamma_homotopy(G, F, t; mode=gamma_mode), roots
end

function benchmark_specs(; degree=2, dimensions=3:6)
    specs = BenchmarkSpec[]
    append!(specs, [BenchmarkSpec(:katsura, n, 2, 10_000 + n) for n in dimensions])
    append!(specs, [BenchmarkSpec(:random_dense, n, degree, 20_000 + n) for n in dimensions])
    return specs
end
