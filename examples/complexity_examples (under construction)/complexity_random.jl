include("../src/certified_monodromy_computation.jl")
include("../src/complexity_test.jl")


CCi = AcbField()

function create_polynomial_ring(CCi, n)
    # Generate variable names x_1, x_2, ..., x_n
    var_names = ["x_$i" for i in 1:n]
    push!(var_names, "η")  # Add η at the end
    
    # Create the polynomial ring with the generated variables
    return CCi[var_names...]
end

n = 6
d = 2
R, xvars = create_polynomial_ring(CCi, n)
x_vars = xvars[1:n]    # Get x variables
η = xvars[end]         # Get η variable
HR, (t) = R["t"]


G = transpose(map(i -> (gens(R))[i]^d-1, 1:n))
function generate_all_roots(CCi, n,d)
    # Generate n-th roots of unity
    single_roots = [
        CCi(cos(2*k*π/d) + im*sin(2*k*π/d))
        for k in 0:d-1
    ]
    
    # Generate all possible combinations of length n
    return [[x...] for x in Iterators.product(fill(single_roots, n)...)]
end

all_solutions = vec(generate_all_roots(CCi,n,d))
F = transpose(hcat(generate_random_dense_system(n, d)))
G = transpose(map(i -> (gens(R))[i]^d-1, 1:n))

    H = ([(1-t)*(CCi(rand(ComplexF64))); (1-t)*(CCi(rand(ComplexF64))); (1-t)*(CCi(rand(ComplexF64))); (1-t)*(CCi(rand(ComplexF64)))]*G+[t; t; t; t]*F)[1,:]
    H1 = map(i -> derivative(i), H)
    H1 = evaluate_matrix(Matrix(transpose(hcat(H1))), 0)
    
info_list = []
for root in all_solutions[1:2]    
    sol1, it1 = tracking_certified_predictor(H, H1, root, .1; iterations_count = true)
    push!(info_list, it1)
end

# 1. Vector{Any}를 Matrix로 변환 (각 내부 벡터를 행으로 간주)
    data_matrix = vcat(transpose.(info_list)...);

# 2. 열 방향(차원 1)으로 평균 계산
column_means = mean(data_matrix, dims=1)
#[iter, minimum(dt_list), median(dt_list),maximum(dt_list),  quantile(dt_list, 0.25), quantile(dt_list, 0.75), minimum(r_list), median(r_list), mean(b_list), maximum(eta_list), mean(ratio_list[2:end])]
map(i -> log(10,i), data_matrix[end,:])
println("각 열의 평균 (Column Mean):")
# 결과를 보기 쉽게 출력
for (i, mean_val) in enumerate(column_means[1, :])
    println("Column $i: $mean_val")
end

all_solutions = vec(generate_all_roots(CCi,n,d))

iters1 = Int[]
iters2 = Int[]
iters3 = Int[]
iters4 = Int[] 
iters5 = Int[]

for i in 1:100

    random_system = transpose(hcat(generate_random_dense_system(n, d)))
    H = ([(1-t)*(1+onei(CCi)); (1-t)*(CCi(rand(ComplexF64))); (1-t)*(CCi(rand(ComplexF64))); (1-t)*(CCi(rand(ComplexF64)))]*G+[t; t; t; t]*random_system)[1,:]
    point = rand(all_solutions)
    H1 = map(i -> derivative(i), H)
    H1 = evaluate_matrix(Matrix(transpose(hcat(H1))), 0)


    sol1, it1 = tracking_certified_predictor(H, H1, point, .1; iterations_count = true)
    sol2, it2 = tracking_without_predictor(H,point,.1; iterations_count = true)
    sol3, it3 = track(H,point,.1; iterations_count = true)
    sol4, it4 = track(H,point,.1; iterations_count = true, tracking = "non-truncate")
    sol5, it5 = tracking_certified_hermite_predictor(H,point,.1; iterations_count = true)
    push!(iters1, it1)
    push!(iters2, it2)
    push!(iters3, it3)
    push!(iters4, it4)
    push!(iters5, it5)
end
mean(iters1), mean(iters2), mean(iters3), mean(iters4), mean(iters5)


