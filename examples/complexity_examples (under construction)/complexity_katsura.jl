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
R, (a,b,c,d,e,f,η) = CCi["a","b","c","d","e","f","η"]
vvars = (a,b,c,d,e,f,η)
#R, vars = create_polynomial_ring(CCi, n)
x_vars = vvars[1:n]    # Get x variables
η = vvars[end]         # Get η variable
HR, (t) = R["t"]




function katsura_start_system(n::Int)
    start_polys = []
    
    # First equation: x[1] - 1 = 0
    push!(start_polys, x_vars[1] - 1)
    
    # Rest: x[i]^2 - 1 = 0 for i = 2, ..., n
    for i in 2:n
        push!(start_polys, x_vars[i]^2 - 1)
    end
    
    return start_polys, x_vars
end
function katsura_start_roots(CCi, n)
    solutions = []
    
    # Generate 2^(n-1) solutions
    for bits in 0:2^(n-1)-1
        sol = Vector{typeof(CCi(1.0))}(undef, n)
        
        # First coordinate: always 1
        sol[1] = CCi(1.0)
        
        # Remaining coordinates: ±1 based on bit pattern
        for i in 2:n
            if (bits >> (i-2)) & 1 == 1
                sol[i] = CCi(1.0)
            else
                sol[i] = CCi(-1.0)
            end
        end
        
        push!(solutions, sol)
    end
    
    return solutions
end

# Example usage
# assuming CCi is your complex interval type
roots = katsura_start_roots(CCi, n)
F = [a+2*b+2*c+2*d+2*e+2*f-1, a^2+2*b^2+2*c^2+2*d^2+2*e^2+2*f^2-a, 2*a*b+2*b*c+2*c*d+2*d*e+2*e*f-b, b^2+2*a*c+2*b*d+2*c*e+2*d*f-c, 2*b*c+2*a*d+2*b*e+2*c*f-d, c^2+2*b*d+2*a*e+2*b*f-e]#[a+2*b+2*c+2*d+2*e-1, a^2+2*b^2+2*c^2+2*d^2+2*e^2-a, 2*a*b+2*b*c+2*c*d+2*d*e-b, b^2+2*a*c+2*b*d+2*c*e-c, 2*b*c+2*a*d+2*b*e-d]#[a+2*b+2*c+2*d-1, a^2+2*b^2+2*c^2+2*d^2-a, 2*a*b+2*b*c+2*c*d-b, b^2+2*a*c+2*b*d-c]#[a+2*b+2*c-1, a^2+2*b^2+2*c^2-a, 2*a*b+2*b*c-b]
G, x_vars = katsura_start_system(n)

F = transpose(hcat(F))
G = transpose(hcat(G))
    H = ([(1-t)*(CCi(rand(ComplexF64))); (1-t)*(CCi(rand(ComplexF64))); (1-t)*(CCi(rand(ComplexF64))); (1-t)*(CCi(rand(ComplexF64)))]*G+[t; t; t; t]*F)[1,:]
    H1 = map(i -> derivative(i), H)
    H1 = evaluate_matrix(Matrix(transpose(hcat(H1))), 0)
    sol1, info1 = tracking_certified_predictor(H, H1, roots[1], .1; iterations_count = true)
    sol2, info2 = tracking_certified_predictor(H, H1, roots[2], .1; iterations_count = true)
    sol3, info3 = tracking_certified_predictor(H, H1, roots[3], .1; iterations_count = true)
    sol4, info4 = tracking_certified_predictor(H, H1, roots[4], .1; iterations_count = true)

info_list = []
for root in roots    
    sol1, it1 = tracking_certified_predictor(H, H1, root, .1; iterations_count = true)
    push!(info_list, it1)
end

# 1. Vector{Any}를 Matrix로 변환 (각 내부 벡터를 행으로 간주)
data_matrix = vcat(transpose.(info_list)...)

# 2. 열 방향(차원 1)으로 평균 계산
column_means = mean(data_matrix, dims=1)
#[iter, minimum(dt_list), median(dt_list),maximum(dt_list),  quantile(dt_list, 0.25), quantile(dt_list, 0.75), minimum(r_list), median(r_list), mean(b_list), maximum(eta_list), mean(ratio_list[2:end])]
map(i -> log(10,i), data_matrix[1,:])
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


