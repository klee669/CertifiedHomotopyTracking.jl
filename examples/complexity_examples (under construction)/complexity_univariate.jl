include("../src/certified_monodromy_computation.jl")
include("../src/complexity_test.jl")


CCi = AcbField()
R, (a,η) = CCi["a","η"]
vvars = (a,η)
#R, vars = create_polynomial_ring(CCi, n)
HR, (t) = R["t"]



m=10000

F = [a^2-m]
G = [a^2-1]

F = transpose(hcat(F))
G = transpose(hcat(G))
    H = ([(1-t); (1-t); (1-t); (1-t)]*G+[t; t; t; t]*F)[1,:]
    H1 = map(i -> derivative(i), H)
    H1 = evaluate_matrix(Matrix(transpose(hcat(H1))), 0)

x,it = tracking_certified_predictor(H, H1, [CCi(1.0)], .1; iterations_count = true)
it[1]/(sqrt(m)-1)
#m = 10 => 27.286042438881594
#m= 100 => 26.88888888888889
#m= 1000 => 26.74479883561464
#m=10000 => 26.696969696969695
#m=20000 => 26.691096713707484
#m=30000 => 26.683300979295947