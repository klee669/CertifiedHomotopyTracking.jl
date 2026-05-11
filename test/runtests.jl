using Test
using CertifiedHomotopyTracking

@testset "CertifiedHomotopyTracking" begin
    @test true
end

@testset "Homogenization utilities" begin
    @variables x y h

    f = x^2 + x * y + y + 3
    F = homogenize_expr(f, [x, y], h)
    @test isequal(Symbolics.simplify(F - (x^2 + x * y + y * h + 3h^2)), Num(0))

    exprs = [x + 1, x^3 + y]
    H = homogenize_system(exprs, [x, y], h)
    @test isequal(Symbolics.simplify(H[1] - (x + h)), Num(0))
    @test isequal(Symbolics.simplify(H[2] - (x^3 + y * h^2)), Num(0))

    c = homogenize_expr(7, [x, y], h)
    @test isequal(Symbolics.simplify(c - 7), Num(0))

    complex_hom = homogenize_expr((1 + 2im) * x^2 + y, [x, y], h)
    complex_residual = Symbolics.simplify(complex_hom - ((1 + 2im) * x^2 + y * h))
    @test isequal(Symbolics.simplify(real(complex_residual)), Num(0))
    @test isequal(Symbolics.simplify(imag(complex_residual)), Num(0))

    rational_hom = homogenize_expr(1 / x + y, [x, y], h)
    @test isequal(Symbolics.simplify(rational_hom - (h^2 + x * y)), Num(0))

    cancelled_hom = homogenize_expr((x + 1) / (x + 1), [x, y], h)
    @test isequal(Symbolics.simplify(cancelled_hom - 1), Num(0))

    @test_throws ArgumentError homogenize_expr(sin(x) + y, [x, y], h)
    @test_throws ArgumentError homogenize_expr(sin(x) / y + 1, [x, y], h)
end

@testset "Affine patch utilities" begin
    a = random_patch_vector(4)
    @test length(a) == 4
    @test isapprox(norm(a), 1.0; atol=1e-12)

    b = random_patch_vector(3; normalize=false)
    @test length(b) == 3

    @variables X0 X1 X2
    patch = patch_equation((X0, X1, X2), (1, 2, 3))
    @test isequal(Symbolics.simplify(patch - (X0 + 2X1 + 3X2 - 1)), Num(0))

    x_affine = [2.0, -1.0]
    a_real = [0.5, -0.25, 0.75]
    X = lift_to_patch(x_affine, a_real)
    @test isapprox(sum(a_real .* X), 1.0; atol=1e-12)
    @test isapprox(X[2] / X[1], x_affine[1]; atol=1e-12)
    @test isapprox(X[3] / X[1], x_affine[2]; atol=1e-12)

    X0_val = [2 + im, -3 + 2im, 4 - im]
    a_complex = [1 - im, 0.25 + 0.5im, -0.5im]
    X_repatched = repatch(X0_val, a_complex)
    @test isapprox(sum(a_complex .* X_repatched), 1 + 0im; atol=1e-12)
    @test isapprox(X_repatched[2] / X_repatched[1], X0_val[2] / X0_val[1]; atol=1e-12)
    @test isapprox(X_repatched[3] / X_repatched[1], X0_val[3] / X0_val[1]; atol=1e-12)

    @test_throws ArgumentError lift_to_patch([1.0], [-1.0, 1.0])
    @test_throws ArgumentError repatch([1.0, 2.0], [2.0, -1.0])
end

@testset "Projective tracking" begin
    @variables x
    CC = AcbField(128)

    F = [x - 2]
    G = [x - 1]
    start = [CC(1)]

    H_affine = straight_line_homotopy(F, G, [x]; CCRing=CC)
    res_affine = track_path(H_affine, start; h_init=0.05)

    H_projective = straight_line_homotopy(F, G, [x]; CCRing=CC, projective=true, patch_vector=[1, 1])
    res_projective = track_path(H_projective, start; h_init=0.05)

    compiled_projective = compile_edge_homotopy(F, [x], Num[]; projective=true, patch_vector=[1, 1])
    sys_projective = HCSystem(compiled_projective, CC)
    direct_projective = track_path(sys_projective, start; h_init=0.05)

    @test succeeded(res_affine)
    @test succeeded(res_projective)
    @test succeeded(direct_projective)
    @test length(input_start(res_projective)) == 1
    @test length(refined_start(res_projective)) == 1
    @test length(projective_input_start(res_projective)) == 2
    @test length(projective_refined_start(res_projective)) == 2
    @test isapprox(solution(res_affine)[1], solution(res_projective)[1]; atol=1e-8)
    @test isapprox(solution(res_affine)[1], solution(direct_projective)[1]; atol=1e-8)

    F_large = [x - 100]
    H_large_affine = straight_line_homotopy(F_large, G, [x]; CCRing=CC)
    large_affine = track_path(H_large_affine, start; h_init=0.05)

    H_large_projective = straight_line_homotopy(F_large, G, [x]; CCRing=CC, projective=true, patch_vector=[1, 1])
    large_projective = track_path(H_large_projective, start; h_init=0.05)

    @test succeeded(large_affine)
    @test succeeded(large_projective)
    @test abs(solution(large_affine)[1]) > 50
    @test maximum(abs.(convert_to_double_int.(projective_solution(large_projective)))) <= 2
    @test isapprox(solution(large_projective)[1], 100 + 0im; atol=1e-6)
end
