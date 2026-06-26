using Test
using CertifiedHomotopyTracking
import HomotopyContinuation

@testset "CertifiedHomotopyTracking" begin
    @test true
end

@testset "Edge homotopy complex constants" begin
    @variables x a
    CC = AcbField(128)

    F = [x + a + (1 + 2im)]
    compiled = compile_edge_homotopy(F, [x], [a])
    sys = make_edge_system(compiled, [CC(0)], [CC(1)])
    start = [CC(-1, -2)]

    @test length(compiled.fixed_const_values) == 1
    @test length(sys.p_const) == 1
    @test iszero(evaluate_H(sys, start, CC(0))[1])
    @test iszero(evaluate_Jac(sys, start, CC(0))[1, 1] - CC(1))
    @test iszero(evaluate_dt(sys, start, CC(0))[1] - CC(1))
end

@testset "Taylor model real domain interval" begin
    RR = ArbField(128)
    CC = AcbField(128)
    h = RR(0.1)
    interval = CertifiedHomotopyTracking._tm_real_interval(CC, h)

    @test contains(real(interval), RR(0))
    @test contains(real(interval), h)
    @test contains(imag(interval), RR(0))
    @test Nemo.radius(imag(interval)) == 0
end

@testset "Direct compiled homotopy tracking" begin
    @variables x t
    CC = AcbField(128)

    compiled = compile_homotopy([(1 - t) * (x^2 - 1) + t * (x^2 - 4)], [x], t)
    res = track_path(compiled, [CC(1)]; h_init=0.05)

    @test success(res)
    @test abs(solution(res)[1] - 2) < 1e-8

    @variables y
    constant_path = compile_homotopy(
        [
            (1 - t) * (x^2 - 1) + t * (x^2 + y - 2),
            (1 - t) * (y^2 - 1) + t * (y^2 + x - 2),
        ],
        [x, y],
        t,
    )
    constant_res = track_path(constant_path, [CC(1), CC(1)])

    @test success(constant_res)
    @test all(abs.(solution(constant_res) .- [1, 1]) .< 1e-8)

    gamma_path = compile_homotopy(
        [
            (1 - t) * (1 + im) * (x^3 - 1) + t * (x^3 + 2y - 3),
            (1 - t) * (y^2 - 1) + t * (y^2 + x - 2),
        ],
        [x, y],
        t,
    )
    gamma_res = track_path(gamma_path, [CC(1), CC(1)])

    @test success(gamma_res)
    @test all(abs.(solution(gamma_res) .- [1, 1]) .< 1e-8)
end

@testset "Static variety system wrapper" begin
    @variables x y
    CC = AcbField(128)

    compiled = compile_system([x^2 + y^2 - 1], [x, y])
    sys = SpecializedHomotopy(compiled, CC)
    p = [CC(1), CC(0)]

    @test evaluate_H(sys, p, CC(0))[1] == 0
    @test size(evaluate_Jac(sys, p, CC(0))) == (1, 2)

    curve = variety_system([x^2 + (1 + im) * y^2 - 1], [x, y]; CCRing=CC)
    @test system(curve) isa SpecializedHomotopy
    @test evaluate_system(curve, p)[1] == 0
    @test size(jacobian_system(curve, p)) == (1, 2)

    frame = CertifiedHomotopyTracking.local_tangent_normal_frame(curve, p)
    @test frame.dim == 1
    @test frame.rank == 1
    passed, _ = krawczyk_test(curve, p, frame, 1e-3, 1e-3)
    @test passed

    box = refine_moore_box(curve, [CC(1.0), CC(1e-8)], 1e-3)
    @test box.success
    @test box.frame.dim == 1
    @test box.normal_radius > 1e-3

    square = variety_system([x^2 - 1], [x]; CCRing=CC)
    root_box = refine_moore_box(square, [CC(1.01)], 0.1)
    @test root_box.success
    @test root_box.frame.dim == 0
    @test abs(convert_to_double_int(root_box.center[1]) - 1) < 1e-8

    curve_approx = certified_variety_approximation(
        variety_system([x^2 + y^2 - 1], [x, y]; CCRing=CC),
        [CC(1), CC(0)];
        tangent_radius=1e-3,
        normal_radius=1e-3,
        max_boxes=5,
    )
    @test length(curve_approx.boxes) == 5
    @test all(box -> box.success && box.frame.dim == 1, curve_approx.boxes)

    @variables z
    surface_approx = certified_variety_approximation(
        variety_system([x^2 + y^2 + z^2 - 1], [x, y, z]; CCRing=CC),
        [CC(1), CC(0), CC(0)];
        tangent_radius=1e-3,
        normal_radius=1e-3,
        max_boxes=7,
    )
    @test length(surface_approx.boxes) == 7
    @test all(box -> box.success && box.frame.dim == 2, surface_approx.boxes)
end

@testset "HomotopyContinuation numerical trace" begin
    HomotopyContinuation.@var x
    G = HomotopyContinuation.System([x^2 - 1])
    F = HomotopyContinuation.System([x^2 - 2])
    H = HomotopyContinuation.StraightLineHomotopy(G, F)
    x_start = [1.0 + 0im]

    out = collect_hc_trace(H, x_start; t_start = 1.0, t_target = 0.0)

    @test length(out.trace) >= 2
    @test isapprox(first(out.trace).t, 1.0; atol = 1e-12)
    if out.success
        @test isapprox(last(out.trace).t, 0.0; atol = 1e-12)
    end
    @test all(length(point.x) == 1 for point in out.trace)
    @test out.success == HomotopyContinuation.is_success(out.status)
    @test out.accepted_steps isa Union{Int,Missing}
    @test out.rejected_steps isa Union{Int,Missing}
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

    H_projective = straight_line_homotopy(F, G, [x]; CCRing=CC, projective=true)
    res_projective = track_path(H_projective, start; h_init=0.05)

    compiled_projective = compile_edge_homotopy(F, [x], Num[]; projective=true)
    sys_projective = SpecializedHomotopy(compiled_projective, CC)
    direct_projective = track_path(sys_projective, start; h_init=0.05)

    @test success(res_affine)
    @test success(res_projective)
    @test success(direct_projective)
    @test length(input_start(res_projective)) == 1
    @test length(refined_start(res_projective)) == 1
    @test length(projective_input_start(res_projective)) == 2
    @test length(projective_refined_start(res_projective)) == 2
    @test isapprox(solution(res_affine)[1], solution(res_projective)[1]; atol=1e-8)
    @test isapprox(solution(res_affine)[1], solution(direct_projective)[1]; atol=1e-8)

    F_large = [x - 100]
    H_large_affine = straight_line_homotopy(F_large, G, [x]; CCRing=CC)
    large_affine = track_path(H_large_affine, start; h_init=0.05)

    H_large_projective = straight_line_homotopy(F_large, G, [x]; CCRing=CC, projective=true)
    large_projective = track_path(H_large_projective, start; h_init=0.05)

    @test success(large_affine)
    @test success(large_projective)
    @test abs(solution(large_affine)[1]) > 50
    @test large_projective.patch_idx == 2
    @test maximum(abs.(convert_to_double_int.(projective_solution(large_projective)))) <= 2
    @test isapprox(solution(large_projective)[1], 100 + 0im; atol=1e-6)
end

@testset "Adaptive precision tracking" begin
    @variables x
    CC = AcbField(128)
    H = straight_line_homotopy([x - 2], [x - 1], [x]; CCRing=CC)
    start = [CC(1)]

    sys53 = system_with_precision(H, 53)
    @test precision(sys53.CC) == 53
    @test precision(H.CC) == 128

    adaptive = track_path(H, start; h_init=0.05)
    fixed = track_path(H, start; h_init=0.05, adaptive_precision=false)
    promoted = track_path(H, start; h_init=0.05, rho=0.1, max_precision=106, precision_rejection_threshold=1)

    @test success(adaptive)
    @test adaptive.initial_precision == 53
    @test adaptive.final_precision == 53
    @test precision(parent(adaptive.root[1])) == 128

    @test success(fixed)
    @test fixed.initial_precision == 128
    @test fixed.final_precision == 128

    @test !success(promoted)
    @test promoted.status == :step_too_small
    @test promoted.initial_precision == 53
    @test promoted.final_precision == 106
    @test precision(parent(promoted.root[1])) == 128
end

@testset "Path visualization export" begin
    @variables x
    CC = AcbField(128)
    H = straight_line_homotopy([x - 2], [x - 1], [x]; CCRing=CC)
    res = track_path(H, [CC(1)]; h_init=0.05, adaptive_precision=false, visualize=true)

    @test success(res)
    @test !isempty(path_boxes(res))

    tikz_file = tempname() * ".tex"
    tikz3_file = tempname() * ".tex"
    obj_file = tempname() * ".obj"
    @test export_path_tikz(res, tikz_file) == tikz_file
    @test export_path_tikz(res, tikz3_file; axes=(:t, 1, (1, :imag))) == tikz3_file
    @test export_path_obj(res, obj_file) == obj_file
    @test isfile(tikz_file)
    @test isfile(tikz3_file)
    @test isfile(obj_file)
    @test occursin("\\documentclass[tikz,border=4pt]{standalone}", read(tikz_file, String))
    @test occursin("\\begin{tikzpicture}", read(tikz_file, String))
    @test occursin("z={(0cm,1cm)}", read(tikz3_file, String))
    @test startswith(read(obj_file, String), "# CertifiedHomotopyTracking.jl path boxes")

    auto_file = tempname() * ".tex"
    auto_res = track_path(
        H,
        [CC(1)];
        h_init=0.05,
        adaptive_precision=false,
        visualize=true,
        visualize_options=(; filename=auto_file, axes=(:t, 1), color="red"),
    )
    @test success(auto_res)
    @test isfile(auto_file)

    trace_file = tempname() * ".tex"
    accepted_boxes = [box for box in path_boxes(res) if get(box.metadata, :stage, nothing) == :accepted_step]
    box = first(accepted_boxes)
    trace_viz = PathVisualization(
        [box],
        (:t, 1),
        :posteriori,
        (; trace_points = [box, last(accepted_boxes)]),
    )
    @test export_path_tikz(trace_viz, trace_file; show_trace=true) == trace_file
    @test occursin("\\draw[black, line width=0.55pt]", read(trace_file, String))
end
