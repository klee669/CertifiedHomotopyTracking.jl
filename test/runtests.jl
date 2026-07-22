using Test
using Random
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

    K = krawczyk_operator(sys, start, CC(0), 1e-3)
    passed, k_norm = krawczyk_test(sys, start, CC(0), 1e-3)
    @test passed
    @test norm_inf(K) == k_norm
end

@testset "ACB inverse coercion and shape checks" begin
    CC = AcbField(128)

    mixed = Any[1 0; 0 CC(2)]
    inv_mixed = inv_acb(mixed, CC)
    @test inv_mixed isa Matrix{AcbFieldElem}
    @test iszero(inv_mixed[1, 1] - CC(1))
    @test iszero(inv_mixed[2, 2] - inv(CC(2)))

    inferred = inv_acb(Any[CC(3) 0; 0 1])
    @test abs(convert_to_double_int(inferred[1, 1]) - 1 / 3) < 1e-12
    @test abs(convert_to_double_int(inferred[2, 2]) - 1) < 1e-12

    @test_throws DimensionMismatch inv_acb(Any[CC(1) 0], CC)
    @test_throws ArgumentError inv_acb(Any[1 0; 0 1])
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
    K = krawczyk_operator(curve, p, frame, 1e-3, 1e-3)
    passed, k_norm = krawczyk_test(curve, p, frame, 1e-3, 1e-3)
    @test passed
    @test norm_inf(K) == k_norm

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

@testset "Hot-path arithmetic is bit-identical to the unfused formulas" begin
    CHT = CertifiedHomotopyTracking
    CC = AcbField(128)
    RR = ArbField(128)

    # Two balls are the same enclosure only if midpoint *and* radius agree on
    # both components. `==` on AcbFieldElem is a containment-style predicate,
    # so it would not catch a widened remainder.
    same(a, b) =
        Nemo.midpoint(real(a)) == Nemo.midpoint(real(b)) &&
        Nemo.radius(real(a)) == Nemo.radius(real(b)) &&
        Nemo.midpoint(imag(a)) == Nemo.midpoint(imag(b)) &&
        Nemo.radius(imag(a)) == Nemo.radius(imag(b))

    ref_interval(h) = (upper = Nemo.midpoint(h) + Nemo.radius(h);
                       half = RR(upper / 2);
                       CC(Nemo.ball(half, half), RR(0)))
    ref_poly(c0, c1, c2, c3, t) = c0 + t * (c1 + t * (c2 + t * c3))

    function ref_mul(a, b)
        t = ref_interval(a.h)
        C0 = a.c0 * b.c0
        C1 = a.c0 * b.c1 + a.c1 * b.c0
        C2 = a.c0 * b.c2 + a.c1 * b.c1 + a.c2 * b.c0
        C3 = a.c0 * b.c3 + a.c1 * b.c2 + a.c2 * b.c1 + a.c3 * b.c0
        term4 = a.c1 * b.c3 + a.c2 * b.c2 + a.c3 * b.c1
        term5 = a.c2 * b.c3 + a.c3 * b.c2
        term6 = a.c3 * b.c3
        t2 = t * t; t4 = t2 * t2; t5 = t4 * t; t6 = t4 * t2
        trunc_error = term4 * t4 + term5 * t5 + term6 * t6
        pA = ref_poly(a.c0, a.c1, a.c2, a.c3, t)
        pB = ref_poly(b.c0, b.c1, b.c2, b.c3, t)
        return (C0, C1, C2, C3, trunc_error + (pA * b.rem + pB * a.rem + a.rem * b.rem))
    end

    rng = MersenneTwister(20260722)
    mk(exact, h) = CHT.TaylorModel3(
        CC(randn(rng), randn(rng)), CC(randn(rng), randn(rng)),
        CC(randn(rng), randn(rng)), CC(randn(rng), randn(rng)),
        exact ? CC(0) : CC(RR("0.01 +/- 0.001"), RR("-0.02 +/- 0.001")),
        h,
    )

    # Both remainder-zero branches and the general branch must all agree.
    for exact_a in (true, false), exact_b in (true, false)
        h = RR(0.05)
        a = mk(exact_a, h)
        b = mk(exact_b, h)
        got = a * b
        want = ref_mul(a, b)
        @test same(got.c0, want[1])
        @test same(got.c1, want[2])
        @test same(got.c2, want[3])
        @test same(got.c3, want[4])
        @test same(got.rem, want[5])

        s = CC(randn(rng), randn(rng))
        m = CHT.get_mid(s)
        scaled = a * s
        want_rem = a.rem * s + ref_poly(a.c0, a.c1, a.c2, a.c3, ref_interval(h)) * (s - m)
        @test same(scaled.c0, a.c0 * m)
        @test same(scaled.rem, want_rem)

        @test same(CHT.evaluate_taylor(a), ref_poly(a.c0, a.c1, a.c2, a.c3, ref_interval(h)) + a.rem)
    end

    @variables x y t
    sys = straight_line_homotopy([x^2 + 3y - 4, y^2 + 3], [x^2 - 1, y^2 - 1], [x, y];
                                 CCRing=CC, gamma=CC(0.5, 0.5))

    # krawczyk_operator fuses -(A*fx)/r + (I - A*Jx)*B; compare against the
    # unfused expression built from the same pieces.
    for trial in 1:20
        pt = [CC(randn(rng), randn(rng)), CC(randn(rng), randn(rng))]
        tv = CC(rand(rng))
        r = 10.0^(-rand(rng, 1:6))
        A = CHT.compute_preconditioner(sys, pt, tv)

        B = CHT._acb_unit_box_vector(CC, RR, 2)
        fx = evaluate_H(sys, pt, tv)
        Jx = evaluate_Jac(sys, pt .+ (B .* CC(r)), tv)
        want = (-(A * fx) ./ CC(r)) + (CHT._acb_identity_matrix(CC, 2) - A * Jx) * B

        got = krawczyk_operator(sys, pt, tv, r, A)
        @test all(same(got[i], want[i]) for i in 1:2)
    end

    # The fused validation path must match the instrumented one exactly.
    cache = CHT.KrawczykValidationCache(CC, RR, 2)
    for trial in 1:20
        h = 0.01
        X_tm = [CHT.TaylorModel3(CC(randn(rng), randn(rng)), CC(randn(rng), randn(rng)),
                                 CC(randn(rng), randn(rng)), CC(randn(rng), randn(rng)),
                                 CC(0), RR(h)) for _ in 1:2]
        t_start = 0.3
        r = 1e-4
        A = CHT.compute_preconditioner(sys, [tm.c0 for tm in X_tm], CC(t_start))
        fused = CHT.validate_step_taylor3_diagnostics(sys, X_tm, t_start, h, r, A;
                                                      cache=cache, profile_validation=false)
        instrumented = CHT.validate_step_taylor3_diagnostics(sys, X_tm, t_start, h, r, A;
                                                            cache=cache, profile_validation=true)
        @test fused.passed == instrumented.passed
        @test fused.norm_K == instrumented.norm_K
        @test fused.Y == instrumented.Y
        @test fused.Z == instrumented.Z
        @test fused.yz_bound == instrumented.yz_bound
    end
end
