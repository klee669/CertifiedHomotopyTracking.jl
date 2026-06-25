using CertifiedHomotopyTracking

CC = AcbField(128)

# Example 1: a plane curve, x^2 + y^2 - 1 = 0.
@variables x y
curve = variety_system([x^2 + y^2 - 1], [x, y]; CCRing=CC)

p_curve = [CC(1), CC(0)]
F_curve = evaluate_system(curve, p_curve)
J_curve = jacobian_system(curve, p_curve)

println("curve = ", curve)
println("F_curve = ", F_curve)
println("J_curve = ", J_curve)

# You can still access the underlying SpecializedHomotopy if you want the lower-level API.
sys_curve = system(curve)
F_curve_low_level = evaluate_H(sys_curve, p_curve, CC(0))
J_curve_low_level = evaluate_Jac(sys_curve, p_curve, CC(0))

println("F_curve_low_level = ", F_curve_low_level)
println("J_curve_low_level = ", J_curve_low_level)

# Example 2: a surface, x^2 + y^2 + z^2 - 1 = 0.
@variables z
surface = variety_system([x^2 + y^2 + z^2 - 1], [x, y, z]; CCRing=CC)

p_surface = [CC(1), CC(0), CC(0)]
F_surface = evaluate_system(surface, p_surface)
J_surface = jacobian_system(surface, p_surface)

println("surface = ", surface)
println("F_surface = ", F_surface)
println("J_surface = ", J_surface)

# Example 3: complex coefficients are parameterized using CHT's existing path.
complex_curve = variety_system([x^2 + (1 + 2im) * y^2 - 1], [x, y]; CCRing=CC)
p_complex = [CC(1), CC(0)]

println("complex_curve = ", complex_curve)
println("F_complex = ", evaluate_system(complex_curve, p_complex))
println("J_complex = ", jacobian_system(complex_curve, p_complex))

# Example 4: local tangent-normal frame and one certified local box.
frame = CertifiedHomotopyTracking.local_tangent_normal_frame(curve, p_curve)
println("curve rank = ", frame.rank)
println("curve dimension = ", frame.dim)
println("curve tangent frame = ", frame.tangent)
println("curve normal frame = ", frame.normal)

passed, k_norm = krawczyk_test(curve, p_curve, frame, 1e-3, 1e-3)
println("normal Krawczyk passed = ", passed)
println("normal Krawczyk norm = ", k_norm)

normal_box = refine_moore_box(curve, [CC(1.0), CC(1e-8)], 1e-3)
println("normal Moore box success = ", normal_box.success)
println("normal Moore box center = ", normal_box.center)
println("normal Moore box radius = ", normal_box.normal_radius)
println("normal Moore box Krawczyk norm = ", normal_box.krawczyk_norm)

# Example 5: dimension-generic paving loop.
curve_approx = certified_variety_approximation(
    curve,
    p_curve;
    tangent_radius=1e-3,
    normal_radius=1e-3,
    max_boxes=5,
)
println("curve approximation boxes = ", length(curve_approx.boxes))
println("curve approximation attempted = ", curve_approx.attempted)
println("curve approximation rejected = ", curve_approx.rejected)
println("curve approximation skipped overlap = ", curve_approx.skipped_overlap)

surface_approx = certified_variety_approximation(
    surface,
    p_surface;
    tangent_radius=1e-3,
    normal_radius=1e-3,
    max_boxes=7,
)
println("surface approximation boxes = ", length(surface_approx.boxes))
println("surface approximation attempted = ", surface_approx.attempted)
println("surface approximation rejected = ", surface_approx.rejected)
println("surface approximation skipped overlap = ", surface_approx.skipped_overlap)
