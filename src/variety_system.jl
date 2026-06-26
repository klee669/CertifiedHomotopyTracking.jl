export AlgebraicVarietySystem, VarietyBox, VarietyApproximation,
       variety_system, system, evaluate_system, jacobian_system,
       certified_variety_approximation, export_variety_obj

"""
    VarietyFrame

Local tangent-normal decomposition for an [`AlgebraicVarietySystem`](@ref).

Fields:

- `tangent`: basis for the local tangent directions.
- `normal`: basis for the normal directions.
- `equation_basis`: selected equation-space basis from the Jacobian SVD.
- `preconditioner`: inverse compressed normal Jacobian used by Krawczyk tests.
- `rank`: numerical Jacobian rank.
- `dim`: local variety dimension.
"""
struct VarietyFrame
    tangent::Matrix{AcbFieldElem}
    normal::Matrix{AcbFieldElem}
    equation_basis::Matrix{AcbFieldElem}
    preconditioner::Matrix{AcbFieldElem}
    rank::Int
    dim::Int
end

"""
    VarietyBox

Certified local box for a variety.

Fields:

- `center`: center point.
- `frame`: internal local tangent-normal data used by the certification boxes.
- `tangent_radius`: radius in tangent directions.
- `normal_radius`: radius in normal directions.
- `krawczyk_norm`: final Krawczyk norm.
- `success`: whether certification/refinement succeeded.
"""
struct VarietyBox
    center::Vector{AcbFieldElem}
    frame::VarietyFrame
    tangent_radius::Float64
    normal_radius::Float64
    krawczyk_norm::Float64
    success::Bool
end

"""
    VarietyApproximation

Collection of certified variety boxes returned by
[`certified_variety_approximation`](@ref).

Fields:

- `boxes`: accepted [`VarietyBox`](@ref)s.
- `attempted`: number of candidate boxes attempted.
- `rejected`: number of failed candidates.
- `skipped_overlap`: number of candidates skipped because they overlapped existing boxes.
- `max_boxes`: requested maximum number of boxes.
"""
struct VarietyApproximation
    boxes::Vector{VarietyBox}
    attempted::Int
    rejected::Int
    skipped_overlap::Int
    max_boxes::Int
end

"""
    AlgebraicVarietySystem

Wrapper for a static polynomial system treated as an algebraic variety.

Use [`variety_system`](@ref) to build one from Symbolics expressions. The
underlying [`SpecializedHomotopy`](@ref) is available through [`system`](@ref).
"""
struct AlgebraicVarietySystem
    system::SpecializedHomotopy
    t::AcbFieldElem
end

function AlgebraicVarietySystem(sys::SpecializedHomotopy)
    return AlgebraicVarietySystem(sys, sys.CC(0))
end

"""
    system(variety::AlgebraicVarietySystem) -> SpecializedHomotopy

Return the underlying [`SpecializedHomotopy`](@ref).
"""
system(variety::AlgebraicVarietySystem) = variety.system

function Base.show(io::IO, variety::AlgebraicVarietySystem)
    print(io, "AlgebraicVarietySystem(", variety.system, ")")
end

function _as_acb_vector(CC::AcbField, x)
    return AcbFieldElem[CC(xi) for xi in x]
end

"""
    evaluate_system(variety, x)

Evaluate the defining equations of `variety` at `x`.

# Example

```julia
@variables x y z;
CC = AcbField(128);
surface = variety_system([x^2 + y^2 + z^2 - 1], [x, y, z]; CCRing=CC);
evaluate_system(surface, [CC(1), CC(0), CC(0)])
```
"""
function evaluate_system(variety::AlgebraicVarietySystem, x)
    sys = variety.system
    return evaluate_H(sys, _as_acb_vector(sys.CC, x), variety.t)
end

"""
    jacobian_system(variety, x)

Evaluate the Jacobian matrix of the defining equations of `variety` at `x`.
"""
function jacobian_system(variety::AlgebraicVarietySystem, x)
    sys = variety.system
    return evaluate_Jac(sys, _as_acb_vector(sys.CC, x), variety.t)
end

function _acb_matrix(CC::AcbField, M::AbstractMatrix)
    out = Matrix{AcbFieldElem}(undef, size(M, 1), size(M, 2))
    for j in axes(M, 2), i in axes(M, 1)
        value = M[i, j]
        out[i, j] = value isa Complex ? CC(real(value), imag(value)) : CC(value)
    end
    return out
end

function _rank_from_singular_values(s::AbstractVector; rank_tol)
    return count(σ -> σ > rank_tol, s)
end

"""
    local_tangent_normal_frame(variety, x; rank_tol=1e-10) -> VarietyFrame

Compute a numerical tangent-normal frame for `variety` near `x` using an SVD of
the Jacobian.

`rank_tol` controls the singular-value threshold used to determine the local
rank and dimension.
"""
function local_tangent_normal_frame(
    variety::AlgebraicVarietySystem,
    x;
    rank_tol=1e-10,
)
    sys = variety.system
    x_mid = get_mid_vec(_as_acb_vector(sys.CC, x))
    J = jacobian_system(variety, x_mid)
    J_double = convert_to_double_matrix(J)
    fact = svd(J_double; full=true)
    rank = _rank_from_singular_values(fact.S; rank_tol=rank_tol)
    m, n = size(J)
    rank <= min(m, n) || throw(ArgumentError("computed rank is inconsistent with the Jacobian size."))

    V = fact.V
    U = fact.U
    tangent_double = V[:, rank+1:n]
    normal_double = V[:, 1:rank]
    equation_basis_double = U[:, 1:rank]

    tangent = _acb_matrix(sys.CC, tangent_double)
    normal = _acb_matrix(sys.CC, normal_double)
    equation_basis = _acb_matrix(sys.CC, equation_basis_double)

    normal_jacobian = transpose(equation_basis) * J * normal
    preconditioner = rank == 0 ? Matrix{AcbFieldElem}(undef, 0, 0) : inv_acb(normal_jacobian)
    return VarietyFrame(tangent, normal, equation_basis, preconditioner, rank, n - rank)
end

function _local_box(variety::AlgebraicVarietySystem, center, frame::VarietyFrame, tangent_radius, normal_radius)
    sys = variety.system
    u_box = _acb_unit_box_vector(sys.CC, sys.RR, frame.dim)
    w_box = _acb_unit_box_vector(sys.CC, sys.RR, frame.rank)
    center_vec = _as_acb_vector(sys.CC, center)
    tangent_part = frame.dim == 0 ? fill(sys.CC(0), length(center_vec)) : frame.tangent * (sys.CC(tangent_radius) .* u_box)
    normal_part = frame.rank == 0 ? fill(sys.CC(0), length(center_vec)) : frame.normal * (sys.CC(normal_radius) .* w_box)
    return center_vec + tangent_part + normal_part
end

function _local_box(
    variety::AlgebraicVarietySystem,
    center,
    frame::VarietyFrame,
    tangent_offsets::AbstractVector,
    tangent_radii::AbstractVector,
    normal_radius,
)
    sys = variety.system
    length(tangent_offsets) == frame.dim || throw(DimensionMismatch("tangent_offsets length must match the variety dimension."))
    length(tangent_radii) == frame.dim || throw(DimensionMismatch("tangent_radii length must match the variety dimension."))
    u_box = _acb_unit_box_vector(sys.CC, sys.RR, frame.dim)
    w_box = _acb_unit_box_vector(sys.CC, sys.RR, frame.rank)
    center_vec = _as_acb_vector(sys.CC, center)
    tangent_coords = AcbFieldElem[sys.CC(tangent_offsets[i]) + sys.CC(tangent_radii[i]) * u_box[i] for i in 1:frame.dim]
    tangent_part = frame.dim == 0 ? fill(sys.CC(0), length(center_vec)) : frame.tangent * tangent_coords
    normal_part = frame.rank == 0 ? fill(sys.CC(0), length(center_vec)) : frame.normal * (sys.CC(normal_radius) .* w_box)
    return center_vec + tangent_part + normal_part
end

function _tangent_box(variety::AlgebraicVarietySystem, center, frame::VarietyFrame, tangent_radius)
    sys = variety.system
    u_box = _acb_unit_box_vector(sys.CC, sys.RR, frame.dim)
    center_vec = _as_acb_vector(sys.CC, center)
    tangent_part = frame.dim == 0 ? fill(sys.CC(0), length(center_vec)) : frame.tangent * (sys.CC(tangent_radius) .* u_box)
    return center_vec + tangent_part
end

function _tangent_box(
    variety::AlgebraicVarietySystem,
    center,
    frame::VarietyFrame,
    tangent_offsets::AbstractVector,
    tangent_radii::AbstractVector,
)
    sys = variety.system
    length(tangent_offsets) == frame.dim || throw(DimensionMismatch("tangent_offsets length must match the variety dimension."))
    length(tangent_radii) == frame.dim || throw(DimensionMismatch("tangent_radii length must match the variety dimension."))
    u_box = _acb_unit_box_vector(sys.CC, sys.RR, frame.dim)
    center_vec = _as_acb_vector(sys.CC, center)
    tangent_coords = AcbFieldElem[sys.CC(tangent_offsets[i]) + sys.CC(tangent_radii[i]) * u_box[i] for i in 1:frame.dim]
    tangent_part = frame.dim == 0 ? fill(sys.CC(0), length(center_vec)) : frame.tangent * tangent_coords
    return center_vec + tangent_part
end

function _zero_tangent_data(frame::VarietyFrame)
    return zeros(Float64, frame.dim), zeros(Float64, frame.dim)
end

function _symmetric_tangent_data(frame::VarietyFrame, tangent_radius)
    return zeros(Float64, frame.dim), fill(Float64(tangent_radius), frame.dim)
end

function _step_tangent_data(frame::VarietyFrame, h, direction_index, direction_sign)
    h_float = Float64(h)
    if frame.dim == 0
        return Float64[], Float64[]
    elseif direction_index === nothing
        return _symmetric_tangent_data(frame, h_float)
    else
        idx = Int(direction_index)
        1 <= idx <= frame.dim || throw(ArgumentError("direction_index must be between 1 and the variety dimension."))
        offsets = zeros(Float64, frame.dim)
        radii = fill(h_float / 2, frame.dim)
        offsets[idx] = direction_sign * h_float / 2
        return offsets, radii
    end
end

function _tangent_midpoint(anchor, frame::VarietyFrame, tangent_offsets::AbstractVector)
    if frame.dim == 0
        return get_mid_vec(anchor)
    end
    return get_mid_vec(anchor + frame.tangent * AcbFieldElem[parent(anchor[1])(offset) for offset in tangent_offsets])
end

function krawczyk_test(
    variety::AlgebraicVarietySystem,
    center,
    frame::VarietyFrame,
    tangent_radius,
    normal_radius;
    rho=7/8,
)
    sys = variety.system
    frame.rank == 0 && return true, 0.0
    normal_radius > 0 || throw(ArgumentError("normal_radius must be positive."))

    tangent_offsets, tangent_radii = _symmetric_tangent_data(frame, tangent_radius)
    return krawczyk_test(variety, center, frame, tangent_offsets, tangent_radii, normal_radius; rho=rho)
end

function krawczyk_test(
    variety::AlgebraicVarietySystem,
    center,
    frame::VarietyFrame,
    tangent_offsets::AbstractVector,
    tangent_radii::AbstractVector,
    normal_radius;
    rho=7/8,
)
    sys = variety.system
    frame.rank == 0 && return true, 0.0
    normal_radius > 0 || throw(ArgumentError("normal_radius must be positive."))

    tangent_box = _tangent_box(variety, center, frame, tangent_offsets, tangent_radii)
    full_box = _local_box(variety, center, frame, tangent_offsets, tangent_radii, normal_radius)
    F_tangent = transpose(frame.equation_basis) * evaluate_system(variety, tangent_box)
    J_full = transpose(frame.equation_basis) * jacobian_system(variety, full_box)
    B = _acb_unit_box_vector(sys.CC, sys.RR, frame.rank)
    I = _acb_identity_matrix(sys.CC, frame.rank)

    term1 = -(frame.preconditioner * F_tangent) ./ sys.CC(normal_radius)
    term2 = (I - frame.preconditioner * J_full * frame.normal) * B
    K = term1 + term2
    k_norm = norm_inf(K)
    return k_norm < Float64(rho), k_norm
end

function _normal_newton_step(variety::AlgebraicVarietySystem, x, frame::VarietyFrame)
    frame.rank == 0 && return _as_acb_vector(variety.system.CC, x)
    F_compressed = transpose(frame.equation_basis) * evaluate_system(variety, x)
    delta_w = -(frame.preconditioner * F_compressed)
    return _as_acb_vector(variety.system.CC, x) + frame.normal * delta_w
end

function _compressed_residual_norm(variety::AlgebraicVarietySystem, x, frame::VarietyFrame)
    frame.rank == 0 && return 0.0
    F_compressed = transpose(frame.equation_basis) * evaluate_system(variety, x)
    return norm_inf(F_compressed)
end

function _polish_variety_center(variety::AlgebraicVarietySystem, x; rank_tol, newton_tol, max_newton_steps)
    y = get_mid_vec(_as_acb_vector(variety.system.CC, x))
    frame = local_tangent_normal_frame(variety, y; rank_tol=rank_tol)
    for _ in 1:max_newton_steps
        _compressed_residual_norm(variety, y, frame) <= newton_tol && return y, frame, true
        y_next = get_mid_vec(_normal_newton_step(variety, y, frame))
        frame_next = local_tangent_normal_frame(variety, y_next; rank_tol=rank_tol)
        if norm_inf(y_next - y) <= newton_tol
            return y_next, frame_next, true
        end
        y, frame = y_next, frame_next
    end
    return y, frame, _compressed_residual_norm(variety, y, frame) <= newton_tol
end

function _variety_distance(x::AbstractVector{AcbFieldElem}, y::AbstractVector{AcbFieldElem})
    length(x) == length(y) || throw(DimensionMismatch("Cannot compare points with different dimensions."))
    isempty(x) && return 0.0
    return maximum(abs(convert_to_double_int(x[i]) - convert_to_double_int(y[i])) for i in eachindex(x))
end

function _is_near_existing_box(x, boxes::AbstractVector{VarietyBox}; overlap_factor)
    for box in boxes
        base_radius = box.frame.dim == 0 ? box.normal_radius : box.tangent_radius
        threshold = overlap_factor * base_radius
        _variety_distance(x, box.center) <= threshold && return true
    end
    return false
end

function _neighbor_candidates(box::VarietyBox; facet_anchor_scale)
    candidates = NamedTuple{(:anchor, :direction_index, :direction_sign), Tuple{Vector{AcbFieldElem}, Union{Nothing, Int}, Int}}[]
    box.frame.dim == 0 && return candidates
    step = box.tangent_radius * facet_anchor_scale
    for i in 1:box.frame.dim
        tangent_direction = box.frame.tangent[:, i]
        push!(candidates, (anchor=get_mid_vec(box.center + step .* tangent_direction), direction_index=i, direction_sign=1))
        push!(candidates, (anchor=get_mid_vec(box.center - step .* tangent_direction), direction_index=i, direction_sign=-1))
    end
    return candidates
end

function _variety_seed_items(variety::AlgebraicVarietySystem, start)
    starts = (!isempty(start) && first(start) isa AbstractVector) ? start : (start,)
    return [
        (anchor = _as_acb_vector(variety.system.CC, s), direction_index = nothing, direction_sign = 1)
        for s in starts
    ]
end

function _process_variety_candidate(
    variety::AlgebraicVarietySystem,
    item;
    tangent_radius,
    normal_radius,
    normal_rho,
    tangent_rho,
    rank_tol,
    min_radius,
    min_tangent_radius,
    max_radius,
    max_tangent_radius,
    tangent_growth_factor,
    radius_shrink_factor,
    max_refinement_iter,
    newton_tol,
    max_newton_steps,
    facet_anchor_scale,
)
    normal_box = refine_moore_box(
        variety,
        item.anchor,
        normal_radius;
        rho=normal_rho,
        rank_tol=rank_tol,
        min_radius=min_radius,
        max_radius=max_radius,
        growth_factor=2.0,
        radius_shrink_factor=radius_shrink_factor,
        max_iter=max_refinement_iter,
        newton_tol=newton_tol,
        max_newton_steps=max_newton_steps,
    )

    normal_box.success || return (status = :rejected, box = nothing, neighbors = ())

    anchor = normal_box.center
    frame = normal_box.frame
    rn = normal_box.normal_radius
    h = Float64(tangent_radius)
    tangent_offsets, tangent_radii = _step_tangent_data(frame, h, item.direction_index, item.direction_sign)
    tangent_passed, tangent_norm = krawczyk_test(variety, anchor, frame, tangent_offsets, tangent_radii, rn; rho=tangent_rho)
    while !tangent_passed && h > min_tangent_radius
        h *= radius_shrink_factor
        tangent_offsets, tangent_radii = _step_tangent_data(frame, h, item.direction_index, item.direction_sign)
        tangent_passed, tangent_norm = krawczyk_test(variety, anchor, frame, tangent_offsets, tangent_radii, rn; rho=tangent_rho)
    end
    tangent_passed || return (status = :rejected, box = nothing, neighbors = ())

    while tangent_growth_factor * h <= max_tangent_radius
        h_next = tangent_growth_factor * h
        offsets_next, radii_next = _step_tangent_data(frame, h_next, item.direction_index, item.direction_sign)
        grown, grown_norm = krawczyk_test(variety, anchor, frame, offsets_next, radii_next, rn; rho=tangent_rho)
        grown || break
        h = h_next
        tangent_offsets, tangent_radii = offsets_next, radii_next
        tangent_norm = grown_norm
    end

    box = VarietyBox(
        _tangent_midpoint(anchor, frame, tangent_offsets),
        frame,
        isempty(tangent_radii) ? 0.0 : maximum(tangent_radii),
        rn,
        tangent_norm,
        true,
    )
    return (
        status = :accepted,
        box = box,
        neighbors = _neighbor_candidates(box; facet_anchor_scale=facet_anchor_scale),
    )
end

"""
    certified_variety_approximation(variety, start; kwargs...) -> VarietyApproximation

Build a small certified paving of an algebraic variety starting from one or more
seed points.

# Options

- `max_boxes=100`: maximum accepted boxes.
- `tangent_radius=0.1`: initial tangent radius.
- `normal_radius=tangent_radius`: initial normal radius.
- `normal_rho=1/8`: Krawczyk threshold for normal refinement.
- `tangent_rho=7/8`: Krawczyk threshold for tangent-box validation.
- `rank_tol=1e-10`: Jacobian rank threshold.
- `min_radius=1e-12`, `min_tangent_radius=min_radius`: minimum allowed radii.
- `max_radius=1.0`, `max_tangent_radius=max_radius`: maximum allowed radii.
- `tangent_growth_factor=2.0`: growth factor when expanding tangent radii.
- `radius_shrink_factor=0.5`: shrink factor after failed refinement.
- `max_refinement_iter=20`: maximum local radius adjustment attempts.
- `newton_tol=1e-24`, `max_newton_steps=10`: center polishing controls.
- `facet_anchor_scale=1.0`: scale used to generate neighboring candidate anchors.
- `overlap_factor=0.75`: distance threshold for skipping overlapping boxes.
- `threading=false`, `ntasks=Threads.nthreads()`: process batches concurrently.

# Example

```julia
@variables x y z;
CC = AcbField(128);
surface = variety_system([x^2 + y^2 + z^2 - 1], [x, y, z]; CCRing=CC);
approx = certified_variety_approximation(
    surface,
    [CC(1), CC(0), CC(0)];
    tangent_radius=1e-3,
    normal_radius=1e-3,
    max_boxes=5,
)
length(approx.boxes)
```
"""
function certified_variety_approximation(
    variety::AlgebraicVarietySystem,
    start;
    max_boxes=100,
    tangent_radius=0.1,
    normal_radius=tangent_radius,
    normal_rho=1/8,
    tangent_rho=7/8,
    rank_tol=1e-10,
    min_radius=1e-12,
    min_tangent_radius=min_radius,
    max_radius=1.0,
    max_tangent_radius=max_radius,
    tangent_growth_factor=2.0,
    radius_shrink_factor=0.5,
    max_refinement_iter=20,
    newton_tol=1e-24,
    max_newton_steps=10,
    facet_anchor_scale=1.0,
    overlap_factor=0.75,
    threading=false,
    ntasks=Base.Threads.nthreads(),
)
    max_boxes >= 0 || throw(ArgumentError("max_boxes must be nonnegative."))
    thread_count = _thread_count(threading, ntasks)
    boxes = VarietyBox[]
    queue = NamedTuple{(:anchor, :direction_index, :direction_sign), Tuple{Vector{AcbFieldElem}, Union{Nothing, Int}, Int}}[]
    append!(queue, _variety_seed_items(variety, start))
    attempted = 0
    rejected = 0
    skipped_overlap = 0

    while !isempty(queue) && length(boxes) < max_boxes
        batch_size = thread_count > 1 ? min(length(queue), max(thread_count, max_boxes - length(boxes))) : 1
        batch = [popfirst!(queue) for _ in 1:batch_size]
        attempted += length(batch)
        results = Vector{Any}(undef, length(batch))
        if thread_count > 1 && length(batch) > 1
            chunks = _thread_chunks(batch, thread_count)
            Base.Threads.@sync for chunk in chunks
                Base.Threads.@spawn begin
                    for pos in chunk
                        results[pos] = _process_variety_candidate(
                            variety,
                            batch[pos];
                            tangent_radius=tangent_radius,
                            normal_radius=normal_radius,
                            normal_rho=normal_rho,
                            tangent_rho=tangent_rho,
                            rank_tol=rank_tol,
                            min_radius=min_radius,
                            min_tangent_radius=min_tangent_radius,
                            max_radius=max_radius,
                            max_tangent_radius=max_tangent_radius,
                            tangent_growth_factor=tangent_growth_factor,
                            radius_shrink_factor=radius_shrink_factor,
                            max_refinement_iter=max_refinement_iter,
                            newton_tol=newton_tol,
                            max_newton_steps=max_newton_steps,
                            facet_anchor_scale=facet_anchor_scale,
                        )
                    end
                end
            end
        else
            results[1] = _process_variety_candidate(
                variety,
                batch[1];
                tangent_radius=tangent_radius,
                normal_radius=normal_radius,
                normal_rho=normal_rho,
                tangent_rho=tangent_rho,
                rank_tol=rank_tol,
                min_radius=min_radius,
                min_tangent_radius=min_tangent_radius,
                max_radius=max_radius,
                max_tangent_radius=max_tangent_radius,
                tangent_growth_factor=tangent_growth_factor,
                radius_shrink_factor=radius_shrink_factor,
                max_refinement_iter=max_refinement_iter,
                newton_tol=newton_tol,
                max_newton_steps=max_newton_steps,
                facet_anchor_scale=facet_anchor_scale,
            )
        end

        for result in results
            length(boxes) >= max_boxes && break
            if result.status !== :accepted
                rejected += 1
                continue
            end
            box = result.box
            if _is_near_existing_box(box.center, boxes; overlap_factor=overlap_factor)
                skipped_overlap += 1
                continue
            end

            push!(boxes, box)
            append!(queue, result.neighbors)
        end
    end

    return VarietyApproximation(boxes, attempted, rejected, skipped_overlap, max_boxes)
end

function _mid_real_float(z::AcbFieldElem)
    return Float64(Nemo.midpoint(real(z)))
end

function _obj_point(x::AbstractVector{AcbFieldElem})
    n = length(x)
    n <= 3 || throw(ArgumentError("OBJ export currently supports ambient dimension at most 3."))
    coords = [_mid_real_float(x[i]) for i in 1:n]
    while length(coords) < 3
        push!(coords, 0.0)
    end
    return coords
end

function _box_basis_and_radii(box::VarietyBox)
    basis = hcat(box.frame.tangent, box.frame.normal)
    radii = [fill(box.tangent_radius, box.frame.dim); fill(box.normal_radius, box.frame.rank)]
    return basis, radii
end

function _box_corners(box::VarietyBox)
    basis, radii = _box_basis_and_radii(box)
    n = length(box.center)
    n <= 3 || throw(ArgumentError("OBJ export currently supports ambient dimension at most 3."))
    size(basis, 2) == n || throw(DimensionMismatch("Variety box frame must span the ambient space."))
    corners = Vector{Vector{AcbFieldElem}}()
    for mask in 0:(2^n - 1)
        point = copy(box.center)
        for j in 1:n
            sign = ((mask >> (j - 1)) & 1) == 1 ? 1.0 : -1.0
            point += sign * radii[j] .* basis[:, j]
        end
        push!(corners, get_mid_vec(point))
    end
    return corners
end

function _write_obj_box!(io, box::VarietyBox, vertex_offset::Int)
    corners = _box_corners(box)
    n = length(box.center)
    for corner in corners
        x, y, z = _obj_point(corner)
        println(io, "v $x $y $z")
    end

    if n == 1
        println(io, "l $(vertex_offset + 1) $(vertex_offset + 2)")
    elseif n == 2
        println(io, "f $(vertex_offset + 1) $(vertex_offset + 2) $(vertex_offset + 4) $(vertex_offset + 3)")
    elseif n == 3
        faces = (
            (1, 2, 4, 3),
            (5, 7, 8, 6),
            (1, 5, 6, 2),
            (3, 4, 8, 7),
            (1, 3, 7, 5),
            (2, 6, 8, 4),
        )
        for face in faces
            println(io, "f ", join((vertex_offset + i for i in face), " "))
        end
    end
    return vertex_offset + length(corners)
end

"""
    export_variety_obj(boxes, filename)
    export_variety_obj(approx::VarietyApproximation, filename)
    export_variety_obj(box::VarietyBox, filename)

Export certified variety boxes to a simple Wavefront OBJ file.

Currently supports ambient dimension at most 3.
"""
function export_variety_obj(boxes::AbstractVector{VarietyBox}, filename::AbstractString)
    open(filename, "w") do io
        println(io, "# CertifiedHomotopyTracking.jl certified variety boxes")
        vertex_offset = 0
        for (i, box) in enumerate(boxes)
            println(io, "o variety_box_$i")
            vertex_offset = _write_obj_box!(io, box, vertex_offset)
        end
    end
    return filename
end

export_variety_obj(approx::VarietyApproximation, filename::AbstractString) =
    export_variety_obj(approx.boxes, filename)

export_variety_obj(box::VarietyBox, filename::AbstractString) =
    export_variety_obj([box], filename)

"""
    variety_system(F_eqs, x_vars; CCRing=AcbField(256), projective=false) -> AlgebraicVarietySystem

Compile symbolic equations as a static algebraic variety system.

Complex numeric coefficients are supported and are internally parameterized as
fixed constants.

# Options

- `CCRing=AcbField(256)`: ACB field used by the underlying system.
- `projective=false`: compile a projective/homogeneous version.

# Example

```julia
CC = AcbField(128);
@variables x y z;
surface = variety_system([x^2 + y^2 + z^2 - 1], [x, y, z]; CCRing=CC);

p = [CC(1), CC(0), CC(0)];
evaluate_system(surface, p)
jacobian_system(surface, p)
```
"""
function variety_system(
    F_eqs::AbstractVector{<:Union{Num, Complex{Num}}},
    x_vars::AbstractVector{Num};
    CCRing=AcbField(256),
    projective=false,
    patch_vector=nothing,
    patch_rng=nothing,
)
    coeff_vars = Num[]
    coeff_values = Any[]
    coeff_map = Dict{Any, Num}()
    F_param = _parameterize_complex_coefficients!(F_eqs, x_vars, coeff_vars, coeff_values, coeff_map)
    compiled = compile_system(
        F_param,
        x_vars;
        projective=projective,
        patch_vector=patch_vector,
        patch_rng=patch_rng,
        const_vars=coeff_vars,
    )
    coeff_vals = [_coefficient_value(CCRing, coeff) for coeff in coeff_values]
    sys = isempty(coeff_vals) ? SpecializedHomotopy(compiled, CCRing) : SpecializedHomotopy(compiled, AcbFieldElem[], AcbFieldElem[], coeff_vals)
    return AlgebraicVarietySystem(sys)
end
