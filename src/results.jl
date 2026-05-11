export TrackResult, solution, certified_region, approximate_solution, succeeded,
    projective_solution, near_infinity, input_start, refined_start,
    projective_input_start, projective_refined_start

struct TrackResult
    root::Vector{AcbFieldElem}
    projective_root::Vector{AcbFieldElem}
    input_start::Vector{AcbFieldElem}
    refined_start::Vector{AcbFieldElem}
    projective_input_start::Vector{AcbFieldElem}
    projective_refined_start::Vector{AcbFieldElem}
    success::Bool
    status::Symbol
    iterations::Int
    accepted_steps::Int
    rejected_steps::Int
    final_t::Float64
    final_h::Float64
    final_radius::Float64
    final_krawczyk_norm::Float64
    rho::Float64
    patch_idx::Int
    near_infinity::Bool
    message::String
end

TrackResult(
    root::Vector{AcbFieldElem},
    success::Bool,
    status::Symbol,
    iterations::Int,
    accepted_steps::Int,
    rejected_steps::Int,
    final_t::Float64,
    final_h::Float64,
    final_radius::Float64,
    final_krawczyk_norm::Float64,
    rho::Float64,
    patch_idx::Int,
    message::String,
) = TrackResult(
    root,
    AcbFieldElem[],
    AcbFieldElem[],
    AcbFieldElem[],
    AcbFieldElem[],
    AcbFieldElem[],
    success,
    status,
    iterations,
    accepted_steps,
    rejected_steps,
    final_t,
    final_h,
    final_radius,
    final_krawczyk_norm,
    rho,
    patch_idx,
    false,
    message,
)

certified_region(res::TrackResult) = res.root
solution(res::TrackResult) = convert_to_double_int.(res.root)
succeeded(res::TrackResult) = res.success
projective_solution(res::TrackResult) = isempty(res.projective_root) ? res.root : res.projective_root
near_infinity(res::TrackResult) = res.near_infinity
input_start(res::TrackResult) = res.input_start
refined_start(res::TrackResult) = res.refined_start
projective_input_start(res::TrackResult) = res.projective_input_start
projective_refined_start(res::TrackResult) = res.projective_refined_start

function approximate_solution(res::TrackResult; digits=8)
    return [round(z; digits=digits) for z in solution(res)]
end

function Base.iterate(res::TrackResult, state=1)
    state == 1 && return (res.root, 2)
    state == 2 && return (res.success, 3)
    return nothing
end

Base.length(::TrackResult) = 2

function Base.show(io::IO, res::TrackResult)
    print(io, "TrackResult(")
    print(io, res.success ? "success" : "failure")
    print(io, ", status=", res.status)
    print(io, ", iterations=", res.iterations)
    print(io, ", accepted=", res.accepted_steps)
    print(io, ", rejected=", res.rejected_steps)
    print(io, ", final_radius=", res.final_radius)
    if res.near_infinity
        print(io, ", near_infinity=true")
    end
    if isfinite(res.final_krawczyk_norm)
        print(io, ", final_krawczyk_norm=", res.final_krawczyk_norm)
    end
    print(io, ", approx=", approximate_solution(res))
    print(io, ")")
end

function _projective_to_affine(sys::HCSystem, X; affine_chart_atol=1e-10)
    X0 = X[1]
    if mag_complex(X0) <= affine_chart_atol
        return copy(X), true
    end

    affine = Vector{AcbFieldElem}(undef, length(X) - 1)
    for i in eachindex(affine)
        affine[i] = X[i + 1] / X0
    end
    return affine, false
end

function _coordinates_for_result(sys::HCSystem, x; affine_chart_atol=1e-10)
    if has_projective_patch(sys)
        affine, is_near_infinity = _projective_to_affine(sys, x; affine_chart_atol=affine_chart_atol)
        return affine, is_near_infinity
    end
    return copy(x), false
end

function _track_result(
    sys::HCSystem,
    x,
    success,
    status,
    iterations,
    accepted_steps,
    rejected_steps,
    t,
    h,
    r,
    k_norm,
    rho,
    message;
    affine_chart_atol=1e-10,
    input_start=AcbFieldElem[],
    tracking_input_start=AcbFieldElem[],
    tracking_refined_start=AcbFieldElem[],
)
    root = copy(x)
    projective_root = AcbFieldElem[]
    result_input_start = copy(input_start)
    result_refined_start = copy(tracking_refined_start)
    result_projective_input_start = AcbFieldElem[]
    result_projective_refined_start = AcbFieldElem[]
    result_success = success
    result_status = status
    result_message = message
    result_near_infinity = false

    if has_projective_patch(sys)
        projective_root = copy(x)
        result_projective_input_start = copy(tracking_input_start)
        result_projective_refined_start = copy(tracking_refined_start)
        if !isempty(tracking_refined_start)
            result_refined_start, _ = _coordinates_for_result(sys, tracking_refined_start; affine_chart_atol=affine_chart_atol)
        end
        if success
            root, result_near_infinity = _projective_to_affine(sys, x; affine_chart_atol=affine_chart_atol)
            if result_near_infinity
                result_success = false
                result_status = :near_infinity
                result_message = "Path tracked in projective coordinates, but the result is near infinity in the original affine chart."
            end
        end
    end

    return TrackResult(
        root,
        projective_root,
        result_input_start,
        result_refined_start,
        result_projective_input_start,
        result_projective_refined_start,
        result_success,
        result_status,
        iterations,
        accepted_steps,
        rejected_steps,
        Float64(t),
        Float64(h),
        Float64(r),
        Float64(k_norm),
        Float64(rho),
        sys.patch_idx,
        result_near_infinity,
        result_message,
    )
end
