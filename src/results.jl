export TrackResult, solution, certified_region, approximate_solution,
    projective_solution, near_infinity, input_start, refined_start,
    projective_input_start, projective_refined_start, path_boxes

import Base: success

"""
    TrackResult

Result object returned by [`track_path`](@ref).

The most commonly used fields are:

- `success::Bool`: whether the path was certified and finished in the requested affine chart.
- `status::Symbol`: machine-readable status such as `:success`, `:step_too_small`,
  `:near_infinity`, or `:final_refinement_failed`.
- `root`: certified endpoint region in affine coordinates when possible.
- `projective_root`: certified endpoint region in homogeneous coordinates for projective tracking.
- `iterations`, `accepted_steps`, `rejected_steps`: path-tracking counters.
- `final_t`, `final_h`, `final_radius`, `final_krawczyk_norm`: final tracking diagnostics. For inspecting the last attempted step when tracking fails.
- `initial_precision`, `final_precision`: precision used by adaptive tracking.

`TrackResult` also iterates as `(certified_region, success)` for compatibility
with older code, but new code should prefer the named accessors.

# Example

```julia
@variables x y;
CC = AcbField(256);
F = [x^2 + 3y - 4, y^2 + 3];
G = [x^2 - 1, y^2 - 1];
H = straight_line_homotopy(F, G, [x, y]; CCRing=CC);

res = track_path(H, [CC(1), CC(-1)]);
success(res)
solution(res)
certified_region(res)
```
"""
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
    initial_precision::Int
    final_precision::Int
    message::String
    boxes::Vector{Any}
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
    0,
    0,
    message,
    Any[],
)

"""
    certified_region(res::TrackResult)

Return the certified endpoint center stored in `res.root`.

The certification radius is stored separately as `res.final_radius`; the entries
returned here may have zero Arb/ACB radius. For projective tracking this is the
endpoint converted back to the original affine chart when possible. If the path
ends near infinity in that affine chart, inspect [`projective_solution`](@ref).
"""
certified_region(res::TrackResult) = res.root

"""
    solution(res::TrackResult)

Return a `Vector{ComplexF64}` obtained from the midpoint of
[`certified_region`](@ref). This is convenient for display and downstream
numerical checks, but it is not itself an interval certificate.
"""
solution(res::TrackResult) = convert_to_double_int.(res.root)

"""
    success(res::TrackResult)

Return `res.success`.

Use this before trusting the endpoint as a completed certified path. For failure
diagnostics inspect `res.status`.
"""
success(res::TrackResult) = res.success

"""
    succeeded(res::TrackResult)

Compatibility alias for [`success`](@ref).
"""
succeeded(res::TrackResult) = success(res)

"""
    projective_solution(res::TrackResult)

Return the certified endpoint in projective coordinates when available.

For affine tracking, this falls back
to `res.root`.
"""
projective_solution(res::TrackResult) = isempty(res.projective_root) ? res.root : res.projective_root

"""
    near_infinity(res::TrackResult)

Return whether a projectively tracked endpoint is near infinity in the original
affine chart.
"""
near_infinity(res::TrackResult) = res.near_infinity

"""
    input_start(res::TrackResult)

Return the start point supplied by the user, converted to the tracker's field.
"""
input_start(res::TrackResult) = res.input_start

"""
    refined_start(res::TrackResult)

Return the start point after initial certified Moore-box refinement.
"""
refined_start(res::TrackResult) = res.refined_start

"""
    projective_input_start(res::TrackResult)

Return the user-supplied start point represented in projective coordinates when
projective tracking is used.
"""
projective_input_start(res::TrackResult) = res.projective_input_start

"""
    projective_refined_start(res::TrackResult)

Return the initially refined start point represented in projective coordinates
when projective tracking is used.
"""
projective_refined_start(res::TrackResult) = res.projective_refined_start
path_boxes(res::TrackResult) = res.boxes

"""
    approximate_solution(res::TrackResult; digits=8)

Return `solution(res)` rounded to `digits` decimal digits.
"""
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
    if res.initial_precision > 0
        print(io, ", precision=", res.initial_precision, "->", res.final_precision)
    end
    if isfinite(res.final_krawczyk_norm)
        print(io, ", final_krawczyk_norm=", res.final_krawczyk_norm)
    end
    print(io, ", approx=", approximate_solution(res))
    print(io, ")")
end

function _result_with_field(res::TrackResult, CC::AcbField)
    isempty(res.root) || parent(res.root[1]) !== CC || return res
    convert_vec(values) = _convert_acb_vector(CC, values)
    return TrackResult(
        convert_vec(res.root),
        convert_vec(res.projective_root),
        convert_vec(res.input_start),
        convert_vec(res.refined_start),
        convert_vec(res.projective_input_start),
        convert_vec(res.projective_refined_start),
        res.success,
        res.status,
        res.iterations,
        res.accepted_steps,
        res.rejected_steps,
        res.final_t,
        res.final_h,
        res.final_radius,
        res.final_krawczyk_norm,
        res.rho,
        res.patch_idx,
        res.near_infinity,
        res.initial_precision,
        res.final_precision,
        res.message,
        Any[_path_box_with_field(box, CC) for box in res.boxes],
    )
end

function _chart_solution_to_projective(sys::SpecializedHomotopy, x; patch_idx::Int=sys.patch_idx)
    X = Vector{eltype(x)}(undef, sys.compiled.n_vars)
    for i in eachindex(X)
        if i == patch_idx
            X[i] = one(x[1])
        elseif i < patch_idx
            X[i] = x[i]
        else
            X[i] = x[i - 1]
        end
    end
    return X
end

function _projective_to_affine(sys::SpecializedHomotopy, x; affine_chart_atol=1e-10)
    X = uses_projective_charts(sys) ? _chart_solution_to_projective(sys, x) : x
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

function _coordinates_for_result(sys::SpecializedHomotopy, x; affine_chart_atol=1e-10)
    if uses_projective_charts(sys) || has_projective_patch(sys)
        affine, is_near_infinity = _projective_to_affine(sys, x; affine_chart_atol=affine_chart_atol)
        return affine, is_near_infinity
    end
    return copy(x), false
end

function _track_result(
    sys::SpecializedHomotopy,
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
    message::AbstractString;
    affine_chart_atol=1e-10,
    input_start=AcbFieldElem[],
    tracking_input_start=AcbFieldElem[],
    tracking_refined_start=AcbFieldElem[],
    tracking_start_patch_idx=sys.patch_idx,
    initial_precision=precision(sys.CC),
    boxes=Any[],
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

    if uses_projective_charts(sys) || has_projective_patch(sys)
        projective_root = uses_projective_charts(sys) ? _chart_solution_to_projective(sys, x) : copy(x)
        result_projective_input_start = uses_projective_charts(sys) && !isempty(tracking_input_start) ?
            _chart_solution_to_projective(sys, tracking_input_start; patch_idx=tracking_start_patch_idx) :
            copy(tracking_input_start)
        result_projective_refined_start = uses_projective_charts(sys) && !isempty(tracking_refined_start) ?
            _chart_solution_to_projective(sys, tracking_refined_start; patch_idx=tracking_start_patch_idx) :
            copy(tracking_refined_start)
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
        initial_precision,
        precision(sys.CC),
        result_message,
        boxes,
    )
end
