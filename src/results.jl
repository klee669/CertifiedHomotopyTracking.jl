export TrackResult, solution, certified_region, approximate_solution, succeeded

struct TrackResult
    root::Vector{AcbFieldElem}
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
    message::String
end

certified_region(res::TrackResult) = res.root
solution(res::TrackResult) = convert_to_double_int.(res.root)
succeeded(res::TrackResult) = res.success

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
    if isfinite(res.final_krawczyk_norm)
        print(io, ", final_krawczyk_norm=", res.final_krawczyk_norm)
    end
    print(io, ", approx=", approximate_solution(res))
    print(io, ")")
end

function _track_result(sys::HCSystem, x, success, status, iterations, accepted_steps, rejected_steps, t, h, r, k_norm, rho, message)
    return TrackResult(
        copy(x),
        success,
        status,
        iterations,
        accepted_steps,
        rejected_steps,
        Float64(t),
        Float64(h),
        Float64(r),
        Float64(k_norm),
        Float64(rho),
        sys.patch_idx,
        message,
    )
end

