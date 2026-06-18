import HomotopyContinuation

const HCTracePoint = NamedTuple{(:t, :x), Tuple{Float64, Vector{ComplexF64}}}

_hc_trace_time(t) = Float64(real(t))
_hc_trace_x(x) = ComplexF64.(copy(x))

function _hc_tracker_steps(tracker, field::Symbol)
    tracker_state = HomotopyContinuation.state(tracker)
    hasproperty(tracker_state, field) || return missing
    return getproperty(tracker_state, field)
end

"""
    collect_hc_trace(H, x_start; kwargs...)

Collect a numerical trace along one HomotopyContinuation.jl homotopy path using the
low-level `Tracker` iterator API.

The returned points are the numerical approximations produced by HC.jl at the start
and at accepted tracker steps. They are **not certified** points, and this function
does not perform interval, Krawczyk, or a posteriori certification.

Keyword arguments are passed into `HomotopyContinuation.TrackerOptions` where
applicable:

* `t_start = 1.0`, `t_target = 0.0`: tracking direction.
* `parameters = :default`
* `max_steps = 10_000`
* `max_step_size = Inf`
* `throw_on_failure = false`: throw an `ErrorException` if HC.jl does not report
  success.

Returns a named tuple
`(trace, status, success, accepted_steps, rejected_steps)`, where `trace` is a
vector of named tuples `(t = ..., x = ...)`.
"""
function collect_hc_trace(
    H,
    x_start;
    t_start = 1.0,
    t_target = 0.0,
    parameters = :default,
    max_steps = 10_000,
    max_step_size = Inf,
    throw_on_failure = false,
)
    options = HomotopyContinuation.TrackerOptions(;
        automatic_differentiation = 1,
        max_steps = max_steps,
        max_step_size = Float64(max_step_size),
        max_initial_step_size = Inf,
        extended_precision = true,
        parameters = parameters,
    )
    tracker = HomotopyContinuation.Tracker(H; options = options)

    trace = HCTracePoint[]

    for (i, (x, t)) in enumerate(
        HomotopyContinuation.iterator(tracker, x_start, t_start, t_target)
    )
        push!(trace, (t = _hc_trace_time(t), x = _hc_trace_x(x)))
    end

    code = HomotopyContinuation.status(tracker)
    success = HomotopyContinuation.is_success(code)
    if throw_on_failure && !success
        error("HomotopyContinuation tracker failed with status $(code)")
    end

    return (
        trace = trace,
        status = code,
        success = success,
        accepted_steps = _hc_tracker_steps(tracker, :accepted_steps),
        rejected_steps = _hc_tracker_steps(tracker, :rejected_steps),
    )
end
