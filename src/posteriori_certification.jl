# A posteriori path certification core.
#
# This file owns the trace-certification machinery:
# - Certified trace nodes and diagnostics
# - adaptive local-parameter certification, including endpoint-box and Hermite-tube variants
# - recursive subdivision helpers and progress reporting

struct CertifiedTraceNode
    t::Float64
    x::Vector{AcbFieldElem}
    radius::Float64
    A::Matrix{AcbFieldElem}
    v::Vector{AcbFieldElem}
end

mutable struct TraceCertDiagnostics
    mode::Symbol
    node_refinements::Int
    midpoint_refinements::Int
    segment_attempts::Int
    krawczyk_validation_calls::Int
    velocity_computations::Int
    preconditioner_computations::Int
    local_parameter_choices::Dict{Symbol,Int}
    time_in_refinement::Float64
    time_in_validation::Float64
    time_in_local_parameter_search::Float64
    max_Y::Float64
    max_Z::Float64
    max_Y_over_r::Float64
end

function _trace_cert_diagnostics(mode::Symbol = :off)
    mode in (:off, :basic, :timing) ||
        throw(ArgumentError("diagnostics must be one of :off, :basic, :timing."))
    return TraceCertDiagnostics(
        mode,
        0,
        0,
        0,
        0,
        0,
        0,
        Dict{Symbol,Int}(),
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )
end

_diag_basic(d::TraceCertDiagnostics) = d.mode != :off
_diag_timing(d::TraceCertDiagnostics) = d.mode == :timing
_diag_basic(::Nothing) = false
_diag_timing(::Nothing) = false

function _freeze_trace_cert_diagnostics(diagnostics)
    snapshot = (
        node_refinements = diagnostics.node_refinements,
        midpoint_refinements = diagnostics.midpoint_refinements,
        segment_attempts = diagnostics.segment_attempts,
        krawczyk_validation_calls = diagnostics.krawczyk_validation_calls,
        velocity_computations = diagnostics.velocity_computations,
        preconditioner_computations = diagnostics.preconditioner_computations,
        local_parameter_choices = copy(diagnostics.local_parameter_choices),
        max_Y = diagnostics.max_Y,
        max_Z = diagnostics.max_Z,
        max_Y_over_r = diagnostics.max_Y_over_r,
    )
    _diag_timing(diagnostics) || return snapshot
    return merge(
        snapshot,
        (
            time_in_refinement = diagnostics.time_in_refinement,
            time_in_validation = diagnostics.time_in_validation,
            time_in_local_parameter_search = diagnostics.time_in_local_parameter_search,
        ),
    )
end

function _record_local_parameter_choice!(diagnostics::TraceCertDiagnostics, local_parameter)
    _diag_basic(diagnostics) || return nothing
    name = local_parameter.name
    diagnostics.local_parameter_choices[name] = get(diagnostics.local_parameter_choices, name, 0) + 1
    return nothing
end
_record_local_parameter_choice!(::Nothing, local_parameter) = nothing

function _record_yz_max!(diagnostics::TraceCertDiagnostics, Y, Z, Y_over_r)
    _diag_basic(diagnostics) || return nothing
    diagnostics.max_Y = max(diagnostics.max_Y, Float64(Y))
    diagnostics.max_Z = max(diagnostics.max_Z, Float64(Z))
    diagnostics.max_Y_over_r = max(diagnostics.max_Y_over_r, Float64(Y_over_r))
    return nothing
end
_record_yz_max!(::Nothing, Y, Z, Y_over_r) = nothing

# converting HC.jl trace points to the format needed for CHT.jl certification
# converting number types
# reverting t values from [1,0] back to [0,1]
function _trace_point_for_cht(sys::HCSystem, point; time_map)
    return (t = Float64(time_map(point.t)), x = sys.CC.(point.x))
end

function _prepare_trace_points(sys::HCSystem, trace; time_map)
    points = [_trace_point_for_cht(sys, point; time_map = time_map) for point in trace]
    length(points) >= 2 || throw(ArgumentError("trace must contain at least two points."))

    if points[end].t < points[1].t
        reverse!(points)
    end

    for i in 1:length(points)-1
        points[i].t < points[i+1].t ||
            throw(ArgumentError("trace times must be strictly monotone after applying time_map."))
    end
    return points
end

function _compute_preconditioner_counted(sys::HCSystem, x, t, diagnostics)
    _diag_basic(diagnostics) && (diagnostics.preconditioner_computations += 1)
    return compute_preconditioner(sys, x, t)
end

function _compute_velocity_counted(sys::HCSystem, x, t, A, diagnostics)
    _diag_basic(diagnostics) && (diagnostics.velocity_computations += 1)
    return compute_velocity(sys, x, t, A)
end

function _certify_trace_node(
    sys::HCSystem,
    point;
    node_refinement_radius,
    tau,
    diagnostics,
    refinement_counter = :node_refinements,
)
    A0 = _compute_preconditioner_counted(sys, point.x, sys.CC(point.t), diagnostics)
    if _diag_basic(diagnostics)
        setproperty!(diagnostics, refinement_counter, getproperty(diagnostics, refinement_counter) + 1)
    end
    x_refined = r_refined = A_refined = nothing
    success = false
    if _diag_timing(diagnostics)
        elapsed = @elapsed begin
            x_refined, r_refined, A_refined, success =
                refine_moore_box(sys, point.x, sys.CC(point.t), node_refinement_radius, A0; tau = tau)
        end
        diagnostics.time_in_refinement += elapsed
    else
        x_refined, r_refined, A_refined, success =
            refine_moore_box(sys, point.x, sys.CC(point.t), node_refinement_radius, A0; tau = tau)
    end
    success || return nothing

    v = _compute_velocity_counted(sys, x_refined, sys.CC(point.t), A_refined, diagnostics)
    return CertifiedTraceNode(point.t, x_refined, Float64(r_refined), A_refined, v)
end

function _prepare_certified_trace_nodes(
    sys::HCSystem,
    points;
    node_refinement_radius,
    tau,
    diagnostics,
)
    nodes = CertifiedTraceNode[]
    for point in points
        node = _certify_trace_node(
            sys,
            point;
            node_refinement_radius = node_refinement_radius,
            tau = tau,
            diagnostics = diagnostics,
            refinement_counter = :node_refinements,
        )
        node === nothing &&
            throw(ArgumentError("failed to refine trace node at t=$(point.t)"))
        push!(nodes, node)
    end
    return nodes
end

function _node_refinement_failure_result(points, err, method, diagnostics_data)
    return (
        success = false,
        status = :node_refinement_failed,
        boxes = NamedTuple[],
        segments = NamedTuple[],
        failed_segments = [(status = :node_refinement_failed, error = err)],
        total_boxes = 0,
        original_segments = max(length(points) - 1, 0),
        max_depth = 0,
        max_krawczyk_norm = Inf,
        method = method,
        boxes_by_parent = Dict{Any,Int}(),
        boxes_by_depth = Dict{Any,Int}(),
        diagnostics = _freeze_trace_cert_diagnostics(diagnostics_data),
    )
end

function _segment_radius_candidates(base_radius, growth, max_radius)
    candidates = Float64[]
    for factor in growth
        r = Float64(base_radius) * Float64(factor)
        if r <= max_radius
            push!(candidates, r)
        end
    end
    isempty(candidates) && push!(candidates, Float64(base_radius))
    return unique(candidates)
end

function _count_by_field(items, field::Symbol)
    counts = Dict{Any,Int}()
    for item in items
        key = getproperty(item, field)
        counts[key] = get(counts, key, 0) + 1
    end
    return counts
end

"""
    _certify_hc_path_a_posteriori(H_hc, sys, x_start; kwargs...)

Internal helper: collect an HC.jl low-level numerical trace and certify it a posteriori with CHT.
If certification fails, the trace is recollected with the next value in
`max_step_size_schedule`.

By default this uses HC.jl's adaptive step choice first:

```julia
max_step_size_schedule = (Inf, 0.05, 0.02, 0.01)
```

The supported certification mode is `:adaptive_local_parameter`.
"""
function _certify_hc_path_a_posteriori(
    H_hc,
    sys::HCSystem,
    x_start;
    time_map = identity,
    certification = :adaptive_local_parameter,
    max_step_size_schedule = (Inf, 0.05, 0.02, 0.01),
    parameters = :default,
    max_steps = 10_000,
    throw_on_hc_failure = false,
    max_depth_schedule = nothing,
    certification_kwargs...,
)
    isempty(max_step_size_schedule) &&
        throw(ArgumentError("max_step_size_schedule must contain at least one value."))
    attempts = NamedTuple[]
    last_cert = nothing
    last_trace = nothing
    source = sys.source === nothing ? sys.compiled.source : sys.source
    source === nothing &&
        throw(ArgumentError("a posteriori HC tracing requires homotopy source metadata."))
    t_start = _posteriori_default_t_start(source)
    t_target = _posteriori_default_t_target(source)

    for max_step_size in max_step_size_schedule
        trace_out = collect_hc_trace(
            H_hc,
            x_start;
            t_start = t_start,
            t_target = t_target,
            parameters = parameters,
            max_steps = max_steps,
            max_step_size = max_step_size,
            throw_on_failure = throw_on_hc_failure,
        )
        last_trace = trace_out

        if !trace_out.success
            push!(
                attempts,
                (
                    max_step_size = max_step_size,
                    hc_success = false,
                    hc_status = trace_out.status,
                    trace_points = length(trace_out.trace),
                    cert_success = false,
                    total_boxes = 0,
                    max_depth = 0,
                ),
            )
            continue
        end

        certification_options = NamedTuple(certification_kwargs)
        depth_values = max_depth_schedule === nothing ? (get(certification_options, :max_depth, nothing),) : max_depth_schedule
        for depth_value in depth_values
            opts = depth_value === nothing ? certification_options : merge(certification_options, (; max_depth = depth_value))
            cert = try
                certification == :adaptive_local_parameter ||
                    throw(ArgumentError("unsupported certification mode $(certification)"))
                certify_hc_trace_adaptive_local_parameter(
                    sys,
                    trace_out.trace;
                    time_map = time_map,
                    opts...,
                )
            catch err
                err isa InterruptException && rethrow(err)
                (
                    success = false,
                    status = :certification_error,
                    error = err,
                    segments = NamedTuple[],
                    boxes = NamedTuple[],
                    failed_segments = NamedTuple[],
                    total_boxes = 0,
                    original_segments = max(length(trace_out.trace) - 1, 0),
                    max_depth = depth_value === nothing ? 0 : depth_value,
                    max_krawczyk_norm = Inf,
                    method = certification,
                    boxes_by_parent = Dict{Any,Int}(),
                    boxes_by_depth = Dict{Any,Int}(),
                    diagnostics = missing,
                )
            end
            last_cert = cert
            push!(
                attempts,
                (
                    max_step_size = max_step_size,
                    hc_success = true,
                    hc_status = trace_out.status,
                    trace_points = length(trace_out.trace),
                    cert_success = cert.success,
                    cert_status = haskey(cert, :status) ? cert.status : missing,
                    total_boxes = cert.total_boxes,
                    max_depth = cert.max_depth,
                ),
            )
            cert.success && return merge(
                cert,
                (
                    hc_trace = trace_out,
                    max_step_size = max_step_size,
                    attempts = attempts,
                    certification = certification,
                ),
            )
        end
    end

    if last_cert === nothing
        return (
            success = false,
            hc_trace = last_trace,
            max_step_size = last(attempts).max_step_size,
            attempts = attempts,
            certification = certification,
            status = :hc_failed,
            total_boxes = 0,
            original_segments = 0,
            max_depth = 0,
            failed_segments = NamedTuple[],
        )
    end

    return merge(
        last_cert,
        (
            hc_trace = last_trace,
            max_step_size = last(attempts).max_step_size,
            attempts = attempts,
            certification = certification,
        ),
    )
end


# ------------------------------------------------------------------------------
# Adaptive local-parameter certification
# ------------------------------------------------------------------------------

export certify_hc_trace_adaptive_local_parameter

function _local_parameter_show_progress(i, total, boxes, failures, max_depth; label, attempts = nothing, final = false)
    pct = total == 0 ? 100.0 : 100 * i / total
    box_depth = isempty(boxes) ? 0 : maximum(b -> haskey(b, :depth) ? b.depth : 0, boxes)
    failure_depth = isempty(failures) ? 0 : maximum(f -> haskey(f, :depth) ? f.depth : 0, failures)
    msg = "certifying trace segment $(i)/$(total) ($(round(pct; digits = 1))%) | " *
          "boxes=$(length(boxes)) failed=$(length(failures)) depth=$(max(max_depth, box_depth, failure_depth))"
    attempts === nothing || (msg *= " attempts=$(attempts)")
    if final
        println("\r", msg, "\033[K")
    else
        print("\r", msg, "\033[K")
    end
    flush(stdout)
end

function _local_parameter_progress_state(i, total)
    return (i = i, total = total, attempts = Ref(0), max_depth = Ref(0), last_ns = Ref(UInt64(0)))
end

function _local_parameter_maybe_show_progress!(state, boxes, failures, depth; force = false, final = false)
    state === nothing && return
    state.attempts[] += 1
    state.max_depth[] = max(state.max_depth[], depth)
    now = time_ns()
    if force || now - state.last_ns[] >= UInt64(500_000_000)
        _local_parameter_show_progress(
            state.i,
            state.total,
            boxes,
            failures,
            state.max_depth[];
            label = "local_parameter",
            attempts = state.attempts[],
            final = final,
        )
        state.last_ns[] = now
    end
end

_local_parameter_mid_float(a::AcbFieldElem) = Float64(Nemo.midpoint(real(a)))
_local_parameter_mid_float(a) = Float64(Nemo.midpoint(a))
_local_parameter_real_acb(sys::HCSystem, a) = sys.CC(a, sys.RR(0))
_local_parameter_real_part(sys::HCSystem, z) = _local_parameter_real_acb(sys, real(z))
_local_parameter_imag_part(sys::HCSystem, z) = _local_parameter_real_acb(sys, imag(z))

function _local_parameter_tm_part(sys::HCSystem, tm::TaylorModel3, part::Symbol)
    f = part === :real ? real : imag
    return TaylorModel3(
        _local_parameter_real_acb(sys, f(tm.c0)),
        _local_parameter_real_acb(sys, f(tm.c1)),
        _local_parameter_real_acb(sys, f(tm.c2)),
        _local_parameter_real_acb(sys, f(tm.c3)),
        _local_parameter_real_acb(sys, f(tm.rem)),
        tm.h,
    )
end

function _local_parameter_real_rows(sys::HCSystem, Fval)
    rows = Vector{AcbFieldElem}(undef, 2 * length(Fval))
    for i in eachindex(Fval)
        rows[2i - 1] = _local_parameter_real_part(sys, Fval[i])
        rows[2i] = _local_parameter_imag_part(sys, Fval[i])
    end
    return rows
end

function _local_parameter_real_rows_tm(sys::HCSystem, F_tm)
    rows = Vector{TaylorModel3{AcbFieldElem,ArbFieldElem}}(undef, 2 * length(F_tm))
    for i in eachindex(F_tm)
        rows[2i - 1] = _local_parameter_tm_part(sys, F_tm[i], :real)
        rows[2i] = _local_parameter_tm_part(sys, F_tm[i], :imag)
    end
    return rows
end

function _local_parameter_fill_real_column!(Jr, col, row0, D, imag_coordinate, sys)
    if imag_coordinate
        Jr[row0, col] = -_local_parameter_imag_part(sys, D)
        Jr[row0 + 1, col] = _local_parameter_real_part(sys, D)
    else
        Jr[row0, col] = _local_parameter_real_part(sys, D)
        Jr[row0 + 1, col] = _local_parameter_imag_part(sys, D)
    end
    return nothing
end

function _local_parameter_name(local_parameter_index, n)
    local_parameter_index == 2n + 1 && return :t
    var_idx = (local_parameter_index + 1) ÷ 2
    part = isodd(local_parameter_index) ? "re" : "im"
    return Symbol("x$(var_idx)_$(part)")
end

function _local_parameter_name(sys::HCSystem, local_parameter_index, n)
    local_parameter_index == 2n + 1 && return :t
    source = sys.source === nothing ? sys.compiled.source : sys.source
    source === nothing && return _local_parameter_name(local_parameter_index, n)
    vars = collect(source.variables)
    var_idx = (local_parameter_index + 1) ÷ 2
    var_idx <= length(vars) || return _local_parameter_name(local_parameter_index, n)
    part = isodd(local_parameter_index) ? "re" : "im"
    return Symbol(string(Symbol(vars[var_idx])), "_", part)
end

function _local_parameter_value(sys::HCSystem, x, t, local_parameter_index)
    n = length(x)
    local_parameter_index == 2n + 1 && return _local_parameter_real_part(sys, t)
    var_idx = (local_parameter_index + 1) ÷ 2
    isodd(local_parameter_index) && return _local_parameter_real_part(sys, x[var_idx])
    return _local_parameter_imag_part(sys, x[var_idx])
end

function _local_parameter_unknown_endpoint(sys::HCSystem, x, t, local_parameter_index)
    values = AcbFieldElem[]
    for i in eachindex(x)
        re = _local_parameter_real_part(sys, x[i])
        im = _local_parameter_imag_part(sys, x[i])
        2i - 1 == local_parameter_index || push!(values, re)
        2i == local_parameter_index || push!(values, im)
    end
    local_parameter_index == 2 * length(x) + 1 || push!(values, _local_parameter_real_part(sys, t))
    return values
end

function _local_parameter_reconstruct_x_t(sys::HCSystem, u, local_parameter_index, local_parameter_value, n::Integer)
    coords = Vector{typeof(local_parameter_value)}(undef, 2n)
    k = 1
    if local_parameter_index == 2n + 1
        for j in 1:2n
            coords[j] = u[k]
            k += 1
        end
        t = local_parameter_value
    else
        for j in 1:2n
            if j == local_parameter_index
                coords[j] = local_parameter_value
            else
                coords[j] = u[k]
                k += 1
            end
        end
        t = u[end]
    end
    im_unit = sys.CC(0, 1)
    x = [coords[2i - 1] + im_unit * coords[2i] for i in 1:n]
    return x, t
end

function _local_parameter_candidates(sys::HCSystem, node0, node1; min_local_parameter_motion = 1e-14)
    n = length(node0.x)
    candidates = NamedTuple[]
    for i in 1:n
        for (idx, v0, v1) in (
            (2i - 1, real(node0.x[i]), real(node1.x[i])),
            (2i, imag(node0.x[i]), imag(node1.x[i])),
        )
            a0 = _local_parameter_mid_float(v0)
            a1 = _local_parameter_mid_float(v1)
            delta = a1 - a0
            scale = max(1.0, abs(a0), abs(a1))
            abs(delta) > min_local_parameter_motion || continue
            push!(
                candidates,
                (
                    index = idx,
                    name = _local_parameter_name(sys, idx, n),
                    score = abs(delta) / scale,
                    delta = delta,
                    valid = true,
                ),
            )
        end
    end
    t_delta = node1.t - node0.t
    abs(t_delta) > min_local_parameter_motion && push!(
        candidates,
        (index = 2n + 1, name = :t, score = abs(t_delta), delta = t_delta, valid = true),
    )
    sort!(candidates; by = c -> c.score, rev = true)
    return candidates
end

function _local_parameter_real_bounds(value)
    a = real(value)
    mid = Float64(Nemo.midpoint(a))
    rad = Float64(Nemo.radius(a))
    return mid - rad, mid + rad
end

function _local_parameter_interval_hull(sys::HCSystem, a, b)
    a_low, a_high = _local_parameter_real_bounds(a)
    b_low, b_high = _local_parameter_real_bounds(b)
    low = min(a_low, b_low)
    high = max(a_high, b_high)
    mid = (low + high) / 2
    rad = (high - low) / 2
    return sys.CC(sys.RR("$mid +/- $rad"), sys.RR(0))
end

function _local_parameter_interval_contains(interval, value; slack = 1e-12)
    interval_low, interval_high = _local_parameter_real_bounds(interval)
    value_low, value_high = _local_parameter_real_bounds(value)
    scale = max(1.0, abs(interval_low), abs(interval_high), abs(value_low), abs(value_high))
    tol = slack * scale
    return interval_low - tol <= value_low && value_high <= interval_high + tol
end

function _local_parameter_endpoint_in_box(U_box, local_parameter_interval, u_endpoint, local_parameter_endpoint)
    unknowns_ok = all(_local_parameter_interval_contains(U_box[i], u_endpoint[i]) for i in eachindex(u_endpoint))
    local_parameter_ok = _local_parameter_interval_contains(local_parameter_interval, local_parameter_endpoint)
    return (unknowns = unknowns_ok, local_parameter = local_parameter_ok, all = unknowns_ok && local_parameter_ok)
end

function _local_parameter_orientation(c0, c1; tol = 1e-14)
    delta = _local_parameter_mid_float(real(c1 - c0))
    orientation = abs(delta) <= tol ? :flat : (delta > 0 ? :increasing : :decreasing)
    return (orientation = orientation, delta = delta, valid = orientation != :flat)
end

function _local_parameter_radius_to_cover(center, value)
    center_mid = _local_parameter_mid_float(real(center))
    value_low, value_high = _local_parameter_real_bounds(value)
    return max(abs(value_low - center_mid), abs(value_high - center_mid))
end

function _local_parameter_endpoint_box_radii(sys::HCSystem, u_mid, u0, u1, inflation)
    return [
        max(
            _local_parameter_radius_to_cover(u_mid[i], u0[i]),
            _local_parameter_radius_to_cover(u_mid[i], u1[i]),
        ) + Float64(inflation)
        for i in eachindex(u_mid)
    ]
end

function _local_parameter_endpoint_box_inflation_candidates(base_radius, growth, max_radius, u_mid, u0, u1, Y, rho)
    candidates = _segment_radius_candidates(base_radius, growth, max_radius)

    half_widths = [
        max(
            _local_parameter_radius_to_cover(u_mid[i], u0[i]),
            _local_parameter_radius_to_cover(u_mid[i], u1[i]),
        )
        for i in eachindex(u_mid)
    ]
    min_half_width = minimum(half_widths)
    rho_float = Float64(rho)
    residual_bound = Float64(Y)

    if isfinite(residual_bound) && residual_bound > 0 && rho_float > 0
        target_min_radius = 4 * residual_bound / rho_float
        residual_inflation = max(0.0, target_min_radius - min_half_width)
        if residual_inflation > 0
            for factor in (0.25, 0.5, 1.0, 2.0, 5.0, 10.0)
                r = residual_inflation * factor
                r <= max_radius && push!(candidates, r)
            end
        end
    end

    sort!(candidates)
    return unique(candidates)
end

function _local_parameter_store_mode(mode, name)
    mode in (:full, :summary) ||
        throw(ArgumentError("$(name) must be one of :full or :summary."))
    return mode
end

function _local_parameter_row_value(row, key, default = missing)
    return haskey(row, key) ? getproperty(row, key) : default
end

function _local_parameter_stored_row(row; parent_index, depth, method = missing, mode = :full)
    mode = _local_parameter_store_mode(mode, :store_boxes)
    full = merge(row, (parent_index = parent_index, depth = depth))
    method !== missing && (full = merge(full, (method = method,)))
    mode === :full && return full

    return (
        parent_index = parent_index,
        depth = depth,
        success = _local_parameter_row_value(row, :success, false),
        status = _local_parameter_row_value(row, :status, missing),
        local_parameter = _local_parameter_row_value(row, :local_parameter, missing),
        krawczyk_norm = _local_parameter_row_value(row, :krawczyk_norm, Inf),
        radius = _local_parameter_row_value(row, :radius, missing),
        max_unknown_radius = _local_parameter_row_value(row, :max_unknown_radius, missing),
        min_unknown_radius = _local_parameter_row_value(row, :min_unknown_radius, missing),
        Y = _local_parameter_row_value(row, :Y, missing),
        Z = _local_parameter_row_value(row, :Z, missing),
        Y_over_r = _local_parameter_row_value(row, :Y_over_r, missing),
        yz_bound = _local_parameter_row_value(row, :yz_bound, missing),
        t_start = _local_parameter_row_value(row, :t_start, missing),
        t_end = _local_parameter_row_value(row, :t_end, missing),
        local_parameter_orientation = _local_parameter_row_value(row, :local_parameter_orientation, missing),
        local_parameter_delta = _local_parameter_row_value(row, :local_parameter_delta, missing),
        method = method === missing ? _local_parameter_row_value(row, :method, missing) : method,
    )
end

_trace_point_radius(point) = 0.0
_trace_point_radius(point::CertifiedTraceNode) = point.radius
_trace_point_radius(point::NamedTuple) = haskey(point, :radius) ? Float64(point.radius) : 0.0

function _endpoint_enclosure_x(sys::HCSystem, point)
    radius = _trace_point_radius(point)
    radius <= 0 && return point.x
    unit = sys.CC(sys.RR("0 +/- 1"), sys.RR("0 +/- 1"))
    return [xi + unit * sys.CC(radius) for xi in point.x]
end

function _validate_local_parameter_endpoint_box_segment_once(
    sys::HCSystem,
    p0,
    p1,
    local_parameter;
    rho,
    tube_radius_floor,
    radius_growth,
    max_radius,
    diagnostics,
)
    n = length(p0.x)
    t0 = sys.CC(p0.t)
    t1 = sys.CC(p1.t)
    x0_enclosure = _endpoint_enclosure_x(sys, p0)
    x1_enclosure = _endpoint_enclosure_x(sys, p1)
    c0 = _local_parameter_value(sys, p0.x, t0, local_parameter.index)
    c1 = _local_parameter_value(sys, p1.x, t1, local_parameter.index)
    u0 = _local_parameter_unknown_endpoint(sys, p0.x, t0, local_parameter.index)
    u1 = _local_parameter_unknown_endpoint(sys, p1.x, t1, local_parameter.index)
    c0_enclosure = _local_parameter_value(sys, x0_enclosure, t0, local_parameter.index)
    c1_enclosure = _local_parameter_value(sys, x1_enclosure, t1, local_parameter.index)
    u0_enclosure = _local_parameter_unknown_endpoint(sys, x0_enclosure, t0, local_parameter.index)
    u1_enclosure = _local_parameter_unknown_endpoint(sys, x1_enclosure, t1, local_parameter.index)
    u_mid = [(u0[i] + u1[i]) / sys.CC(2) for i in eachindex(u0)]
    c_mid = (c0 + c1) / sys.CC(2)
    local_parameter_interval = _local_parameter_interval_hull(sys, c0_enclosure, c1_enclosure)
    orientation = _local_parameter_orientation(c0, c1)

    A = inv_acb(_local_parameter_jacobian(sys, u_mid, local_parameter.index, c_mid, n))
    x_mid, t_mid = _local_parameter_reconstruct_x_t(sys, u_mid, local_parameter.index, local_parameter_interval, n)
    F_val = _local_parameter_real_rows(sys, evaluate_H(sys, x_mid, t_mid))

    I = Matrix{AcbFieldElem}(undef, 2n, 2n)
    for i in 1:2n, j in 1:2n
        I[i, j] = i == j ? sys.CC(1) : sys.CC(0)
    end
    unit = sys.CC(sys.RR("0 +/- 1"), sys.RR(0))
    B = fill(unit, 2n)

    base_radius = max(Float64(tube_radius_floor), 0.0)
    AH = A * F_val
    Y = norm_inf(AH)
    best = nothing
    for inflation in _local_parameter_endpoint_box_inflation_candidates(base_radius, radius_growth, max_radius, u_mid, u0_enclosure, u1_enclosure, Y, rho)
        _diag_basic(diagnostics) && (diagnostics.krawczyk_validation_calls += 1)
        u_radii = _local_parameter_endpoint_box_radii(sys, u_mid, u0_enclosure, u1_enclosure, inflation)
        U_expanded = [u_mid[i] + B[i] * sys.CC(u_radii[i]) for i in eachindex(u_mid)]
        J = _local_parameter_jacobian(sys, U_expanded, local_parameter.index, local_parameter_interval, n)
        linear_defect = I - A * J
        K = Vector{AcbFieldElem}(undef, 2n)
        for i in 1:2n
            row_term = sys.CC(0)
            for j in 1:2n
                row_term += linear_defect[i, j] * B[j] * sys.CC(u_radii[j])
            end
            K[i] = (-AH[i] + row_term) / sys.CC(u_radii[i])
        end
        start_enclosure = _local_parameter_endpoint_in_box(U_expanded, local_parameter_interval, u0_enclosure, c0_enclosure)
        end_enclosure = _local_parameter_endpoint_in_box(U_expanded, local_parameter_interval, u1_enclosure, c1_enclosure)
        endpoint_enclosure = (
            start = start_enclosure,
            finish = end_enclosure,
            both = start_enclosure.all && end_enclosure.all,
        )
        k_norm = norm_inf(K)
        Z = _matrix_inf_norm_bound(linear_defect)
        min_radius = minimum(u_radii)
        _record_yz_max!(diagnostics, Y, Z, Y / min_radius)
        row = (
            success = k_norm < rho && endpoint_enclosure.both && orientation.valid,
            status = !(k_norm < rho) ? :krawczyk_failed :
                !endpoint_enclosure.both ? :endpoint_enclosure_failed :
                !orientation.valid ? :local_parameter_orientation_failed :
                :success,
            local_parameter = local_parameter,
            krawczyk_norm = k_norm,
            radius = Float64(inflation),
            max_unknown_radius = maximum(u_radii),
            min_unknown_radius = min_radius,
            Y = Y,
            Z = Z,
            Y_over_r = Y / min_radius,
            yz_bound = Y / min_radius + Z,
            unknown_box = U_expanded,
            local_parameter_interval = local_parameter_interval,
            t_start = p0.t,
            t_end = p1.t,
            x_start = p0.x,
            x_end = p1.x,
            endpoint_enclosure = endpoint_enclosure,
            local_parameter_orientation = orientation.orientation,
            local_parameter_delta = orientation.delta,
            orientation_valid = orientation.valid,
            method = :local_parameter_endpoint_box,
        )
        row.success && begin
            _record_local_parameter_choice!(diagnostics, local_parameter)
            return row
        end
        if best === nothing || row.krawczyk_norm < best.krawczyk_norm
            best = row
        end
    end
    return best
end

function _local_parameter_midpoint_point(sys::HCSystem, p0, p1; midpoint_policy, node_refinement_radius, tau)
    t_mid = (p0.t + p1.t) / 2
    x_guess = (p0.x .+ p1.x) ./ sys.CC(2)
    if midpoint_policy in (:krawczyk_polish, :krawczyk, :moore)
        midpoint = _certify_trace_node(
            sys,
            (t = t_mid, x = x_guess);
            node_refinement_radius = node_refinement_radius,
            tau = tau,
            diagnostics = nothing,
            refinement_counter = :midpoint_refinements,
        )
        midpoint !== nothing && return midpoint
        return nothing
    elseif midpoint_policy != :newton && midpoint_policy != :average
        throw(ArgumentError("unsupported midpoint_policy $(midpoint_policy)"))
    end
    return (t = t_mid, x = x_guess)
end

function _certify_local_parameter_endpoint_box_segment!(
    boxes,
    failures,
    sys::HCSystem,
    p0,
    p1;
    parent_index,
    depth,
    max_depth,
    rho,
    tau,
    node_refinement_radius,
    tube_radius_floor,
    radius_growth,
    max_radius,
    local_parameter_variables,
    midpoint_policy,
    store_boxes,
    store_failures,
    diagnostics,
    progress_state = nothing,
    fail_fast = true,
)
    _diag_basic(diagnostics) && (diagnostics.segment_attempts += 1)
    _local_parameter_maybe_show_progress!(progress_state, boxes, failures, depth)
    candidates = nothing
    if _diag_timing(diagnostics)
        elapsed = @elapsed begin
            candidates = _local_parameter_candidates(sys, p0, p1)
            variable_indices = _selected_local_parameter_vars(sys, local_parameter_variables)
            candidates = _filter_local_parameter_candidates(candidates, variable_indices)
        end
        diagnostics.time_in_local_parameter_search += elapsed
    else
        candidates = _local_parameter_candidates(sys, p0, p1)
        variable_indices = _selected_local_parameter_vars(sys, local_parameter_variables)
        candidates = _filter_local_parameter_candidates(candidates, variable_indices)
    end
    if isempty(candidates)
        push!(failures, (parent_index = parent_index, depth = depth, status = :no_local_parameter_candidates))
        return false
    end

    best = nothing
    for local_parameter in candidates
        row = try
            if _diag_timing(diagnostics)
                timed_row = nothing
                elapsed = @elapsed begin
                    timed_row = _validate_local_parameter_endpoint_box_segment_once(
                        sys,
                        p0,
                        p1,
                        local_parameter;
                        rho = rho,
                        tube_radius_floor = tube_radius_floor,
                        radius_growth = radius_growth,
                        max_radius = max_radius,
                        diagnostics = diagnostics,
                    )
                end
                diagnostics.time_in_validation += elapsed
                timed_row
            else
                _validate_local_parameter_endpoint_box_segment_once(
                    sys,
                    p0,
                    p1,
                    local_parameter;
                    rho = rho,
                    tube_radius_floor = tube_radius_floor,
                    radius_growth = radius_growth,
                    max_radius = max_radius,
                    diagnostics = diagnostics,
                )
            end
        catch err
            err isa InterruptException && rethrow(err)
            nothing
        end
        row === nothing && continue
        row.success && begin
            push!(
                boxes,
                _local_parameter_stored_row(
                    row;
                    parent_index = parent_index,
                    depth = depth,
                    mode = store_boxes,
                ),
            )
            return true
        end
        if best === nothing || row.krawczyk_norm < best.krawczyk_norm
            best = row
        end
    end

    if depth >= max_depth
        push!(
            failures,
            best === nothing ?
                (parent_index = parent_index, depth = depth, status = :all_local_parameter_attempts_failed) :
                _local_parameter_stored_row(
                    best;
                    parent_index = parent_index,
                    depth = depth,
                    mode = store_failures,
                ),
        )
        return false
    end

    midpoint = _local_parameter_midpoint_point(
        sys,
        p0,
        p1;
        midpoint_policy = midpoint_policy,
        node_refinement_radius = node_refinement_radius,
        tau = tau,
    )
    if midpoint === nothing
        push!(failures, (parent_index = parent_index, depth = depth, status = :split_refinement_failed))
        return false
    end
    left = _certify_local_parameter_endpoint_box_segment!(
        boxes,
        failures,
        sys,
        p0,
        midpoint;
        parent_index = parent_index,
        depth = depth + 1,
        max_depth = max_depth,
        rho = rho,
        tau = tau,
        node_refinement_radius = node_refinement_radius,
        tube_radius_floor = tube_radius_floor,
        radius_growth = radius_growth,
        max_radius = max_radius,
        local_parameter_variables = local_parameter_variables,
        midpoint_policy = midpoint_policy,
        store_boxes = store_boxes,
        store_failures = store_failures,
        diagnostics = diagnostics,
        progress_state = progress_state,
        fail_fast = fail_fast,
    )
    !left && fail_fast && return false
    right = _certify_local_parameter_endpoint_box_segment!(
        boxes,
        failures,
        sys,
        midpoint,
        p1;
        parent_index = parent_index,
        depth = depth + 1,
        max_depth = max_depth,
        rho = rho,
        tau = tau,
        node_refinement_radius = node_refinement_radius,
        tube_radius_floor = tube_radius_floor,
        radius_growth = radius_growth,
        max_radius = max_radius,
        local_parameter_variables = local_parameter_variables,
        midpoint_policy = midpoint_policy,
        store_boxes = store_boxes,
        store_failures = store_failures,
        diagnostics = diagnostics,
        progress_state = progress_state,
        fail_fast = fail_fast,
    )
    return left && right
end

_local_parameter_variable_symbol(v) = Symbol(v)
_local_parameter_variable_symbol(v::Symbol) = v
_local_parameter_variable_symbol(v::AbstractString) = Symbol(v)

function _selected_local_parameter_vars(sys::HCSystem, local_parameter_variables)
    n = sys.compiled.n_vars == 0 ? length(sys.p_start) : sys.compiled.n_vars
    isempty(local_parameter_variables) && return Int[]

    source = sys.source === nothing ? sys.compiled.source : sys.source
    name_to_index = Dict{Symbol,Int}()
    if source !== nothing
        source_vars = collect(source.variables)
        for (i, var) in enumerate(source_vars)
            name_to_index[_local_parameter_variable_symbol(var)] = i
        end
    end

    indices = Int[]
    for var in local_parameter_variables
        if var isa Integer
            1 <= Int(var) <= n || throw(ArgumentError("local_parameter variable index $(var) is outside 1:$n."))
            push!(indices, Int(var))
        else
            source === nothing &&
                throw(ArgumentError("named local_parameter_variables require source metadata; use integer indices instead."))
            sym = _local_parameter_variable_symbol(var)
            haskey(name_to_index, sym) ||
                throw(ArgumentError("local_parameter variable $(sym) is not one of the source variables $(keys(name_to_index))."))
            push!(indices, name_to_index[sym])
        end
    end
    return unique(indices)
end

_is_t_local_parameter_candidate(candidate) = candidate.name === :t

function _filter_local_parameter_candidates(candidates, variable_indices)
    isempty(variable_indices) && return candidates
    allowed = Set{Int}()
    for i in variable_indices
        push!(allowed, 2i - 1)
        push!(allowed, 2i)
    end
    return filter(c -> c.index in allowed || _is_t_local_parameter_candidate(c), candidates)
end

function _local_parameter_interval(sys::HCSystem, c0, c1)
    a0 = _local_parameter_mid_float(c0)
    a1 = _local_parameter_mid_float(c1)
    mid = (a0 + a1) / 2
    rad = abs(a1 - a0) / 2
    return sys.CC(sys.RR("$mid +/- $rad"), sys.RR(0))
end

function _local_parameter_jacobian(sys::HCSystem, u, local_parameter_index, local_parameter_value, n::Integer)
    x, t = _local_parameter_reconstruct_x_t(sys, u, local_parameter_index, local_parameter_value, n)
    Jx = evaluate_Jac(sys, x, t)
    dt = evaluate_dt(sys, x, t)
    Jr = Matrix{AcbFieldElem}(undef, 2n, 2n)
    col = 1
    for j in 1:n
        if 2j - 1 != local_parameter_index
            for i in 1:n
                _local_parameter_fill_real_column!(Jr, col, 2i - 1, Jx[i, j], false, sys)
            end
            col += 1
        end
        if 2j != local_parameter_index
            for i in 1:n
                _local_parameter_fill_real_column!(Jr, col, 2i - 1, Jx[i, j], true, sys)
            end
            col += 1
        end
    end
    if local_parameter_index != 2n + 1
        for i in 1:n
            Jr[2i - 1, col] = _local_parameter_real_part(sys, dt[i])
            Jr[2i, col] = _local_parameter_imag_part(sys, dt[i])
        end
    end
    return Jr
end

function _local_parameter_column(sys::HCSystem, x, t, local_parameter_index)
    n = length(x)
    col = Vector{AcbFieldElem}(undef, 2n)
    if local_parameter_index == 2n + 1
        dt = evaluate_dt(sys, x, t)
        for i in 1:n
            col[2i - 1] = _local_parameter_real_part(sys, dt[i])
            col[2i] = _local_parameter_imag_part(sys, dt[i])
        end
        return col
    end

    var_idx = (local_parameter_index + 1) ÷ 2
    imag_coordinate = iseven(local_parameter_index)
    Jx = evaluate_Jac(sys, x, t)
    for i in 1:n
        if imag_coordinate
            col[2i - 1] = -_local_parameter_imag_part(sys, Jx[i, var_idx])
            col[2i] = _local_parameter_real_part(sys, Jx[i, var_idx])
        else
            col[2i - 1] = _local_parameter_real_part(sys, Jx[i, var_idx])
            col[2i] = _local_parameter_imag_part(sys, Jx[i, var_idx])
        end
    end
    return col
end

function _local_parameter_unknown_velocity(sys::HCSystem, node::CertifiedTraceNode, local_parameter_index)
    n = length(node.x)
    u = _local_parameter_unknown_endpoint(sys, node.x, sys.CC(node.t), local_parameter_index)
    c = _local_parameter_value(sys, node.x, sys.CC(node.t), local_parameter_index)
    Ju = _local_parameter_jacobian(sys, u, local_parameter_index, c, n)
    Js = _local_parameter_column(sys, node.x, sys.CC(node.t), local_parameter_index)
    return -(inv_acb(Ju) * Js)
end

function _construct_local_parameter_hermite_tm(sys::HCSystem, u0, u1, du0, du1, c0, c1)
    CC = sys.CC
    RR = sys.RR
    h_cc = c1 - c0
    h_mid = _local_parameter_mid_float(real(h_cc))
    abs(h_mid) > 0 || throw(ArgumentError("local_parameter interval has zero width."))
    orientation = h_mid > 0 ? 1 : -1
    h2 = h_cc^2
    h3 = h_cc^3
    h_arb = RR(abs(h_mid))
    direction = CC(orientation)

    tm_vec = Vector{TaylorModel3{AcbFieldElem,ArbFieldElem}}(undef, length(u0))
    for i in eachindex(u0)
        dx = u1[i] - u0[i]
        a0 = u0[i]
        a1 = du0[i]
        a2 = 3 * dx / h2 - (2 * du0[i] + du1[i]) / h_cc
        a3 = (du0[i] + du1[i]) / h2 - 2 * dx / h3
        tm_vec[i] = TaylorModel3(a0, a1 * direction, a2, a3 * direction, CC(0), h_arb)
    end
    return tm_vec
end

function _validate_local_parameter_segment_once(
    sys::HCSystem,
    node0::CertifiedTraceNode,
    node1::CertifiedTraceNode,
    local_parameter;
    rho,
    tube_radius_floor,
    radius_growth,
    max_radius,
    diagnostics,
)
    n = length(node0.x)
    c0 = _local_parameter_value(sys, node0.x, sys.CC(node0.t), local_parameter.index)
    c1 = _local_parameter_value(sys, node1.x, sys.CC(node1.t), local_parameter.index)
    u0 = _local_parameter_unknown_endpoint(sys, node0.x, sys.CC(node0.t), local_parameter.index)
    u1 = _local_parameter_unknown_endpoint(sys, node1.x, sys.CC(node1.t), local_parameter.index)
    du0 = _local_parameter_unknown_velocity(sys, node0, local_parameter.index)
    du1 = _local_parameter_unknown_velocity(sys, node1, local_parameter.index)
    U_tm = _construct_local_parameter_hermite_tm(sys, u0, u1, du0, du1, c0, c1)

    c_mid = (c0 + c1) / sys.CC(2)
    local_parameter_delta = _local_parameter_mid_float(real(c1 - c0))
    local_parameter_abs = abs(local_parameter_delta)
    local_parameter_orientation = local_parameter_delta > 0 ? 1 : -1
    local_parameter_tm = TaylorModel3(c0, sys.CC(local_parameter_orientation), sys.CC(0), sys.CC(0), sys.CC(0), sys.RR(local_parameter_abs))
    x_tm, t_tm = _local_parameter_reconstruct_x_t(sys, U_tm, local_parameter.index, local_parameter_tm, n)

    U_mid = [(u0[i] + u1[i]) / sys.CC(2) for i in eachindex(u0)]
    A = inv_acb(_local_parameter_jacobian(sys, U_mid, local_parameter.index, c_mid, n))
    F_tm = _local_parameter_real_rows_tm(sys, evaluate_H(sys, x_tm, t_tm))

    F_val = Vector{AcbFieldElem}(undef, 2n)
    U_bound = Vector{AcbFieldElem}(undef, 2n)
    for i in 1:2n
        F_val[i] = evaluate_taylor(F_tm[i])
        U_bound[i] = evaluate_taylor(U_tm[i])
    end

    base_radius = max(Float64(tube_radius_floor), 2 * max(node0.radius, node1.radius))
    unit = sys.CC(sys.RR("0 +/- 1"), sys.RR("0 +/- 1"))
    B = fill(unit, 2n)
    I = Matrix{AcbFieldElem}(undef, 2n, 2n)
    for i in 1:2n, j in 1:2n
        I[i, j] = i == j ? sys.CC(1) : sys.CC(0)
    end

    best = nothing
    for radius in _segment_radius_candidates(base_radius, radius_growth, max_radius)
        _diag_basic(diagnostics) && (diagnostics.krawczyk_validation_calls += 1)
        U_expanded = U_bound .+ (sys.CC(radius) .* B)
        local_parameter_interval = _local_parameter_interval(sys, c0, c1)
        J = _local_parameter_jacobian(sys, U_expanded, local_parameter.index, local_parameter_interval, n)
        AH = A * F_val
        linear_defect = I - A * J
        K = -AH ./ sys.CC(radius) + linear_defect * B
        Y = norm_inf(AH)
        Z = _matrix_inf_norm_bound(linear_defect)
        Y_over_r = Y / Float64(radius)
        _record_yz_max!(diagnostics, Y, Z, Y_over_r)
        row = (
            success = norm_inf(K) < rho,
            status = norm_inf(K) < rho ? :success : :krawczyk_failed,
            local_parameter = local_parameter,
            krawczyk_norm = norm_inf(K),
            radius = Float64(radius),
            Y = Y,
            Z = Z,
            Y_over_r = Y_over_r,
            yz_bound = Y_over_r + Z,
            unknown_box = U_expanded,
            local_parameter_interval = local_parameter_interval,
            t_start = node0.t,
            t_end = node1.t,
            x_start = node0.x,
            x_end = node1.x,
        )
        row.success && begin
            _record_local_parameter_choice!(diagnostics, local_parameter)
            return row
        end
        if best === nothing || row.krawczyk_norm < best.krawczyk_norm
            best = row
        end
    end
    return best
end

function _certify_local_parameter_segment!(
    boxes,
    failures,
    sys::HCSystem,
    node0::CertifiedTraceNode,
    node1::CertifiedTraceNode;
    parent_index,
    depth,
    max_depth,
    rho,
    tau,
    node_refinement_radius,
    tube_radius_floor,
    radius_growth,
    max_radius,
    diagnostics,
    local_parameter_variables,
    fail_fast = true,
)
    _diag_basic(diagnostics) && (diagnostics.segment_attempts += 1)
    candidates = nothing
    if _diag_timing(diagnostics)
        elapsed = @elapsed begin
            candidates = _local_parameter_candidates(sys, node0, node1)
            variable_indices = _selected_local_parameter_vars(sys, local_parameter_variables)
            candidates = _filter_local_parameter_candidates(candidates, variable_indices)
        end
        diagnostics.time_in_local_parameter_search += elapsed
    else
        candidates = _local_parameter_candidates(sys, node0, node1)
        variable_indices = _selected_local_parameter_vars(sys, local_parameter_variables)
        candidates = _filter_local_parameter_candidates(candidates, variable_indices)
    end
    isempty(candidates) && begin
        push!(failures, (parent_index = parent_index, depth = depth, status = :no_local_parameter_candidates))
        return false
    end

    best = nothing
    for local_parameter in candidates
        row = try
            if _diag_timing(diagnostics)
                timed_row = nothing
                elapsed = @elapsed begin
                    timed_row = _validate_local_parameter_segment_once(
                        sys,
                        node0,
                        node1,
                        local_parameter;
                        rho = rho,
                        tube_radius_floor = tube_radius_floor,
                        radius_growth = radius_growth,
                        max_radius = max_radius,
                        diagnostics = diagnostics,
                    )
                end
                diagnostics.time_in_validation += elapsed
                timed_row
            else
                _validate_local_parameter_segment_once(
                    sys,
                    node0,
                    node1,
                    local_parameter;
                    rho = rho,
                    tube_radius_floor = tube_radius_floor,
                    radius_growth = radius_growth,
                    max_radius = max_radius,
                    diagnostics = diagnostics,
                )
            end
        catch err
            err isa InterruptException && rethrow(err)
            nothing
        end
        row === nothing && continue
        row.success && begin
            push!(boxes, merge(row, (parent_index = parent_index, depth = depth, method = :adaptive_local_parameter_hermite)))
            return true
        end
        if best === nothing || row.krawczyk_norm < best.krawczyk_norm
            best = row
        end
    end

    if depth >= max_depth
        push!(
            failures,
            best === nothing ?
                (parent_index = parent_index, depth = depth, status = :all_local_parameter_attempts_failed) :
                merge(best, (parent_index = parent_index, depth = depth, method = :adaptive_local_parameter_hermite)),
        )
        return false
    end

    t_mid = (node0.t + node1.t) / 2
    x_guess = (node0.x .+ node1.x) ./ sys.CC(2)
    midpoint = _certify_trace_node(
        sys,
        (t = t_mid, x = x_guess);
        node_refinement_radius = node_refinement_radius,
        tau = tau,
        diagnostics = diagnostics,
        refinement_counter = :midpoint_refinements,
    )
    if midpoint === nothing
        push!(failures, (parent_index = parent_index, depth = depth, status = :split_refinement_failed))
        return false
    end

    left = _certify_local_parameter_segment!(
        boxes,
        failures,
        sys,
        node0,
        midpoint;
        parent_index = parent_index,
        depth = depth + 1,
        max_depth = max_depth,
        rho = rho,
        tau = tau,
        node_refinement_radius = node_refinement_radius,
        tube_radius_floor = tube_radius_floor,
        radius_growth = radius_growth,
        max_radius = max_radius,
        diagnostics = diagnostics,
        local_parameter_variables = local_parameter_variables,
        fail_fast = fail_fast,
    )
    !left && fail_fast && return false
    right = _certify_local_parameter_segment!(
        boxes,
        failures,
        sys,
        midpoint,
        node1;
        parent_index = parent_index,
        depth = depth + 1,
        max_depth = max_depth,
        rho = rho,
        tau = tau,
        node_refinement_radius = node_refinement_radius,
        tube_radius_floor = tube_radius_floor,
        radius_growth = radius_growth,
        max_radius = max_radius,
        diagnostics = diagnostics,
        local_parameter_variables = local_parameter_variables,
        fail_fast = fail_fast,
    )
    return left && right
end

"""
    certify_hc_trace_adaptive_local_parameter(sys, trace; kwargs...)

Experimental a posteriori trace verifier which may use any real coordinate
`Re(x_i)`, `Im(x_i)`, or `t` as the local local_parameter on each segment.

The default endpoint-box method validates each segment with an interval
Krawczyk test in an adaptively chosen local parameter. The optional Hermite
method uses cached endpoint refinement data and evaluates a Hermite predictor
residual with `TaylorModel3`.
"""
function certify_hc_trace_adaptive_local_parameter(
    sys::HCSystem,
    trace;
    time_map = identity,
    rho = 7//8,
    tau = 0.125,
    node_refinement_radius = 1e-8,
    tube_radius_floor = 1e-12,
    radius_growth = (1.0, 10.0),
    max_radius = 1e-1,
    max_depth = 8,
    show_progress = false,
    diagnostics = :off,
    local_parameter_variables = (),
    local_parameter_method = :endpoint_box,
    midpoint_policy = :krawczyk_polish,
    store_boxes = :full,
    store_failures = :summary,
    fail_fast = true,
)
    points = _prepare_trace_points(sys, trace; time_map = time_map)
    diagnostics_data = _trace_cert_diagnostics(diagnostics)
    active_diagnostics = diagnostics_data.mode === :off ? nothing : diagnostics_data
    store_boxes = _local_parameter_store_mode(store_boxes, :store_boxes)
    store_failures = _local_parameter_store_mode(store_failures, :store_failures)
    if local_parameter_method in (:endpoint_box, :box)
        nodes = try
            _prepare_certified_trace_nodes(
                sys,
                points;
                node_refinement_radius = node_refinement_radius,
                tau = tau,
                diagnostics = active_diagnostics,
            )
        catch err
            err isa InterruptException && rethrow(err)
            return _node_refinement_failure_result(
                points,
                err,
                :adaptive_local_parameter_endpoint_box,
                diagnostics_data,
            )
        end
        boxes = NamedTuple[]
        failures = NamedTuple[]
        total = length(nodes) - 1
        for i in 1:total
            progress_state = show_progress ? _local_parameter_progress_state(i, total) : nothing
            ok = _certify_local_parameter_endpoint_box_segment!(
                boxes,
                failures,
                sys,
                nodes[i],
                nodes[i+1];
                parent_index = i,
                depth = 0,
                max_depth = max_depth,
                rho = rho,
                tau = tau,
                node_refinement_radius = node_refinement_radius,
                tube_radius_floor = tube_radius_floor,
                radius_growth = radius_growth,
                max_radius = max_radius,
                local_parameter_variables = local_parameter_variables,
                midpoint_policy = midpoint_policy,
                store_boxes = store_boxes,
                store_failures = store_failures,
                diagnostics = active_diagnostics,
                progress_state = progress_state,
                fail_fast = fail_fast,
            )
            show_progress && _local_parameter_show_progress(
                i,
                total,
                boxes,
                failures,
                progress_state === nothing ? (isempty(boxes) ? 0 : maximum(b -> b.depth, boxes)) : progress_state.max_depth[];
                label = "adaptive local-parameter endpoint-box",
                attempts = progress_state === nothing ? nothing : progress_state.attempts[],
                final = i == total,
            )
            !ok && fail_fast && break
        end
        return (
            success = isempty(failures),
            boxes = boxes,
            segments = boxes,
            failed_segments = failures,
            total_boxes = length(boxes),
            original_segments = length(points) - 1,
            max_depth = isempty(boxes) ? 0 : maximum(b -> b.depth, boxes),
            max_krawczyk_norm = isempty(boxes) ? Inf : maximum(b -> b.krawczyk_norm, boxes),
            method = :adaptive_local_parameter_endpoint_box,
            boxes_by_parent = _count_by_field(boxes, :parent_index),
            boxes_by_depth = _count_by_field(boxes, :depth),
            diagnostics = _freeze_trace_cert_diagnostics(diagnostics_data),
        )
    elseif local_parameter_method != :hermite
        throw(ArgumentError("unsupported local_parameter_method $(local_parameter_method)"))
    end

    nodes = try
        _prepare_certified_trace_nodes(
            sys,
            points;
            node_refinement_radius = node_refinement_radius,
            tau = tau,
            diagnostics = active_diagnostics,
        )
    catch err
        err isa InterruptException && rethrow(err)
        return _node_refinement_failure_result(
            points,
            err,
            :adaptive_local_parameter_hermite,
            diagnostics_data,
        )
    end

    boxes = NamedTuple[]
    failures = NamedTuple[]
    total = length(nodes) - 1
    for i in 1:total
        ok = _certify_local_parameter_segment!(
            boxes,
            failures,
            sys,
            nodes[i],
            nodes[i+1];
            parent_index = i,
            depth = 0,
            max_depth = max_depth,
            rho = rho,
            tau = tau,
            node_refinement_radius = node_refinement_radius,
            tube_radius_floor = tube_radius_floor,
            radius_growth = radius_growth,
            max_radius = max_radius,
            diagnostics = active_diagnostics,
            local_parameter_variables = local_parameter_variables,
            fail_fast = fail_fast,
        )
        show_progress && _local_parameter_show_progress(
            i,
            total,
            boxes,
            failures,
            isempty(boxes) ? 0 : maximum(b -> b.depth, boxes);
            label = "adaptive local-parameter",
            final = i == total,
        )
        !ok && fail_fast && break
    end

    return (
        success = isempty(failures),
        boxes = boxes,
        segments = boxes,
        failed_segments = failures,
        total_boxes = length(boxes),
        original_segments = length(nodes) - 1,
        max_depth = isempty(boxes) ? 0 : maximum(b -> b.depth, boxes),
        max_krawczyk_norm = isempty(boxes) ? Inf : maximum(b -> b.krawczyk_norm, boxes),
        method = :adaptive_local_parameter_hermite,
        boxes_by_parent = _count_by_field(boxes, :parent_index),
        boxes_by_depth = _count_by_field(boxes, :depth),
        diagnostics = _freeze_trace_cert_diagnostics(diagnostics_data),
    )
end
