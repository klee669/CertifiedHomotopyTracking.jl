export PosterioriTracker, PosterioriPathResult, prepare_posteriori_tracker, posteriori_hc_system,
       posteriori_hc_homotopy, certify_path_a_posteriori

import HomotopyContinuation

_complex_f64_vector(values) = ComplexF64[ComplexF64(value) for value in values]

struct PosterioriPathResult
    data::NamedTuple
end

function Base.getproperty(result::PosterioriPathResult, name::Symbol)
    name === :data && return getfield(result, :data)
    return getproperty(getfield(result, :data), name)
end

Base.propertynames(result::PosterioriPathResult; private::Bool = false) =
    private ? (:data, propertynames(getfield(result, :data))...) : propertynames(getfield(result, :data))
Base.haskey(result::PosterioriPathResult, key::Symbol) = haskey(getfield(result, :data), key)
Base.keys(result::PosterioriPathResult) = keys(getfield(result, :data))
Base.getindex(result::PosterioriPathResult, key::Symbol) = getproperty(result, key)

function Base.show(io::IO, result::PosterioriPathResult)
    data = getfield(result, :data)
    success = get(data, :success, missing)
    method = get(data, :method, missing)
    total_boxes = get(data, :total_boxes, missing)
    original_segments = get(data, :original_segments, missing)
    max_depth = get(data, :max_depth, missing)
    failed_count = haskey(data, :failed_segments) ? length(data.failed_segments) : missing
    trace_points = haskey(data, :hc_trace) && haskey(data.hc_trace, :trace) ? length(data.hc_trace.trace) : missing
    hc_status = haskey(data, :hc_trace) && haskey(data.hc_trace, :status) ? data.hc_trace.status : missing

    print(io, "PosterioriPathResult(")
    print(io, "success=", success)
    print(io, ", method=", method)
    print(io, ", trace_points=", trace_points)
    print(io, ", original_segments=", original_segments)
    print(io, ", boxes=", total_boxes)
    print(io, ", failed=", failed_count)
    print(io, ", max_depth=", max_depth)
    print(io, ", hc_status=", hc_status)
    print(io, ")")
end

Base.show(io::IO, ::MIME"text/plain", result::PosterioriPathResult) = show(io, result)

mutable struct PosterioriTracker
    compiled::CompiledHomotopy
    source::HomotopySourceData
    cht_system_cache::Any
    hc_system_cache::Any
    hc_homotopy_cache::Any
    hc_variables_cache::Any
    hc_parameters_cache::Any
    hc_parameter_values_cache::Any
    hc_t_cache::Any
end

function prepare_posteriori_tracker(compiled::CompiledHomotopy)
    compiled.source === nothing &&
        throw(ArgumentError("compiled homotopy does not contain source metadata. Recompile it with compile_edge_homotopy or compile_homotopy."))
    return PosterioriTracker(compiled, compiled.source, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
end

function prepare_posteriori_tracker(sys::HCSystem)
    source = sys.source === nothing ? sys.compiled.source : sys.source
    source === nothing &&
        throw(ArgumentError("specialized homotopy does not contain source metadata. Rebuild it with source-aware constructors."))
    tracker = PosterioriTracker(sys.compiled, source, sys, nothing, nothing, nothing, nothing, nothing, nothing)
    tracker.cht_system_cache = sys
    if source.kind === :direct && (!isempty(sys.p_start) || !isempty(sys.p_const))
        tracker.hc_parameter_values_cache = ComplexF64[
            _complex_f64_vector(sys.p_start);
            _complex_f64_vector(sys.p_const)
        ]
    end
    return tracker
end

_source_variable_name(v) = String(Symbol(v))

function _make_hc_variable_map(source::HomotopySourceData)
    vars = HomotopyContinuation.Variable[HomotopyContinuation.Variable(_source_variable_name(v)) for v in source.variables]
    params = HomotopyContinuation.Variable[HomotopyContinuation.Variable(_source_variable_name(p)) for p in source.parameters]
    consts = HomotopyContinuation.Variable[HomotopyContinuation.Variable(_source_variable_name(c)) for c in source.const_vars]
    t_var = source.t_var === nothing ? nothing : HomotopyContinuation.Variable(_source_variable_name(source.t_var))

    mapping = Dict{Any,Any}()
    for (src, dst) in zip(source.variables, vars)
        mapping[src] = dst
        mapping[Symbolics.unwrap(src)] = dst
    end
    for (src, dst) in zip(source.parameters, params)
        mapping[src] = dst
        mapping[Symbolics.unwrap(src)] = dst
    end
    for (src, dst) in zip(source.const_vars, consts)
        mapping[src] = dst
        mapping[Symbolics.unwrap(src)] = dst
    end
    if source.t_var !== nothing
        mapping[source.t_var] = t_var
        mapping[Symbolics.unwrap(source.t_var)] = t_var
    end

    return mapping, vars, params, consts, t_var
end

function _symbolics_to_hc_expr(expr, mapping::Dict{Any,Any})
    if haskey(mapping, expr)
        return mapping[expr]
    elseif expr isa Complex
        return _symbolics_to_hc_expr(real(expr), mapping) +
               im * _symbolics_to_hc_expr(imag(expr), mapping)
    elseif expr isa Number && !(expr isa Num)
        return expr
    end

    unwrapped = try
        Symbolics.unwrap(expr)
    catch
        expr
    end
    haskey(mapping, unwrapped) && return mapping[unwrapped]
    unwrapped isa Number && return unwrapped
    if !Symbolics.istree(unwrapped)
        literal = try
            Symbolics.value(unwrapped)
        catch
            nothing
        end
        literal isa Number && return literal
    end

    if Symbolics.istree(unwrapped)
        op = Symbolics.operation(unwrapped)
        args = [_symbolics_to_hc_expr(arg, mapping) for arg in Symbolics.arguments(unwrapped)]
        if op === +
            return sum(args)
        elseif op === *
            return prod(args)
        elseif op === ^
            return args[1] ^ args[2]
        elseif op === /
            return args[1] / args[2]
        else
            return op(args...)
        end
    end

    throw(ArgumentError("cannot convert symbolic expression to HomotopyContinuation.jl expression: $expr"))
end

function _build_posteriori_hc_system!(tracker::PosterioriTracker)
    source = tracker.source
    source.projective &&
        throw(ArgumentError("posteriori HC.jl source reconstruction for projective compiled homotopies is not implemented yet."))

    mapping, vars, params, consts, t_var = _make_hc_variable_map(source)
    equations = [_symbolics_to_hc_expr(eq, mapping) for eq in source.equations]

    tracker.hc_variables_cache = vars
    tracker.hc_parameters_cache = HomotopyContinuation.Variable[params; consts]
    tracker.hc_t_cache = t_var

    if source.kind === :edge
        tracker.hc_system_cache = HomotopyContinuation.System(
            equations;
            variables = vars,
            parameters = tracker.hc_parameters_cache,
        )
    elseif source.kind === :direct
        tracker.hc_homotopy_cache = HomotopyContinuation.Homotopy(
            equations,
            vars,
            t_var;
            parameters = tracker.hc_parameters_cache,
        )
    else
        throw(ArgumentError("unknown homotopy source kind: $(source.kind)"))
    end

    return tracker
end

function posteriori_hc_system(tracker::PosterioriTracker)
    tracker.source.kind === :edge ||
        throw(ArgumentError("posteriori_hc_system is only defined for edge homotopy sources."))
    tracker.hc_system_cache === nothing && _build_posteriori_hc_system!(tracker)
    return tracker.hc_system_cache
end

function posteriori_hc_homotopy(
    tracker::PosterioriTracker;
    start_parameters = nothing,
    target_parameters = nothing,
    const_parameters = ComplexF64[],
    parameter_values = nothing,
)
    source = tracker.source
    if source.kind === :direct
        tracker.hc_homotopy_cache === nothing && _build_posteriori_hc_system!(tracker)
        values = parameter_values === nothing ? tracker.hc_parameter_values_cache : parameter_values
        values === nothing && return tracker.hc_homotopy_cache
        return HomotopyContinuation.fix_parameters(tracker.hc_homotopy_cache, _complex_f64_vector(values))
    elseif source.kind === :edge
        if tracker.cht_system_cache !== nothing
            start_parameters === nothing && (start_parameters = tracker.cht_system_cache.p_start)
            target_parameters === nothing && (target_parameters = tracker.cht_system_cache.p_end)
            isempty(const_parameters) && (const_parameters = tracker.cht_system_cache.p_const)
        end
        start_parameters === nothing &&
            throw(ArgumentError("start_parameters are required for an edge homotopy source."))
        target_parameters === nothing &&
            throw(ArgumentError("target_parameters are required for an edge homotopy source."))
        tracker.hc_system_cache === nothing && _build_posteriori_hc_system!(tracker)
        start = ComplexF64[_complex_f64_vector(start_parameters); _complex_f64_vector(const_parameters)]
        target = ComplexF64[_complex_f64_vector(target_parameters); _complex_f64_vector(const_parameters)]
        length(start) == length(tracker.hc_parameters_cache) ||
            throw(DimensionMismatch("start_parameters plus const_parameters must match the stored source parameters."))
        length(target) == length(tracker.hc_parameters_cache) ||
            throw(DimensionMismatch("target_parameters plus const_parameters must match the stored source parameters."))
        return HomotopyContinuation.ParameterHomotopy(
            tracker.hc_system_cache;
            start_parameters = start,
            target_parameters = target,
        )
    else
        throw(ArgumentError("unknown homotopy source kind: $(source.kind)"))
    end
end

function _posteriori_default_t_start(source::HomotopySourceData)
    return source.kind === :direct ? 0.0 : 1.0
end

function _posteriori_default_t_target(source::HomotopySourceData)
    return source.kind === :direct ? 1.0 : 0.0
end

function _posteriori_default_time_map(source::HomotopySourceData)
    return source.kind === :edge ? (t -> 1.0 - t) : identity
end

const DEFAULT_POSTERIORI_DEPTH_SCHEDULE = (6, 8, 10, 12)
const DEFAULT_POSTERIORI_NODE_REFINEMENT_RADIUS = 1e-8
const DEFAULT_POSTERIORI_REFINEMENT_TAU = 0.125
const DEFAULT_POSTERIORI_RADIUS_GROWTH = (1.0, 10.0)
const DEFAULT_POSTERIORI_MAX_RADIUS = 1e-1

function _posteriori_max_step_schedule(max_step_size)
    max_step_size === nothing && return (Inf, 0.005, 0.0025)
    if max_step_size isa Tuple || max_step_size isa AbstractVector
        isempty(max_step_size) &&
            throw(ArgumentError("max_step_size schedule must contain at least one value."))
        schedule = Float64.(max_step_size)
        all(x -> isinf(x) || x > 0, schedule) ||
            throw(ArgumentError("max_step_size schedule entries must be positive or Inf."))
        return Tuple(schedule)
    end
    max_step_size isa Real ||
        throw(ArgumentError("max_step_size must be a real number, a schedule, or nothing."))
    (isinf(Float64(max_step_size)) || Float64(max_step_size) > 0) ||
        throw(ArgumentError("max_step_size must be positive or Inf."))
    return (Float64(max_step_size),)
end

function _posteriori_depth_schedule(max_depth)
    max_depth === nothing && return DEFAULT_POSTERIORI_DEPTH_SCHEDULE
    max_depth isa Integer ||
        throw(ArgumentError("max_depth must be an integer or nothing."))
    Int(max_depth) >= 0 ||
        throw(ArgumentError("max_depth must be nonnegative."))
    return (Int(max_depth),)
end

function _posteriori_storage_mode(store_boxes)
    store_boxes in (:summary, :full) ||
        throw(ArgumentError("store_boxes must be :summary or :full."))
    return store_boxes
end

function _refine_start_for_hc_trace(
    sys::HCSystem,
    x_start,
    t_cht;
    radius,
    tau,
)
    x = sys.CC.(x_start)
    t = sys.CC(t_cht)
    A = compute_preconditioner(sys, x, t)
    x_refined, _, _, success = refine_moore_box(sys, x, t, radius, A; tau = tau)
    return success ? x_refined : x
end

"""
    certify_path_a_posteriori(sys, x_start; kwargs...)

Collect an HC.jl numerical trace for a specialized CHT homotopy and certify that
trace a posteriori. This is the high-level convenience wrapper around
`prepare_posteriori_tracker`, `posteriori_hc_homotopy`, `collect_hc_trace`, and
the trace certifiers.

Public options are intentionally compact:

* `max_step_size = nothing`: first use HC.jl's adaptive step choice, then retry
  failed certifications with smaller trace spacing. Give a positive number, `Inf`,
  or a tuple such as `(Inf, 0.005, 0.0025)` to override this schedule.
* `max_depth = nothing`: use the internal depth schedule `(6, 8, 10, 12)`.
  Give an integer to try only that depth.
* `local_parameter_variables = ()`: optionally restrict adaptive local-parameter
  candidates to selected variables by name or index. The homotopy parameter `t`
  is always kept as a candidate.
* `local_parameter_method = :endpoint_box`: use `:hermite` to compare against the
  Hermite local-parameter tube validator.
* `store_boxes = :summary`: use `:full` only for visualization/debugging.
* `show_progress = false`
* `throw_on_failure = false`
"""
function certify_path_a_posteriori(
    sys::HCSystem,
    x_start;
    max_step_size = nothing,
    max_depth = nothing,
    local_parameter_variables = (),
    local_parameter_method = :endpoint_box,
    store_boxes = :summary,
    show_progress = false,
    throw_on_failure = false,
)
    local_parameter_method in (:endpoint_box, :hermite) ||
        throw(ArgumentError("local_parameter_method must be :endpoint_box or :hermite."))
    tracker = prepare_posteriori_tracker(sys)
    source = tracker.source
    H_hc = posteriori_hc_homotopy(tracker)

    actual_t_start = _posteriori_default_t_start(source)
    actual_t_target = _posteriori_default_t_target(source)
    actual_time_map = _posteriori_default_time_map(source)
    refined = _refine_start_for_hc_trace(
        sys,
        x_start,
        actual_time_map(actual_t_start);
        radius = DEFAULT_POSTERIORI_NODE_REFINEMENT_RADIUS,
        tau = DEFAULT_POSTERIORI_REFINEMENT_TAU,
    )
    actual_x_start = ComplexF64.(refined)
    depth_schedule = _posteriori_depth_schedule(max_depth)

    cert = _certify_hc_path_a_posteriori(
        H_hc,
        sys,
        actual_x_start;
        t_start = actual_t_start,
        t_target = actual_t_target,
        time_map = actual_time_map,
        certification = :adaptive_local_parameter,
        max_step_size_schedule = _posteriori_max_step_schedule(max_step_size),
        max_depth_schedule = depth_schedule,
        local_parameter_method = local_parameter_method,
        midpoint_policy = :krawczyk_polish,
        node_refinement_radius = DEFAULT_POSTERIORI_NODE_REFINEMENT_RADIUS,
        tau = DEFAULT_POSTERIORI_REFINEMENT_TAU,
        radius_growth = DEFAULT_POSTERIORI_RADIUS_GROWTH,
        max_radius = DEFAULT_POSTERIORI_MAX_RADIUS,
        local_parameter_variables = local_parameter_variables,
        store_boxes = _posteriori_storage_mode(store_boxes),
        store_failures = :summary,
        show_progress = show_progress,
    )

    if !cert.success && max_step_size === nothing
        @info "A posteriori certification failed with HC.jl's default trace spacing. Try a smaller max_step_size, for example max_step_size = 0.005 or 0.0025." max_depth_tried = last(depth_schedule)
    end

    if throw_on_failure && !cert.success
        status = haskey(cert, :status) ? cert.status : :certification_failed
        error("a posteriori path certification failed with status $(status)")
    end
    return PosterioriPathResult(
        merge(
            cert,
            (
                posteriori_tracker = tracker,
                t_start = actual_t_start,
                t_target = actual_t_target,
            ),
        ),
    )
end
