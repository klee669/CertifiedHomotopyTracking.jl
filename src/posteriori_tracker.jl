export PosterioriTracker, PosterioriPathResult, prepare_posteriori_tracker, posteriori_hc_system,
       posteriori_hc_homotopy, certify_path_a_posteriori

import HomotopyContinuation

_complex_f64_vector(values) = ComplexF64[ComplexF64(value) for value in values]

const _POSTERIORI_PROJECTIVE_COMPILED_CACHE = IdDict{CompiledHomotopy,CompiledHomotopy}()

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

function _posteriori_projective_compiled(compiled::CompiledHomotopy)
    haskey(_POSTERIORI_PROJECTIVE_COMPILED_CACHE, compiled) &&
        return _POSTERIORI_PROJECTIVE_COMPILED_CACHE[compiled]

    source = compiled.source
    source === nothing &&
        throw(ArgumentError("projective posteriori certification requires source metadata."))
    source.projective &&
        throw(ArgumentError("projective posteriori certification expects the original affine compiled homotopy."))
    projective_compiled =
        if source.kind === :edge
            compile_edge_homotopy(
                source.equations,
                source.variables,
                source.parameters;
                projective = true,
                const_vars = source.const_vars,
            )
        elseif source.kind === :direct
            source.t_var === nothing &&
                throw(ArgumentError("direct projective posteriori certification requires a homotopy time variable."))
            compile_edge_homotopy(
                source.equations,
                source.variables,
                Num[source.t_var; source.parameters...; source.const_vars...];
                projective = true,
            )
        else
            throw(ArgumentError("unknown homotopy source kind: $(source.kind)"))
        end
    _POSTERIORI_PROJECTIVE_COMPILED_CACHE[compiled] = projective_compiled
    return projective_compiled
end

function _posteriori_projective_system(sys::HCSystem)
    compiled_projective = _posteriori_projective_compiled(sys.compiled)
    source = sys.source === nothing ? sys.compiled.source : sys.source
    if source.kind === :edge
        return make_edge_system(
            compiled_projective,
            AcbFieldElem[sys.p_start...],
            AcbFieldElem[sys.p_end...],
            AcbFieldElem[sys.p_const...],
        )
    elseif source.kind === :direct
        values = AcbFieldElem[sys.p_start...; sys.p_const...]
        expected = length(source.parameters) + length(source.const_vars)
        length(values) == expected ||
            throw(DimensionMismatch("direct projective posteriori expected $(expected) fixed parameters, got $(length(values))."))
        return make_edge_system(
            compiled_projective,
            AcbFieldElem[sys.CC(0); values],
            AcbFieldElem[sys.CC(1); values],
        )
    else
        throw(ArgumentError("unknown homotopy source kind: $(source.kind)"))
    end
end

function _posteriori_affine_to_projective_chart(x_affine, chart_idx::Int)
    X = ComplexF64[1.0 + 0im; ComplexF64.(x_affine)]
    1 <= chart_idx <= length(X) ||
        throw(ArgumentError("projective chart index $(chart_idx) is outside 1:$(length(X))."))
    scale = X[chart_idx]
    abs(scale) > 0 ||
        throw(ArgumentError("cannot use projective chart $(chart_idx): scale is numerically zero."))
    return ComplexF64[X[i] / scale for i in eachindex(X) if i != chart_idx]
end

function _posteriori_projective_chart_score(trace, chart_idx::Int)
    best = Inf
    for point in trace
        X = ComplexF64[1.0 + 0im; point.x]
        best = min(best, abs(X[chart_idx]))
    end
    return best
end

function _posteriori_select_projective_chart(trace, chart)
    isempty(trace) && throw(ArgumentError("cannot choose a projective chart from an empty trace."))
    n_projective = length(first(trace).x) + 1
    if chart === nothing || chart === :auto
        _, idx = findmax([_posteriori_projective_chart_score(trace, i) for i in 1:n_projective])
        return idx
    end
    chart isa Integer ||
        throw(ArgumentError("projective_chart must be an integer, :auto, or nothing."))
    1 <= Int(chart) <= n_projective ||
        throw(ArgumentError("projective_chart $(chart) is outside 1:$(n_projective)."))
    return Int(chart)
end

function _posteriori_trace_to_projective_chart(trace, chart_idx::Int)
    return [
        (t = point.t, x = _posteriori_affine_to_projective_chart(point.x, chart_idx))
        for point in trace
    ]
end

function _certify_path_a_posteriori_projective(
    H_hc,
    sys::HCSystem,
    x_start;
    time_map,
    max_step_size_schedule,
    max_depth_schedule,
    projective_chart = :auto,
    local_parameter_method,
    local_parameter_variables,
    store_boxes,
    diagnostics,
    show_progress,
)
    source = sys.source === nothing ? sys.compiled.source : sys.source
    source === nothing &&
        throw(ArgumentError("projective posteriori HC tracing requires homotopy source metadata."))

    sys_projective = _posteriori_projective_system(sys)
    t_start = _posteriori_default_t_start(source)
    t_target = _posteriori_default_t_target(source)

    attempts = NamedTuple[]
    last_cert = nothing
    last_trace = nothing
    last_chart = missing

    for max_step_size in max_step_size_schedule
        trace_out = collect_hc_trace(
            H_hc,
            x_start;
            t_start = t_start,
            t_target = t_target,
            max_step_size = max_step_size,
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
                    projective_chart = missing,
                ),
            )
            continue
        end

        chart_idx = _posteriori_select_projective_chart(trace_out.trace, projective_chart)
        last_chart = chart_idx
        sys_projective.patch_idx = chart_idx
        projective_trace = _posteriori_trace_to_projective_chart(trace_out.trace, chart_idx)

        for depth_value in max_depth_schedule
            cert = try
                certify_hc_trace_adaptive_local_parameter(
                    sys_projective,
                    projective_trace;
                    time_map = time_map,
                    max_depth = depth_value,
                    local_parameter_method = local_parameter_method,
                    midpoint_policy = :krawczyk_polish,
                    node_refinement_radius = DEFAULT_POSTERIORI_NODE_REFINEMENT_RADIUS,
                    tau = DEFAULT_POSTERIORI_REFINEMENT_TAU,
                    radius_growth = DEFAULT_POSTERIORI_RADIUS_GROWTH,
                    max_radius = DEFAULT_POSTERIORI_MAX_RADIUS,
                    local_parameter_variables = local_parameter_variables,
                    store_boxes = store_boxes,
                    store_failures = :summary,
                    diagnostics = diagnostics,
                    show_progress = show_progress,
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
                    original_segments = max(length(projective_trace) - 1, 0),
                    max_depth = depth_value,
                    max_krawczyk_norm = Inf,
                    method = :adaptive_local_parameter_projective,
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
                    projective_chart = chart_idx,
                ),
            )
            cert.success && return merge(
                cert,
                (
                    hc_trace = trace_out,
                    affine_endpoint = sys.CC.(last(trace_out.trace).x),
                    max_step_size = max_step_size,
                    attempts = attempts,
                    certification = :adaptive_local_parameter_projective,
                    certification_chart = :projective,
                    projective_chart = chart_idx,
                    projective_system = sys_projective,
                ),
            )
        end
    end

    last_cert === nothing && return (
        success = false,
        hc_trace = last_trace,
        affine_endpoint = last_trace === nothing || isempty(last_trace.trace) ? AcbFieldElem[] : sys.CC.(last(last_trace.trace).x),
        max_step_size = isempty(attempts) ? missing : last(attempts).max_step_size,
        attempts = attempts,
        certification = :adaptive_local_parameter_projective,
        certification_chart = :projective,
        projective_chart = last_chart,
        status = :hc_failed,
        total_boxes = 0,
        original_segments = 0,
        max_depth = 0,
        failed_segments = NamedTuple[],
    )

    return merge(
        last_cert,
        (
            hc_trace = last_trace,
            affine_endpoint = last_trace === nothing || isempty(last_trace.trace) ? AcbFieldElem[] : sys.CC.(last(last_trace.trace).x),
            max_step_size = isempty(attempts) ? missing : last(attempts).max_step_size,
            attempts = attempts,
            certification = :adaptive_local_parameter_projective,
            certification_chart = :projective,
            projective_chart = last_chart,
            projective_system = sys_projective,
        ),
    )
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
* `certification_chart = :affine`: use `:projective` to certify the same HC trace
  in a projective chart. Use `:auto` to try affine certification first and retry
  in a projective chart only if affine certification fails. With
  `projective_chart = :auto`, a projective chart is chosen from the trace.
* `store_boxes = :summary`: use `:full` only for visualization/debugging.
* `diagnostics = :off`: use `:basic` or `:timing` to record certification
  counters/diagnostic summaries.
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
    certification_chart = :affine,
    projective_chart = :auto,
    store_boxes = :summary,
    diagnostics = :off,
    show_progress = false,
    throw_on_failure = false,
)
    local_parameter_method in (:endpoint_box, :hermite) ||
        throw(ArgumentError("local_parameter_method must be :endpoint_box or :hermite."))
    certification_chart in (:affine, :projective, :auto) ||
        throw(ArgumentError("certification_chart must be :affine, :projective, or :auto."))
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
    storage_mode = _posteriori_storage_mode(store_boxes)

    if certification_chart === :projective
        cert = _certify_path_a_posteriori_projective(
            H_hc,
            sys,
            actual_x_start;
            time_map = actual_time_map,
            max_step_size_schedule = _posteriori_max_step_schedule(max_step_size),
            max_depth_schedule = depth_schedule,
            projective_chart = projective_chart,
            local_parameter_method = local_parameter_method,
            local_parameter_variables = local_parameter_variables,
            store_boxes = storage_mode,
            diagnostics = diagnostics,
            show_progress = show_progress,
        )
        if throw_on_failure && !cert.success
            status = haskey(cert, :status) ? cert.status : :certification_failed
            error("projective a posteriori path certification failed with status $(status)")
        end
        return PosterioriPathResult(
            merge(
                cert,
                (
                    posteriori_tracker = tracker,
                    certification_chart = :projective,
                    t_start = actual_t_start,
                    t_target = actual_t_target,
                ),
            ),
        )
    end

    cert = _certify_hc_path_a_posteriori(
        H_hc,
        sys,
        actual_x_start;
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
        store_boxes = storage_mode,
        store_failures = :summary,
        diagnostics = diagnostics,
        show_progress = show_progress,
    )

    if cert.success || certification_chart === :affine
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
                    certification_chart = :affine,
                    t_start = actual_t_start,
                    t_target = actual_t_target,
                ),
            ),
        )
    end

    projective_cert = _certify_path_a_posteriori_projective(
        H_hc,
        sys,
        actual_x_start;
        time_map = actual_time_map,
        max_step_size_schedule = _posteriori_max_step_schedule(max_step_size),
        max_depth_schedule = depth_schedule,
        projective_chart = projective_chart,
        local_parameter_method = local_parameter_method,
        local_parameter_variables = local_parameter_variables,
        store_boxes = storage_mode,
        diagnostics = diagnostics,
        show_progress = show_progress,
    )

    if !projective_cert.success && max_step_size === nothing
        @info "A posteriori certification failed with HC.jl's default trace spacing. Try a smaller max_step_size, for example max_step_size = 0.005 or 0.0025." max_depth_tried = last(depth_schedule)
    end

    if throw_on_failure && !projective_cert.success
        status = haskey(projective_cert, :status) ? projective_cert.status : :certification_failed
        error("a posteriori path certification failed with status $(status)")
    end
    return PosterioriPathResult(
        merge(
            projective_cert,
            (
                posteriori_tracker = tracker,
                certification_chart = :projective,
                affine_cert = cert,
                t_start = actual_t_start,
                t_target = actual_t_target,
            ),
        ),
    )
end
