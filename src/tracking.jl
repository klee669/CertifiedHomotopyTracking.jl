export track, tracking_without_predictor, track_path

function _print_track_progress(iter, t, h, precision_bits)
    msg = "Iter $iter: t=$(round(Float64(t); sigdigits=12)), h=$(round(Float64(h); sigdigits=12)), precision=$(precision_bits)"
    print("\r", msg, "\033[K")
    flush(stdout)
end

function _clear_track_progress()
    print("\r\033[K")
    flush(stdout)
end

# tracking without predictor
function tracking_without_predictor(H, x; r = .1, iterations_count = false)
    try
        ring = parent(H[1]);
        coeff_ring = base_ring(H[1]);
        CCi = coefficient_ring(coeff_ring)
        t = 0;
        h = 1;
        n = length(x);
        G = evaluate_matrix(Matrix(transpose(hcat(H))), 0)
        A = jacobian_inverse(G, x)

        iter = 0;
        while t < 1
            yield()
            
            rt = round(t, digits = 10)
            print("\r(h, progress t): ($h,$rt)")

            Ft = evaluate_matrix(Matrix(transpose(hcat(H))), t);
            x,r,A = refine_moore_box(Ft, x, r, A, 1/8);
            h = 2*h;
            radii = h/2;


            midt = t+h/2;
            T = CCi("$midt +/- $radii");
            FT = evaluate_matrix(Matrix(transpose(hcat(H))), t+h);
            while krawczyk_test(FT, x, r, A, 7/8) == false
                h = 1/2 * h;
                midt = t+h/2;
                radii = h/2;
        
                T = CCi("$midt +/- $radii");
                FT = evaluate_matrix(Matrix(transpose(hcat(H))), t+h);
            end
            t = max_int_norm(T);
            iter = iter+1;
        end

        Ft = evaluate_matrix(Matrix(transpose(hcat(H))), 1);
        x,r,A = refine_moore_box(Ft, x, r, A, 1/8);
        
        if iterations_count
            return x, iter
        else
            return x
        end

    catch e
        if isa(e, InterruptException)
            println("\n[Warning] Tracking interrupted manually.")
            return nothing
        else
            rethrow(e)
        end
    end
end


"""
    track(H, x, r; options...)

Track a solution path. Returns `nothing` if interrupted.
"""
function track(
    H::Union{Matrix, Vector},
    x::Vector{AcbFieldElem};
    r = 0.1, # initial radius
    show_display = true,
    refinement_threshold = 1/8,
    predictor = "hermitian",
    iterations_count = false,
    tracking = "non-truncate",
    projective = false,
)
    try
        if predictor == "without_predictor"
            return tracking_without_predictor(H, x)
        end

        ## --- Initialization --------------------------------------------------
        coeff_ring  = base_ring(H[1])
        CCi         = coefficient_ring(coeff_ring)
        t           = 0
        h           = 1 / 10
        iter        = 0
        tprev       = 0
        n           = length(x)
        max_deg_H   = max_degree(hcat(H))

        G = evaluate_matrix(Matrix(transpose(hcat(H))), 0)
        A = jacobian_inverse(G, x)

        ## --- First Step (Linear predictor) ----------------------------------
        x, v, h, X, r, A = linear_tracking(H, t, x, r, A, h, n, refinement_threshold)

        xprev  = x
        vprev  = v
        hprev  = h
        t     += h

        ## --- Main Loop -------------------------------------------------------
        while t < 1
            yield()
            rt = round(t, digits = 10)

            H, x, v, h, X, r, A = hermite_tracking(
                H, max_deg_H, t, x, r, A, h, n, xprev, vprev, hprev, refinement_threshold;
                tracking, projective
            )

            xprev  = x
            vprev  = v
            hprev  = h
            t     += h

            input        = zeros(CCi, n + 1)
            input[n + 1] = CCi(h)

            x = midpoint_complex_box(Matrix(evaluate_matrix(matrix(X), input))[1:n])
            iter += 1

            show_display && print("\rprogress t: $rt")
        end

        show_display && print("\r ")

        ## --- Final Refinement ------------------------------------------------
        Ft, x, r, A, v, h, radii = refine_step(H, 1, x, r, A, h; threshold = 1 / 100)
        
        if iterations_count
            return x, iter
        else
            return x
        end

    catch e
        if isa(e, InterruptException)
            println("\n[Warning] Tracking interrupted manually.")
            return nothing 
        else
            rethrow(e)
        end
    end
end

function _projective_coordinate_from_chart(x, chart_idx::Int, coord_idx::Int)
    if coord_idx == chart_idx
        return one(x[1])
    elseif coord_idx < chart_idx
        return x[coord_idx]
    else
        return x[coord_idx - 1]
    end
end

function _chart_to_projective_coordinates(sys::HCSystem, x, chart_idx::Int=sys.patch_idx)
    X = Vector{eltype(x)}(undef, sys.compiled.n_vars)
    for i in eachindex(X)
        X[i] = _projective_coordinate_from_chart(x, chart_idx, i)
    end
    return X
end

function _projective_to_chart_coordinates(sys::HCSystem, X, chart_idx::Int)
    scale = X[chart_idx]
    if mag_complex(scale) <= 1e-30
        throw(ArgumentError("Cannot move to projective chart $chart_idx because its coordinate is numerically zero."))
    end
    x = Vector{eltype(X)}(undef, length(X) - 1)
    out_idx = 1
    for i in eachindex(X)
        if i != chart_idx
            x[out_idx] = X[i] / scale
            out_idx += 1
        end
    end
    return x
end

function _projective_to_chart_coordinates(sys::HCSystem, X)
    mags = [mag_complex(xi) for xi in X]
    _, chart_idx = findmax(mags)
    return _projective_to_chart_coordinates(sys, X, chart_idx), chart_idx
end

function _canonical_projective_chart(sys::HCSystem, x)
    X = _chart_to_projective_coordinates(sys, x)
    mags = [mag_complex(xi) for xi in X]
    max_val, chart_idx = findmax(mags)
    return max_val > 1.5 ? chart_idx : sys.patch_idx
end

function _projective_chart_box(sys::HCSystem, x, chart_idx::Int, r)
    CC = sys.CC
    unit = sys.RR("0 +/- 1")
    unit_box = CC(unit, unit)
    x_box = x .+ unit_box .* CC(r)
    return _chart_to_projective_coordinates(sys, x_box, chart_idx)
end

function _chart_box_center_radius(y_box; inflate=1.1, min_radius=1e-30)
    center = get_mid_vec(y_box)
    radius = min_radius
    for i in eachindex(y_box)
        radius = max(radius, inflate * max_int_norm(y_box[i] - center[i]))
    end
    return center, radius
end

function _certified_projective_chart_transfer(
    sys::HCSystem,
    x,
    t,
    r,
    target_chart_idx::Int;
    rho=0.7,
)
    old_chart_idx = sys.patch_idx
    target_chart_idx == old_chart_idx && return x, r, nothing, true, 0.0

    X_box = _projective_chart_box(sys, x, old_chart_idx, r)
    denominator = X_box[target_chart_idx]
    contains(denominator, 0) && return x, r, nothing, false, Inf

    y_box = Vector{AcbFieldElem}(undef, length(X_box) - 1)
    out_idx = 1
    for i in eachindex(X_box)
        if i != target_chart_idx
            y_box[out_idx] = X_box[i] / denominator
            out_idx += 1
        end
    end

    y, r_y = _chart_box_center_radius(y_box)
    sys.patch_idx = target_chart_idx
    A_y = compute_preconditioner(sys, y, t)
    passed, k_norm = krawczyk_test(sys, y, t, r_y, A_y; rho=rho)
    if passed
        return y, r_y, A_y, true, k_norm
    end
    sys.patch_idx = old_chart_idx
    return x, r, nothing, false, k_norm
end


function _track_path_at_precision(
    sys::HCSystem,
    x_start_input::Vector{AcbFieldElem};
    t_start=0.0,
    t_end=1.0,
    h_init=0.1,
    rho=0.7,
    show_progress=false,
    projective=false,
    patch_vector=nothing,
    affine_chart_atol=1e-10,
    request_precision_retry=false,
    precision_rejection_threshold=8,
    initial_precision=precision(sys.CC),
)
    CC = sys.CC; RR = sys.RR
    t = RR(t_start)
    t_target = RR(t_end)
    h = RR(h_init)
    r = 1e-6
    iter = 0
    accepted_steps = 0
    rejected_steps = 0
    last_k_norm = Inf
    user_input_start = copy(x_start_input)

    patch_vector === nothing || throw(ArgumentError("patch_vector is no longer supported; projective tracking uses coordinate charts."))
    if projective && !uses_projective_charts(sys)
        throw(ArgumentError("projective=true requires a system compiled with projective=true."))
    end

    if uses_projective_charts(sys)
        dim = sys.compiled.n_vars - 1
        if length(x_start_input) == dim
            x = copy(x_start_input)
            sys.patch_idx = 1
        elseif length(x_start_input) == sys.compiled.n_vars
            x, chart_idx = _projective_to_chart_coordinates(sys, x_start_input)
            sys.patch_idx = chart_idx
        else
            throw(DimensionMismatch("x_start_input must have affine chart length n or projective-coordinate length n + 1."))
        end
    elseif sys.projective_coordinates
        if sys.compiled.n_vars != 0 && length(x_start_input) == sys.compiled.n_vars - 1
            x = [CC(1); x_start_input]
        elseif sys.compiled.n_vars == 0 || length(x_start_input) == sys.compiled.n_vars
            x = copy(x_start_input)
        else
            throw(DimensionMismatch("x_start_input must have affine length n or projective-coordinate length n + 1."))
        end

        mags = [mag_complex(xi) for xi in x]
        max_val, max_idx = findmax(mags)
        scale = x[max_idx]
        x = x ./ scale
        sys.patch_idx = max_idx
    else
        x = copy(x_start_input)
    end
    tracking_start_patch_idx = sys.patch_idx
    tracking_input_start = copy(x)
    tm_cache = TMCache(CC)
    validation_cache = KrawczykValidationCache(CC, RR, length(x))
    x_next_interval = [CC(0) for _ in 1:length(x)]

    A = compute_preconditioner(sys, x, t)
    x, r, A, success = refine_moore_box(sys, x, t, r, A)
    tracking_refined_start = copy(x)
    if !success
        return _track_result(
            sys,
            x,
            false,
            :initial_refinement_failed,
            iter,
            accepted_steps,
            rejected_steps,
            t,
            h,
            r,
            last_k_norm,
            rho,
            "Initial Moore box refinement failed.";
            affine_chart_atol=affine_chart_atol,
            input_start=user_input_start,
            tracking_input_start=tracking_input_start,
            tracking_refined_start=tracking_refined_start,
            tracking_start_patch_idx=tracking_start_patch_idx,
            initial_precision=initial_precision,
        )
    end

    x_prev = copy(x)
    v_prev = compute_velocity(sys, x, t, A)
    h_prev = h

    while t < t_target
        iter += 1
        dt_remaining = t_target - t
        if h > dt_remaining h = dt_remaining end

        if uses_projective_charts(sys)
            chart_idx = _canonical_projective_chart(sys, x)
            if chart_idx != sys.patch_idx
                x_chart, r_chart, A_chart, transfer_success, transfer_k_norm =
                    _certified_projective_chart_transfer(sys, x, t, r, chart_idx; rho=rho)
                last_k_norm = transfer_k_norm
                if transfer_success
                    x = x_chart
                    r = r_chart
                    A = A_chart
                    x_prev = copy(x)
                    v_prev = compute_velocity(sys, x, t, A)
                    h_prev = h
                end
            end
        elseif sys.projective_coordinates && !has_projective_patch(sys)
            mags = [mag_complex(xi) for xi in x]
            max_val, max_idx = findmax(mags)
            if max_idx != sys.patch_idx && max_val > 1.5
                scale = x[max_idx]
                x = x ./ scale
                x_prev = x_prev ./ scale
                v_prev = v_prev ./ scale
                sys.patch_idx = max_idx
                A = compute_preconditioner(sys, x, t)
            end
        end

        x, r, A, success = refine_moore_box(sys, x, t, r, A)
        if !success
            return _track_result(
                sys,
                x,
                false,
                :refinement_failed,
                iter,
                accepted_steps,
                rejected_steps,
                t,
                h,
                r,
                last_k_norm,
                rho,
                "Moore box refinement failed during tracking.";
                affine_chart_atol=affine_chart_atol,
                input_start=user_input_start,
                tracking_input_start=tracking_input_start,
                tracking_refined_start=tracking_refined_start,
                tracking_start_patch_idx=tracking_start_patch_idx,
                initial_precision=initial_precision,
            )
        end

        v = compute_velocity(sys, x, t, A)
        step_accepted = false
        min_h = RR(1e-20)
        previous_rejected_k_norm = Inf
        stagnant_rejections = 0

        while h > min_h
            local X_tm
            if iter == 1
                X_tm = [TaylorModel3(x[i], v[i], CC(0), CC(0), CC(0), RR(h)) for i in 1:length(x)]
            else
                X_tm = construct_hermite_predictor_tm(sys, x, x_prev, v, v_prev, h_prev, h)
            end

            passed, k_norm = validate_step_taylor3(sys, X_tm, t, h, Float64(r), A; rho=rho, cache=validation_cache)
            last_k_norm = k_norm

            if passed
                step_accepted = true
                accepted_steps += 1
                for i in eachindex(X_tm)
                    evaluate_taylor!(x_next_interval[i], X_tm[i], tm_cache)
                end
                x_new = get_mid_vec(x_next_interval)
                x_prev = x; v_prev = v; h_prev = h
                x = x_new; t += h

                if uses_projective_charts(sys)
                    chart_idx = _canonical_projective_chart(sys, x)
                    if chart_idx != sys.patch_idx
                        x_chart, r_chart, A_chart, transfer_success, transfer_k_norm =
                            _certified_projective_chart_transfer(sys, x, t, r, chart_idx; rho=rho)
                        last_k_norm = transfer_k_norm
                        if transfer_success
                            x = x_chart
                            r = r_chart
                            A = A_chart
                            x_prev = copy(x)
                            v_prev = compute_velocity(sys, x, t, A)
                            h_prev = h
                        end
                    end
                end

                show_progress && _print_track_progress(iter, t, h, precision(sys.CC))
                h = min(h * 2, RR(0.5))
                break
            else
                rejected_steps += 1
                if !isfinite(k_norm) || k_norm >= 0.9 * previous_rejected_k_norm
                    stagnant_rejections += 1
                else
                    stagnant_rejections = 0
                end
                previous_rejected_k_norm = k_norm
                h /= 2
                if request_precision_retry && stagnant_rejections >= precision_rejection_threshold
                    show_progress && _clear_track_progress()
                    return _track_result(
                        sys,
                        x,
                        false,
                        :precision_increase_required,
                        iter,
                        accepted_steps,
                        rejected_steps,
                        t,
                        h,
                        r,
                        last_k_norm,
                        rho,
                        "Repeated certified step validation failures suggest insufficient arithmetic precision.";
                        affine_chart_atol=affine_chart_atol,
                        input_start=user_input_start,
                        tracking_input_start=tracking_input_start,
                        tracking_refined_start=tracking_refined_start,
                        tracking_start_patch_idx=tracking_start_patch_idx,
                        initial_precision=initial_precision,
                    )
                end
            end
        end

        if !step_accepted
            show_progress && _clear_track_progress()
            return _track_result(
                sys,
                x,
                false,
                :step_too_small,
                iter,
                accepted_steps,
                rejected_steps,
                t,
                h,
                r,
                last_k_norm,
                rho,
                "Step size became too small before validation succeeded.";
                affine_chart_atol=affine_chart_atol,
                input_start=user_input_start,
                tracking_input_start=tracking_input_start,
                tracking_refined_start=tracking_refined_start,
                tracking_start_patch_idx=tracking_start_patch_idx,
                initial_precision=initial_precision,
            )
        end
    end

    show_progress && _clear_track_progress()

    if uses_projective_charts(sys)
        chart_idx = _canonical_projective_chart(sys, x)
        if chart_idx != sys.patch_idx
            x_chart, r_chart, A_chart, transfer_success, transfer_k_norm =
                _certified_projective_chart_transfer(sys, x, t_target, r, chart_idx; rho=rho)
            last_k_norm = transfer_k_norm
            if transfer_success
                x = x_chart
                r = r_chart
                A = A_chart
            end
        end
    end

    final_status = :success
    final_message = "Path tracked successfully."
    final_success = true

    if sys.projective_coordinates
        try
            A_final = compute_preconditioner(sys, x, t_target)
            x_polished, r_polished, _, success_polish = refine_moore_box(sys, x, t_target, 1e-8, A_final)
            if success_polish
                x = x_polished
                r = r_polished
            else
                final_status = :final_refinement_failed
                final_message = "Path reached target, but final polishing failed."
                final_success = false
            end
        catch err
            final_status = :final_refinement_error
            final_message = "Path reached target, but final polishing raised an error: $(sprint(showerror, err))"
            final_success = false
        end
    else
        A_final = compute_preconditioner(sys, x, t_target)
        x_polished, r_polished, _, success_polish = refine_moore_box(sys, x, t_target, 1e-8, A_final)
        if success_polish
            x = x_polished
            r = r_polished
        else
            final_status = :final_refinement_failed
            final_message = "Path reached target, but final polishing failed."
            final_success = false
        end
    end

    if final_success && uses_projective_charts(sys)
        chart_idx = _canonical_projective_chart(sys, x)
        if chart_idx != sys.patch_idx
            x_chart, r_chart, A_chart, transfer_success, transfer_k_norm =
                _certified_projective_chart_transfer(sys, x, t_target, r, chart_idx; rho=rho)
            last_k_norm = transfer_k_norm
            if transfer_success
                x = x_chart
                r = r_chart
                x_polished, r_polished, _, success_polish = refine_moore_box(sys, x, t_target, r, A_chart)
                if success_polish
                    x = x_polished
                    r = r_polished
                end
            end
        end
    end

    return _track_result(
        sys,
        x,
        final_success,
        final_status,
        iter,
        accepted_steps,
        rejected_steps,
        t,
        h,
        r,
        last_k_norm,
        rho,
        final_message;
        affine_chart_atol=affine_chart_atol,
        input_start=user_input_start,
        tracking_input_start=tracking_input_start,
        tracking_refined_start=tracking_refined_start,
        tracking_start_patch_idx=tracking_start_patch_idx,
        initial_precision=initial_precision,
    )
end

function track_path(
    sys::HCSystem,
    x_start_input::Vector{AcbFieldElem};
    adaptive_precision=true,
    min_precision=53,
    max_precision=max(100_000, precision(sys.CC)),
    precision_rejection_threshold=8,
    kwargs...,
)
    min_precision >= 53 || throw(ArgumentError("min_precision must be at least 53 bits."))
    max_precision >= min_precision || throw(ArgumentError("max_precision must be at least min_precision."))
    precision_rejection_threshold > 0 || throw(ArgumentError("precision_rejection_threshold must be positive."))

    if !adaptive_precision
        return _track_path_at_precision(sys, x_start_input; kwargs...)
    end

    current_precision = min_precision
    retryable_statuses = (
        :initial_refinement_failed,
        :refinement_failed,
        :precision_increase_required,
        :step_too_small,
        :final_refinement_failed,
        :final_refinement_error,
    )

    while true
        current_sys = system_with_precision(sys, current_precision)
        current_start = current_sys.CC.(x_start_input)
        result = _track_path_at_precision(
            current_sys,
            current_start;
            request_precision_retry=current_precision < max_precision,
            precision_rejection_threshold=precision_rejection_threshold,
            initial_precision=min_precision,
            kwargs...,
        )
        if succeeded(result) || !(result.status in retryable_statuses) || current_precision >= max_precision
            return _result_with_field(result, sys.CC)
        end

        next_precision = min(max_precision, 2 * current_precision)
        @info "Retrying certified path tracking with higher precision" status=result.status precision=current_precision next_precision
        current_precision = next_precision
    end
end
