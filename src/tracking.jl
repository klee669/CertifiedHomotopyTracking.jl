export track, tracking_without_predictor, track_path

function _print_track_progress(iter, t, h)
    msg = "Iter $iter: t=$(round(Float64(t); sigdigits=12)), h=$(round(Float64(h); sigdigits=12))"
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


function track_path(
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

    if projective
        sys.homogeneous || throw(ArgumentError("projective=true requires a homogenized/projective HCSystem. Construct it with projective=true."))
        if patch_vector !== nothing
            length(patch_vector) == sys.compiled.n_vars || throw(DimensionMismatch("patch_vector length must match the number of projective variables."))
            sys.patch_vector = Tuple(CC(a) for a in patch_vector)
        elseif !has_projective_patch(sys)
            sys.compiled.n_vars > 0 || throw(ArgumentError("Cannot create a projective patch because the compiled variable count is unknown."))
            sys.patch_vector = Tuple(CC(a) for a in random_patch_vector(sys.compiled.n_vars))
        end
    end

    if has_projective_patch(sys)
        if length(x_start_input) == length(sys.patch_vector) - 1
            x = lift_to_patch(x_start_input, collect(sys.patch_vector))
        elseif length(x_start_input) == length(sys.patch_vector)
            x = repatch(x_start_input, collect(sys.patch_vector))
        else
            throw(DimensionMismatch("x_start_input must have affine length n or projective length n + 1."))
        end
        sys.patch_idx = 1
    elseif sys.homogeneous
        if sys.compiled.n_vars != 0 && length(x_start_input) == sys.compiled.n_vars - 1
            x = [CC(1); x_start_input]
        elseif sys.compiled.n_vars == 0 || length(x_start_input) == sys.compiled.n_vars
            x = copy(x_start_input)
        else
            throw(DimensionMismatch("x_start_input must have affine length n or homogeneous length n + 1."))
        end

        mags = [mag_complex(xi) for xi in x]
        max_val, max_idx = findmax(mags)
        scale = x[max_idx]
        x = x ./ scale
        sys.patch_idx = max_idx
    else
        x = copy(x_start_input)
    end
    tracking_input_start = copy(x)

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
        )
    end

    x_prev = copy(x)
    v_prev = compute_velocity(sys, x, t, A)
    h_prev = h

    while t < t_target
        iter += 1
        dt_remaining = t_target - t
        if h > dt_remaining h = dt_remaining end

        if sys.homogeneous && !has_projective_patch(sys)
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
            )
        end

        v = compute_velocity(sys, x, t, A)
        step_accepted = false
        min_h = RR(1e-20)

        while h > min_h
            local X_tm
            if iter == 1
                X_tm = [TaylorModel3(x[i], v[i], CC(0), CC(0), CC(0), RR(h)) for i in 1:length(x)]
            else
                X_tm = construct_hermite_predictor_tm(sys, x, x_prev, v, v_prev, h_prev, h)
            end

            passed, k_norm = validate_step_taylor3(sys, X_tm, t, h, Float64(r), A; rho=rho)
            last_k_norm = k_norm

            if passed
                step_accepted = true
                accepted_steps += 1
                cache = TMCache(CC)
                x_next_interval = [evaluate_taylor!(CC(0), tm, cache) for tm in X_tm]
                x_new = get_mid_vec(x_next_interval)
                x_prev = x; v_prev = v; h_prev = h
                x = x_new; t += h

                show_progress && _print_track_progress(iter, t, h)
                h = min(h * 2, RR(0.5))
                break
            else
                rejected_steps += 1
                h /= 2
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
            )
        end
    end

    show_progress && _clear_track_progress()

    final_status = :success
    final_message = "Path tracked successfully."
    final_success = true

    if sys.homogeneous
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
    )
end
