export track, tracking_without_predictor, track_path

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
    tracking = "truncate",
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


function track_path(sys::HCSystem, x_start_input::Vector{AcbFieldElem}; t_start=0.0, t_end=1.0, h_init=0.1)
    t = RR(t_start)
    t_target = RR(t_end)
    h = RR(h_init)
    
    if sys.homogeneous
        if length(x_start_input) == length(sys.vars) - 1
            x = [CC(1); x_start_input]
        else
            x = copy(x_start_input)
        end
        
        mags = [mag_complex(xi) for xi in x]
        max_val, max_idx = findmax(mags)
        
        scale = x[max_idx]
        x = x ./ scale
        sys.patch_idx = max_idx
    else
        x = copy(x_start_input)
    end
    
    A = compute_preconditioner(sys, x, t)
    x, r, A, success = refine_moore_box(sys, x, t, 1e-6, A)
    if !success 
        return x, false 
    end
    
    x_prev = copy(x)
    v_prev = compute_velocity(sys, x, t, A)
    h_prev = h
    
    iter = 0
    
    while t < t_target
        iter += 1
        dt_remaining = t_target - t
        if h > dt_remaining h = dt_remaining end
        
        if sys.homogeneous
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
        if !success return x, false end

        v = compute_velocity(sys, x, t, A)
        
        step_accepted = false
        min_h = RR(1e-20)
        
        while h > min_h
            local X_tm
            if iter == 1
                X_tm = [TaylorModel3(x[i], v[i], CC(0), CC(0), CC(0), h) for i in 1:length(x)]
            else
                X_tm = construct_hermite_predictor_tm(x, x_prev, v, v_prev, h_prev, h)
            end
            
            passed, k_norm = validate_step_taylor3(sys, X_tm, t, h, Float64(r), A)
            
            if passed
                step_accepted = true
                x_next_interval = evaluate_taylor.(X_tm)
                x_new = get_mid_vec(x_next_interval)
                x_prev = x; v_prev = v; h_prev = h
                x = x_new; t += h
                
                print("Iter $iter: t=$(Float64(t)), h=$(Float64(h)) \r")
                h = min(h * 2, RR(0.5))
                break
            else
                h /= 2
            end
        end
        if !step_accepted 
            @info("too small step size!")
            return x, false 
        end
    end
    
    print("\r" * " "^60 * "\r") 

    if sys.homogeneous
        try
            A_final = compute_preconditioner(sys, x, t_target)
            x_polished, _, _, success_polish = refine_moore_box(sys, x, t_target, 1e-8, A_final)
            if success_polish
                x = x_polished
            end
        catch
        end

        return x, true
    else
        A_final = compute_preconditioner(sys, x, t_target)
        x_polished, _, _, success_polish = refine_moore_box(sys, x, t_target, 1e-8, A_final)
        if success_polish
            x = x_polished
        end
        return x, true
    end
end
