include("../src/certified_monodromy_computation.jl")
using Random, Distributions
d = Normal(0, 1)
td = truncated(d, 0.0, Inf)

function dt_bound(H1, x, r, A)
    evalH1 = evaluate_matrix(H1, x)
    rB = CCi("$r+/- $r")
    r_point = map(i -> i + rB, x)
    evalJH1 = evaluate_matrix(jac(H1), r_point)

    norm_evalH1 = max_norm(A*transpose(evalH1))
    norm_evalJH1 = max_norm(A*transpose(evalJH1))
    e = r*norm_evalJH1/norm_evalH1
    return (3/4)*r/(norm_evalH1 + norm_evalJH1*r), e, norm_evalH1
end

function tracking_certified_predictor(H, H1, x, r; iterations_count = false)
    ring = parent(H[1]);
    coeff_ring = base_ring(H[1]);
    t = 0;
    dt = 0;
    h = 1;
    n = length(x);

    G = evaluate_matrix(Matrix(transpose(hcat(H))), 0)
    A = jacobian_inverse(G, x)

    iter = 0;
    r_list = Float64[];
    eta_list = Float64[];
    b_list = Float64[];
    ratio_list = Float64[];
    dt_list = Float64[];
    while t < 1
        rt = round(t, digits = 10)
        print("\r(dt, progress t, radius): ($dt,$rt,$r)")

        Ft = evaluate_matrix(Matrix(transpose(hcat(H))), t);
        x,r,A = refine_moore_box(Ft, x, r, A, 1/8);
        b = max_norm(A*transpose(evaluate_matrix(Ft, x)));


        dt, e, b1 = dt_bound(H1, x, r, A);
    
        t = t+dt;
        iter = iter+1;
        push!(dt_list, dt);
        push!(r_list, r);
        push!(eta_list, e);
        push!(b_list, b/(1+1/8));
        push!(ratio_list, r/(b/(1+1/8)));
    end

    Ft = evaluate_matrix(Matrix(transpose(hcat(H))), 1);
    x,r,A = refine_moore_box(Ft, x, r, A, 1/8);
    if iterations_count
        return x, [iter, minimum(dt_list), median(dt_list),maximum(dt_list),  quantile(dt_list, 0.25), quantile(dt_list, 0.75), minimum(r_list), median(r_list), mean(b_list), maximum(eta_list), mean(ratio_list[2:end])]
    else
        return x
    end
end

function tracking_certified_hermite_predictor(
    H::Union{Matrix, Vector},
    x::Vector{AcbFieldElem},
    r::Number;
    show_display = true,
    refinement_threshold = 1/8,
    predictor = "hermitian",
    iterations_count = false,
    tracking = "truncate"
)

    if predictor == "without_predictor"
        x = CertifiedHomotopyTracking.tracking_without_predictor(H, x; r)
        return x
    end
    ## --- Initialization --------------------------------------------------
    coeff_ring  = base_ring(H[1])
    CCi         = coefficient_ring(coeff_ring)
    t           = 0
    h           = 1 / 10
    iter        = 0
    tprev       = 0
    n           = length(x)

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
        rt = round(t, digits = 10)

        x, v, h, X, r, A = certified_hermite_tracking(
            H, t, x, r, A, h, n, xprev, vprev, hprev, refinement_threshold; tracking
        )

        xprev  = x
        vprev  = v
        hprev  = h
        t     += h

        input        = zeros(CCi, n + 1)
        input[n + 1] = CCi(h)

        x     = midpoint_complex_box(Matrix(evaluate_matrix(matrix(X), input))[1:n])
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
end


function generate_random_dense_system(n::Int, d::Int)
    # Initialize field and rings
    CCi = AcbField()
    
    # Create variable names
    var_names = ["x_$(i)" for i in 1:n]
    push!(var_names, "η")
    
    # Create polynomial ring
    R, vars = CCi[var_names...]
    HR, (t) = R["t"]
    
    # Generate random system
    system = []
    for eq_num in 1:n
        poly = zero(R)
        
        # Generate all possible monomial combinations of degree ≤ d
        for total_deg in 0:d
            # Generate all possible combinations of variables that sum to total_deg
            for part in partitions(total_deg, n)
                coeff = rand(td)  # Random coefficient in [-1,1]
                term = CCi(coeff)
                
                # Multiply by variables with their respective powers
                for (var_idx, power) in enumerate(part)
                    if power > 0
                        term *= vars[var_idx]^power
                    end
                end
                
                poly += term
            end
        end
        
        # Add constant term -1
        poly -= 1
        
        push!(system, poly)
    end
    
    return system
end

# Helper function to generate partitions
function partitions(n::Int, k::Int)
    if k == 1
        return [[n]]
    end
    if n == 0
        return [[0 for _ in 1:k]]
    end
    
    result = []
    for i in 0:n
        for p in partitions(n-i, k-1)
            push!(result, vcat([i], p))
        end
    end
    return result
end
