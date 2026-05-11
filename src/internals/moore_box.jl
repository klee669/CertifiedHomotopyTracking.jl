export refine_moore_box

# refining the Moore box
function refine_moore_box(
    system::Union{Matrix,AbstractAlgebra.Generic.MatSpaceElem}, 
    point::Vector{AcbFieldElem}, 
    r::Number,
    A::AcbMatrix,
    ρ::Number
)
    y = point;
    CR = parent(point[1]);
    n = size(A)[1];
    while krawczyk_test(system, y, r, A, ρ) == false
        d = A * transpose(evaluate_matrix(system, y));
        if max_norm(d) <= (1/64)*ρ*r
            r = (1/2)*r;
        else
            y = midpoint_complex_box(y-d[:,1]);
        end
        A = jacobian_inverse(system, y);
    end
    while 2*r <= 1 && krawczyk_test(system, point, 2*r, A, ρ)
        r = 2*r;
    end

    [y, r, A]
end


function refine_moore_box(sys::HCSystem, x::AbstractVector{AcbFieldElem}, t, r_init, A_init; tau=0.125)
    RR = sys.RR
    y = copy(x)
    s = RR(r_init)
    U = copy(A_init)
    max_iter = 20
    
    for iter in 1:max_iter
        passed, k_norm = krawczyk_test(sys, y, t, s, U; rho=tau)
        if passed
            while 2*s <= 1.0
                p2, _ = krawczyk_test(sys, y, t, 2*s, U; rho=tau)
                if p2 s *= 2 else break end
            end
            return y, s, U, true
        end
        
        fy = evaluate_H(sys, y, t)
        delta = U * fy
        
        if norm_inf(delta) <= (1/64) * tau * Float64(s)
            s /= 2
        else
            y_next = y - delta
            y = get_mid_vec(y_next)

            Jy = evaluate_Jac(sys, y, t)
            U = inv_acb(Jy)
        end
    end
    return y, s, U, false
end
