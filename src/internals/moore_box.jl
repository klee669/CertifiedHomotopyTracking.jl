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
    state = (y=copy(x), A=copy(A_init), t=t)
    return _refine_moore_box(sys, state, sys.RR(r_init); rho=tau)
end

function _refine_moore_box(
    problem,
    state,
    r;
    rho,
    min_radius=0,
    max_radius=1.0,
    growth_factor=2.0,
    radius_shrink_factor=0.5,
    max_iter=20,
    kwargs...,
)
    growth_factor > 1 || throw(ArgumentError("growth_factor must be greater than 1."))
    0 < radius_shrink_factor < 1 || throw(ArgumentError("radius_shrink_factor must lie between 0 and 1."))
    last_k_norm = Inf

    for _ in 1:max_iter
        state = _moore_prepare(problem, state; kwargs...)
        passed, k_norm = _moore_test(problem, state, r, rho)
        last_k_norm = k_norm
        if passed
            while growth_factor * r <= max_radius
                r_next = growth_factor * r
                grown, grown_norm = _moore_test(problem, state, r_next, rho)
                grown || break
                r = r_next
                last_k_norm = grown_norm
            end
            return _moore_success(problem, state, r, last_k_norm)
        end

        correction = _moore_correction(problem, state)
        if _moore_should_shrink(problem, state, correction, r, rho)
            r *= radius_shrink_factor
            r < min_radius && break
        else
            state = _moore_apply_correction(problem, state, correction; kwargs...)
        end
    end

    return _moore_failure(problem, state, r, last_k_norm)
end

_moore_prepare(::HCSystem, state; kwargs...) = state

_moore_test(sys::HCSystem, state, r, rho) =
    krawczyk_test(sys, state.y, state.t, r, state.A; rho=rho)

_moore_correction(sys::HCSystem, state) = state.A * evaluate_H(sys, state.y, state.t)

_moore_should_shrink(::HCSystem, state, delta, r, rho) =
    norm_inf(delta) <= (1/64) * rho * Float64(r)

function _moore_apply_correction(sys::HCSystem, state, delta; kwargs...)
    y = get_mid_vec(state.y - delta)
    A = inv_acb(evaluate_Jac(sys, y, state.t))
    return (y=y, A=A, t=state.t)
end

_moore_success(::HCSystem, state, r, k_norm) = (state.y, r, state.A, true)

_moore_failure(::HCSystem, state, r, k_norm) = (state.y, r, state.A, false)


function refine_moore_box(
    variety::AlgebraicVarietySystem,
    x,
    r_init;
    rho=1/8,
    rank_tol=1e-10,
    min_radius=1e-12,
    max_radius=1.0,
    growth_factor=2.0,
    radius_shrink_factor=0.5,
    max_iter=20,
    newton_tol=1e-24,
    max_newton_steps=10,
)
    anchor = get_mid_vec(_as_acb_vector(variety.system.CC, x))
    frame = local_tangent_normal_frame(variety, anchor; rank_tol=rank_tol)
    state = (
        anchor=anchor,
        frame=frame,
    )
    return _refine_moore_box(
        variety,
        state,
        Float64(r_init);
        rho=rho,
        min_radius=min_radius,
        max_radius=max_radius,
        growth_factor=growth_factor,
        radius_shrink_factor=radius_shrink_factor,
        max_iter=max_iter,
        rank_tol=rank_tol,
        newton_tol=newton_tol,
        max_newton_steps=max_newton_steps,
    )
end

function _moore_prepare(
    variety::AlgebraicVarietySystem,
    state;
    rank_tol,
    newton_tol,
    max_newton_steps,
    kwargs...,
)
    anchor, frame, _ = _polish_variety_center(
        variety,
        state.anchor;
        rank_tol=rank_tol,
        newton_tol=newton_tol,
        max_newton_steps=max_newton_steps,
    )
    return (
        anchor=anchor,
        frame=frame,
    )
end

_moore_test(variety::AlgebraicVarietySystem, state, r, rho) =
    krawczyk_test(variety, state.anchor, state.frame, 0.0, r; rho=rho)

_moore_correction(::AlgebraicVarietySystem, state) = nothing

_moore_should_shrink(::AlgebraicVarietySystem, state, correction, r, rho) = true

_moore_apply_correction(::AlgebraicVarietySystem, state, correction; kwargs...) = state

function _moore_success(::AlgebraicVarietySystem, state, r, k_norm)
    return VarietyBox(state.anchor, state.frame, 0.0, r, k_norm, true)
end

function _moore_failure(::AlgebraicVarietySystem, state, r, k_norm)
    return VarietyBox(state.anchor, state.frame, 0.0, r, k_norm, false)
end
