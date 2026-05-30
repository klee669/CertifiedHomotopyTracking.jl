export krawczyk_operator, krawczyk_test, compute_preconditioner, validate_step_taylor3

mutable struct KrawczykValidationCache
    tm_cache::TMCache
    B::Vector{AcbFieldElem}
    I_mat::Matrix{AcbFieldElem}
    F_val::Vector{AcbFieldElem}
    X_bound::Vector{AcbFieldElem}
    X_expanded::Vector{AcbFieldElem}
end

function _fill_identity!(I_mat::AbstractMatrix{AcbFieldElem}, CC::AcbField)
    n, m = size(I_mat)
    n == m || throw(DimensionMismatch("identity matrix cache must be square."))
    for j in 1:n, i in 1:n
        I_mat[i,j] = (i == j ? CC(1) : CC(0))
    end
    return I_mat
end

function _fill_unit_box!(B::AbstractVector{AcbFieldElem}, CC::AcbField, RR::ArbField)
    one_int = RR("0 +/- 1")
    b_int = CC(one_int, one_int)
    fill!(B, b_int)
    return B
end

function KrawczykValidationCache(CC::AcbField, RR::ArbField, n::Integer)
    B = Vector{AcbFieldElem}(undef, n)
    _fill_unit_box!(B, CC, RR)
    I_mat = Matrix{AcbFieldElem}(undef, n, n)
    _fill_identity!(I_mat, CC)
    return KrawczykValidationCache(
        TMCache(CC),
        B,
        I_mat,
        [CC(0) for _ in 1:n],
        [CC(0) for _ in 1:n],
        [CC(0) for _ in 1:n],
    )
end

function _validation_cache(sys::HCSystem, n::Integer, cache)
    if cache === nothing
        return KrawczykValidationCache(sys.CC, sys.RR, n)
    end
    length(cache.B) == n || throw(DimensionMismatch("validation cache dimension does not match the system."))
    return cache
end


####### functions for Krawczyk test
function krawczyk_operator(
    system::Union{Matrix,AbstractAlgebra.Generic.MatSpaceElem}, 
    point::Vector{AcbFieldElem}, 
    r::Number, 
    A::AcbMatrix
)

    n   = length(system)
    CC  = base_ring(system[1])
    B   = CC("+/- 1", "+/-1")  # Unit interval ball
    I   = identity_matrix(CC, n)

    # Build unit box matrix
    box = fill(B, n, 1)

    # Evaluate system and its Jacobian
    eval_sys = evaluate_matrix(system, point)
    j        = jac(system)
    eval_jac = evaluate_matrix(j, vec(point + r * box))

    # Krawczyk operator
    K = (-1 / r) * A * transpose(eval_sys) + (I - A * eval_jac) * matrix(box)
    return K
end


# Krawczyk test function.
# checking if the norm of the Krawczyk operator is smaller than rho.
function krawczyk_test(
    system::Union{Matrix,AbstractAlgebra.Generic.MatSpaceElem}, 
    point::Vector{AcbFieldElem}, 
    r::Number, 
    A::AcbMatrix,
    ρ::Number
)
    K = krawczyk_operator(system,point,r, A);
    max_norm(K) < ρ
end


function compute_preconditioner(sys::HCSystem, x::AbstractVector{AcbFieldElem}, t)
    CC = sys.CC
    x_mid = get_mid_vec(x)
    t_mid = t isa AcbFieldElem ? get_mid(t) : CC(t)
    J_val = evaluate_Jac(sys, x_mid, t_mid) 
    return inv_acb(J_val)
end

function krawczyk_test(sys::HCSystem, x::AbstractVector{AcbFieldElem}, t, r; rho=0.7)
    A = compute_preconditioner(sys, x, t)
    return krawczyk_test(sys, x, t, r, A; rho=rho)
end

function krawczyk_test(sys::HCSystem, x::AbstractVector{AcbFieldElem}, t, r, A::AbstractMatrix{AcbFieldElem}; rho=0.7)
    CC = sys.CC; RR = sys.RR
    n = length(x)

    one_int = RR("0 +/- 1") 
    b_int = CC(one_int, one_int)
    B = [b_int for _ in 1:n]
    
    fx = evaluate_H(sys, x, t)
    x_expanded = x .+ (B .* CC(r))
    Jx = evaluate_Jac(sys, x_expanded, t)
    
    term1 = -(A * fx) ./ CC(r)
    
    I_mat = Matrix{AcbFieldElem}(undef, n, n)
    for i in 1:n, j in 1:n
        I_mat[i,j] = (i == j ? CC(1) : CC(0))
    end
    
    term2 = (I_mat - A * Jx) * B
    K = term1 + term2
    k_norm = norm_inf(K)
    
    return k_norm < rho, k_norm
end

function validate_step_taylor3(sys::HCSystem, X_tm::Vector{<:TaylorModel3}, t_start, h, r, A; rho=0.7, cache=nothing)
    CC = sys.CC; RR = sys.RR
    n = length(X_tm)
    validation_cache = _validation_cache(sys, n, cache)
    t_tm = TaylorModel3(CC(t_start), CC(1), CC(0), CC(0), CC(0), RR(h))
    
    F_tm = evaluate_H(sys, X_tm, t_tm)
    
    for i in 1:n
        evaluate_taylor!(validation_cache.F_val[i], F_tm[i], validation_cache.tm_cache)
        evaluate_taylor!(validation_cache.X_bound[i], X_tm[i], validation_cache.tm_cache)
    end
    
    r_cc = CC(r)
    for i in 1:n
        validation_cache.X_expanded[i] = validation_cache.X_bound[i] + validation_cache.B[i] * r_cc
    end
    T_expanded = CC(t_start) + CC(RR(0), RR(h))
    
    J_val = evaluate_Jac(sys, validation_cache.X_expanded, T_expanded)
    
    term1 = -(A * validation_cache.F_val) ./ r_cc
    
    term2 = (validation_cache.I_mat - A * J_val) * validation_cache.B
    
    K = term1 + term2
    norm_K = norm_inf(K)
    
    return norm_K < rho, norm_K
end
