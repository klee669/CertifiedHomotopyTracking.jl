export krawczyk_operator, krawczyk_test, compute_preconditioner, validate_step_taylor3,
       validate_step_taylor3_diagnostics

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

function _acb_identity_matrix(CC::AcbField, n::Integer)
    I_mat = Matrix{AcbFieldElem}(undef, n, n)
    return _fill_identity!(I_mat, CC)
end

function _acb_unit_box_vector(CC::AcbField, RR::ArbField, n::Integer)
    B = Vector{AcbFieldElem}(undef, n)
    return _fill_unit_box!(B, CC, RR)
end

function KrawczykValidationCache(CC::AcbField, RR::ArbField, n::Integer)
    B = _acb_unit_box_vector(CC, RR, n)
    I_mat = _acb_identity_matrix(CC, n)
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

    B = _acb_unit_box_vector(CC, RR, n)
    
    fx = evaluate_H(sys, x, t)
    x_expanded = x .+ (B .* CC(r))
    Jx = evaluate_Jac(sys, x_expanded, t)
    
    term1 = -(A * fx) ./ CC(r)
    
    I_mat = _acb_identity_matrix(CC, n)
    
    term2 = (I_mat - A * Jx) * B
    K = term1 + term2
    k_norm = norm_inf(K)
    
    return k_norm < rho, k_norm
end

function _matrix_inf_norm_bound(M::AbstractMatrix{AcbFieldElem})
    row_max = 0.0
    for i in axes(M, 1)
        row_sum = 0.0
        for j in axes(M, 2)
            row_sum += mag_complex(M[i, j])
        end
        row_max = max(row_max, row_sum)
    end
    return row_max
end

function _profiled_validation_step(f, profile_validation::Bool)
    if profile_validation
        timed = @timed f()
        return timed.value, timed.time, timed.bytes
    else
        return f(), 0.0, 0
    end
end

function validate_step_taylor3_diagnostics(
    sys::HCSystem,
    X_tm::Vector{<:TaylorModel3},
    t_start,
    h,
    r,
    A;
    rho=0.7,
    cache=nothing,
    profile_validation=false,
)
    CC = sys.CC; RR = sys.RR
    n = length(X_tm)
    validation_cache = _validation_cache(sys, n, cache)
    t_tm, time_build_t, alloc_build_t = _profiled_validation_step(profile_validation) do
        TaylorModel3(CC(t_start), CC(1), CC(0), CC(0), CC(0), RR(h))
    end

    F_tm, time_evaluate_H, alloc_evaluate_H = _profiled_validation_step(profile_validation) do
        evaluate_H(sys, X_tm, t_tm)
    end

    _, time_evaluate_taylor, alloc_evaluate_taylor = _profiled_validation_step(profile_validation) do
        for i in 1:n
            evaluate_taylor!(validation_cache.F_val[i], F_tm[i], validation_cache.tm_cache)
            evaluate_taylor!(validation_cache.X_bound[i], X_tm[i], validation_cache.tm_cache)
        end
        nothing
    end

    r_cc = CC(r)
    _, time_expand_X, alloc_expand_X = _profiled_validation_step(profile_validation) do
        for i in 1:n
            validation_cache.X_expanded[i] = validation_cache.X_bound[i] + validation_cache.B[i] * r_cc
        end
        nothing
    end

    T_expanded, time_build_T, alloc_build_T = _profiled_validation_step(profile_validation) do
        h_rr = RR(h)
        CC(t_start) + _tm_real_interval(CC, h_rr)
    end

    J_val, time_evaluate_Jac, alloc_evaluate_Jac = _profiled_validation_step(profile_validation) do
        evaluate_Jac(sys, validation_cache.X_expanded, T_expanded)
    end

    AH, time_AH, alloc_AH = _profiled_validation_step(profile_validation) do
        A * validation_cache.F_val
    end

    term1, time_term1, alloc_term1 = _profiled_validation_step(profile_validation) do
        -AH ./ r_cc
    end

    linear_defect, time_AJ, alloc_AJ = _profiled_validation_step(profile_validation) do
        validation_cache.I_mat - A * J_val
    end

    term2, time_term2, alloc_term2 = _profiled_validation_step(profile_validation) do
        linear_defect * validation_cache.B
    end

    K, time_K, alloc_K = _profiled_validation_step(profile_validation) do
        term1 + term2
    end

    norms, time_norms, alloc_norms = _profiled_validation_step(profile_validation) do
        Y = norm_inf(AH)
        Z = _matrix_inf_norm_bound(linear_defect)
        Y_over_r = Y / Float64(r)
        (norm_K = norm_inf(K), Y = Y, Z = Z, Y_over_r = Y_over_r, yz_bound = Y_over_r + Z)
    end
    
    return (
        passed = norms.norm_K < rho,
        norm_K = norms.norm_K,
        Y = norms.Y,
        Z = norms.Z,
        Y_over_r = norms.Y_over_r,
        yz_bound = norms.yz_bound,
        radius = Float64(r),
        profile = (
            time_build_t = time_build_t,
            time_evaluate_H = time_evaluate_H,
            time_evaluate_taylor = time_evaluate_taylor,
            time_expand_X = time_expand_X,
            time_build_T = time_build_T,
            time_evaluate_Jac = time_evaluate_Jac,
            time_AH = time_AH,
            time_term1 = time_term1,
            time_AJ = time_AJ,
            time_term2 = time_term2,
            time_K = time_K,
            time_norms = time_norms,
            alloc_build_t = alloc_build_t,
            alloc_evaluate_H = alloc_evaluate_H,
            alloc_evaluate_taylor = alloc_evaluate_taylor,
            alloc_expand_X = alloc_expand_X,
            alloc_build_T = alloc_build_T,
            alloc_evaluate_Jac = alloc_evaluate_Jac,
            alloc_AH = alloc_AH,
            alloc_term1 = alloc_term1,
            alloc_AJ = alloc_AJ,
            alloc_term2 = alloc_term2,
            alloc_K = alloc_K,
            alloc_norms = alloc_norms,
        ),
    )
end

function validate_step_taylor3(sys::HCSystem, X_tm::Vector{<:TaylorModel3}, t_start, h, r, A; rho=0.7, cache=nothing)
    diagnostic = validate_step_taylor3_diagnostics(
        sys,
        X_tm,
        t_start,
        h,
        r,
        A;
        rho = rho,
        cache = cache,
        profile_validation = false,
    )
    return diagnostic.passed, diagnostic.norm_K
end
