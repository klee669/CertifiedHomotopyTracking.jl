export krawczyk_operator, krawczyk_test, compute_preconditioner

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

"""
    KrawczykWorkspace

Task-local workspace for [`krawczyk_operator`](@ref), keyed by precision and
dimension.

`B` (the unit interval box) and `I_mat` (the identity) depend only on those two
keys and are read-only operands, so they are built once instead of rebuilding
`n + n^2` ball elements per step. The remaining fields are scratch: `x_expanded`
holds the inflated evaluation point, and `acc`/`tmp`/`t1`/`t2`/`mij` accumulate
the operator row by row so that no intermediate matrix is materialized.
"""
struct KrawczykWorkspace
    prec::Int
    n::Int
    B::Vector{AcbFieldElem}
    I_mat::Matrix{AcbFieldElem}
    x_expanded::Vector{AcbFieldElem}
    acc::AcbFieldElem
    tmp::AcbFieldElem
    t1::AcbFieldElem
    t2::AcbFieldElem
    mij::AcbFieldElem
    kk::AcbFieldElem
end

_acb_entry(::AcbField, v::AcbFieldElem) = v
_acb_entry(CC::AcbField, v) = CC(v)

_as_acb_matrix(::AcbField, M::AbstractMatrix{AcbFieldElem}, ::Integer) = M
_as_acb_matrix(CC::AcbField, M::AbstractMatrix, n::Integer) =
    AcbFieldElem[_acb_entry(CC, M[i, j]) for i in 1:n, j in 1:n]

function _krawczyk_workspace(CC::AcbField, RR::ArbField, n::Integer)
    store = task_local_storage()
    cached = get(store, :cht_krawczyk_workspace, nothing)
    if cached isa KrawczykWorkspace && cached.prec == precision(CC) && cached.n == n
        return cached
    end
    ws = KrawczykWorkspace(
        precision(CC),
        Int(n),
        _acb_unit_box_vector(CC, RR, n),
        _acb_identity_matrix(CC, n),
        [CC(0) for _ in 1:n],
        CC(), CC(), CC(), CC(), CC(), CC(),
    )
    store[:cht_krawczyk_workspace] = ws
    return ws
end

# Fused K = -(A*F)/r + (I - A*J)*B, accumulating the infinity norms of K, A*F
# and I - A*J in the same pass so that no intermediate vector or matrix is
# materialized. Returns (norm_K, Y, Z); the operation order matches the
# unfused expressions used when `profile_validation=true`.
function _fused_validation_norms(A, F_val, J_val, I_mat, B, r_cc, ws::KrawczykWorkspace)
    n = length(F_val)
    acc = ws.acc; tmp = ws.tmp; t1 = ws.t1; t2 = ws.t2; mij = ws.mij; kk = ws.kk
    norm_K = 0.0; Y = 0.0; Z = 0.0
    for i in 1:n
        Nemo.mul!(acc, A[i, 1], F_val[1])
        for k in 2:n
            Nemo.mul!(tmp, A[i, k], F_val[k])
            Nemo.add!(acc, acc, tmp)
        end
        Y = max(Y, mag_complex(acc))
        Nemo.neg!(t1, acc)
        Nemo.div!(t1, t1, r_cc)

        row_sum = 0.0
        for j in 1:n
            Nemo.mul!(acc, A[i, 1], J_val[1, j])
            for k in 2:n
                Nemo.mul!(tmp, A[i, k], J_val[k, j])
                Nemo.add!(acc, acc, tmp)
            end
            Nemo.sub!(mij, I_mat[i, j], acc)
            row_sum += mag_complex(mij)
            if j == 1
                Nemo.mul!(t2, mij, B[1])
            else
                Nemo.mul!(tmp, mij, B[j])
                Nemo.add!(t2, t2, tmp)
            end
        end
        Z = max(Z, row_sum)

        Nemo.add!(kk, t1, t2)
        norm_K = max(norm_K, mag_complex(kk))
    end
    return norm_K, Y, Z
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

function _validation_cache(sys::SpecializedHomotopy, n::Integer, cache)
    if cache === nothing
        return KrawczykValidationCache(sys.CC, sys.RR, n)
    end
    length(cache.B) == n || throw(DimensionMismatch("validation cache dimension does not match the system."))
    return cache
end


####### functions for Krawczyk test
"""
    krawczyk_operator(system, point, r, A)
    krawczyk_operator(sys::SpecializedHomotopy, x, t, r)
    krawczyk_operator(sys::SpecializedHomotopy, x, t, r, A)

Compute the legacy polynomial-ring Krawczyk operator for `system` at `point`
with radius `r` and preconditioner `A`.

For current [`SpecializedHomotopy`](@ref)-based code, the operator is the
interval vector whose infinity norm is tested by [`krawczyk_test`](@ref).
"""
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


"""
    krawczyk_test(sys::SpecializedHomotopy, x, t, r; rho=0.7)
    krawczyk_test(sys::SpecializedHomotopy, x, t, r, A; rho=0.7)

Run a Krawczyk inclusion/contraction test.
For a given system and point, returns `(passed, k_norm)`. 
"""
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


"""
    compute_preconditioner(sys::SpecializedHomotopy, x, t)

Evaluate the Jacobian at the midpoint of `(x, t)` and return its inverse as the
preconditioner used by Krawczyk tests.
"""
function compute_preconditioner(sys::SpecializedHomotopy, x::AbstractVector{AcbFieldElem}, t)
    CC = sys.CC
    x_mid = get_mid_vec(x)
    t_mid = t isa AcbFieldElem ? get_mid(t) : CC(t)
    J_val = evaluate_Jac(sys, x_mid, t_mid) 
    return inv_acb(J_val, CC)
end

function krawczyk_test(sys::SpecializedHomotopy, x::AbstractVector{AcbFieldElem}, t, r; rho=0.7)
    A = compute_preconditioner(sys, x, t)
    return krawczyk_test(sys, x, t, r, A; rho=rho)
end

function krawczyk_operator(sys::SpecializedHomotopy, x::AbstractVector{AcbFieldElem}, t, r)
    A = compute_preconditioner(sys, x, t)
    return krawczyk_operator(sys, x, t, r, A)
end

krawczyk_operator(sys::SpecializedHomotopy, x::AbstractVector{AcbFieldElem}, t, r, A::AbstractMatrix{AcbFieldElem}) =
    krawczyk_operator(sys, x, t, r, A, evaluate_H(sys, x, t))

# `fx` is H(x, t). It does not depend on the radius, so callers sweeping r over
# a range (the Moore box growth loop) evaluate it once and pass it in.
function krawczyk_operator(sys::SpecializedHomotopy, x::AbstractVector{AcbFieldElem}, t, r, A::AbstractMatrix{AcbFieldElem}, fx::AbstractVector)
    CC = sys.CC; RR = sys.RR
    n = length(x)

    ws = _krawczyk_workspace(CC, RR, n)
    B = ws.B
    I_mat = ws.I_mat
    r_cc = CC(r)

    x_expanded = ws.x_expanded
    for i in 1:n
        Nemo.mul!(x_expanded[i], B[i], r_cc)
        Nemo.add!(x_expanded[i], x[i], x_expanded[i])
    end
    # Constant Jacobian entries come back unboxed (e.g. Int), which makes the
    # generated matrix `Matrix{Any}` and every entry a dynamic dispatch below.
    # Converting a constant to a ball at this precision is exact.
    Jx = _as_acb_matrix(CC, evaluate_Jac(sys, x_expanded, t), n)

    # K = -(A*fx)/r + (I - A*Jx)*B, accumulated one row at a time so that the
    # intermediate A*Jx and I - A*Jx matrices are never materialized. The
    # operation order matches the unfused expression exactly.
    acc = ws.acc; tmp = ws.tmp; t1 = ws.t1; t2 = ws.t2; mij = ws.mij
    K = Vector{AcbFieldElem}(undef, n)
    for i in 1:n
        Nemo.mul!(acc, A[i, 1], fx[1])
        for k in 2:n
            Nemo.mul!(tmp, A[i, k], fx[k])
            Nemo.add!(acc, acc, tmp)
        end
        Nemo.neg!(t1, acc)
        Nemo.div!(t1, t1, r_cc)

        for j in 1:n
            Nemo.mul!(acc, A[i, 1], Jx[1, j])
            for k in 2:n
                Nemo.mul!(tmp, A[i, k], Jx[k, j])
                Nemo.add!(acc, acc, tmp)
            end
            Nemo.sub!(mij, I_mat[i, j], acc)
            if j == 1
                Nemo.mul!(t2, mij, B[1])
            else
                Nemo.mul!(tmp, mij, B[j])
                Nemo.add!(t2, t2, tmp)
            end
        end

        K[i] = t1 + t2
    end
    return K
end

krawczyk_test(sys::SpecializedHomotopy, x::AbstractVector{AcbFieldElem}, t, r, A::AbstractMatrix{AcbFieldElem}; rho=0.7) =
    krawczyk_test(sys, x, t, r, A, evaluate_H(sys, x, t); rho=rho)

function krawczyk_test(sys::SpecializedHomotopy, x::AbstractVector{AcbFieldElem}, t, r, A::AbstractMatrix{AcbFieldElem}, fx::AbstractVector; rho=0.7)
    K = krawczyk_operator(sys, x, t, r, A, fx)
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

_empty_validation_profile() = (
    time_build_t = 0.0,
    time_evaluate_H = 0.0,
    time_evaluate_taylor = 0.0,
    time_expand_X = 0.0,
    time_build_T = 0.0,
    time_evaluate_Jac = 0.0,
    time_AH = 0.0,
    time_term1 = 0.0,
    time_AJ = 0.0,
    time_term2 = 0.0,
    time_K = 0.0,
    time_norms = 0.0,
    alloc_build_t = 0,
    alloc_evaluate_H = 0,
    alloc_evaluate_taylor = 0,
    alloc_expand_X = 0,
    alloc_build_T = 0,
    alloc_evaluate_Jac = 0,
    alloc_AH = 0,
    alloc_term1 = 0,
    alloc_AJ = 0,
    alloc_term2 = 0,
    alloc_K = 0,
    alloc_norms = 0,
)

function _profiled_validation_step(f, profile_validation::Bool)
    if profile_validation
        timed = @timed f()
        return timed.value, timed.time, timed.bytes
    else
        return f(), 0.0, 0
    end
end

"""
    validate_step_taylor3_diagnostics(sys, X_tm, t_start, h, r, A; kwargs...)

Validate one Taylor-model predictor step and return diagnostic data.

Options include:

- `rho=0.7`: Krawczyk threshold.
- `cache=nothing`: optional `KrawczykValidationCache`-compatible cache.
- `profile_validation=false`: include timing information when enabled.
"""
function validate_step_taylor3_diagnostics(
    sys::SpecializedHomotopy,
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

    if !profile_validation
        # Hot path: one fused pass, no intermediate vectors or matrices.
        ws = _krawczyk_workspace(CC, RR, n)
        norm_K, Y, Z = _fused_validation_norms(
            A,
            validation_cache.F_val,
            _as_acb_matrix(CC, J_val, n),
            validation_cache.I_mat,
            validation_cache.B,
            r_cc,
            ws,
        )
        Y_over_r = Y / Float64(r)
        return (
            passed = norm_K < rho,
            norm_K = norm_K,
            Y = Y,
            Z = Z,
            Y_over_r = Y_over_r,
            yz_bound = Y_over_r + Z,
            radius = Float64(r),
            profile = _empty_validation_profile(),
        )
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

"""
    validate_step_taylor3(sys, X_tm, t_start, h, r, A; rho=0.7, cache=nothing)

Return `(passed, k_norm)` for one Taylor-model path-tracking step.
"""
function validate_step_taylor3(sys::SpecializedHomotopy, X_tm::Vector{<:TaylorModel3}, t_start, h, r, A; rho=0.7, cache=nothing)
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
