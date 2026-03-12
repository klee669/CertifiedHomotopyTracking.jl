export jac,
       evaluate_matrix,
       jacobian_inverse,
       interval_svd,
       inv_acb,
       numerical_kernel,
       create_seed_pair


# --------------------------------------------------------------------------
# Compute the symbolic Jacobian of a system (matrix form)
function jac(system::Union{Matrix, AbstractAlgebra.Generic.MatSpaceElem})
    R    = parent(system[1])
    vars = gens(R)
    n    = length(system)

    J = zero_matrix(R, n, n)
    for j in 1:n
        for i in 1:n
            d = derivative(system[j], vars[i])
            J[j, i] = d == 0 ? zero(R) : d
        end
    end
    return Matrix(J)  # Always return as Julia Matrix
end


# --------------------------------------------------------------------------
# Evaluate a symbolic matrix at a given point (number or vector)
function evaluate_matrix(
    m::Union{Matrix, AbstractAlgebra.Generic.MatSpaceElem},
    vec::Union{Number, AcbFieldElem, Vector{AcbFieldElem}, AbstractAlgebra.Generic.MPoly}
)
    R = base_ring(m[1])
    nr, nc = size(m)
    out = zero_matrix(R, nr, nc)

    for j in 1:nr
        for i in 1:nc
            out[j, i] = AbstractAlgebra.evaluate(m[j, i], vec)
        end
    end

    return out  # Avoid returning AcbMatrix
end


# --------------------------------------------------------------------------
# Inverse of evaluated Jacobian with fallback
function jacobian_inverse(
    system::Union{Matrix, AbstractAlgebra.Generic.MatSpaceElem},
    vec::Vector{AcbFieldElem}
)
    J  = jac(system)
    Jv = evaluate_matrix(J, vec)

    Jinv, factor = pseudo_inv(Jv)
    if radius(real(prod(Jinv))) > 1e300
        return inv(Jv)
    end

    return (1 / factor) * Jinv
end


# --------------------------------------------------------------------------
# Interval-based SVD decomposition
function interval_svd(M::Matrix)
    R = parent(M[1, 1])
    M_real = convert_to_double_matrix(M)
    u, s, vt = svd(M_real; full = true)

    return [
        convert_to_box_matrix(u, R),
        convert_to_box_matrix(hcat(s), R),
        convert_to_box_matrix(vt, R)
    ]
end


# --------------------------------------------------------------------------
# Compute numerical kernel based on interval SVD
function numerical_kernel(M::Matrix; tol = 1e-5)
    R = parent(M[1, 1])
    u, s, vt = interval_svd(M)
    m, n = size(M)

    if m == 0
        return identity_matrix(R, 3)
    elseif n == 0
        return zeros(R, 0, 0)
    end

    cols = findall(i -> max_int_norm(i) < tol, vec(s))
    if isempty(cols)
        return zeros(R, m, 0)
    end

    return matrix(vt[:, cols])  # safe return of kernel vectors
end


# --------------------------------------------------------------------------
# Generate a seed pair (initial condition + perturbation)
function create_seed_pair(
    F::Matrix,
    p::AbstractVector,
)
    CC  = parent(p[1])
    R   = parent(F[1])
    n_p = length(gens(R))
    n   = length(F)

    id_mat = identity_matrix(CC, n_p)

    # b: evaluation of F at p (zero-η-t form)
    b = evaluate_matrix(evaluate_matrix(evaluate_matrix(hcat(F), zeros(CC, n_p)), 0), p)

    # Construct numerical Jacobian A (approximate)
    A = zeros(CC, n, n_p)
    for i in 1:n_p
        shift = id_mat[i, :]
        A[:, i] = Matrix(transpose(
            evaluate_matrix(evaluate_matrix(evaluate_matrix(hcat(F), shift), 0), p) - b
        ))
    end

    # Solve normal equation to get xp
    cA = convert_to_double_matrix(A)
    AtA_box = convert_to_box_matrix(cA' * cA, CC)
    pseudo = pseudo_inv(matrix(AtA_box))

    xp = Matrix((1 / pseudo[2]) * pseudo[1]) *
         convert_to_box_matrix(cA', CC) *
         vec(Matrix(-b))

    # Compute xh in the nullspace of A (perturbation)
    K = Matrix(numerical_kernel(A))
    xh = size(K, 2) == 0 ? zeros(CC, n, 1) :
         K * convert_to_box_matrix(rand(size(K, 2), 1), CC)

    return (Matrix(transpose(xp + xh)), p)
end


function inv_acb(A::AbstractMatrix)
    if eltype(A) !== AcbFieldElem
        n_rows, n_cols = size(A)
        A_typed = Matrix{AcbFieldElem}(undef, n_rows, n_cols)
        for i in 1:n_rows, j in 1:n_cols
            A_typed[i,j] = CC(A[i,j])
        end
        A = A_typed
    end

    R = CC 
    n = size(A, 1)
    
    Mat = matrix(R, A) 
    InvMat = inv(Mat)
    
    res = Matrix{AcbFieldElem}(undef, n, n)
    for i in 1:n, j in 1:n
        res[i,j] = InvMat[i,j]
    end
    return res
end
