export homogenize_poly,
    homogenize_expr,
    homogenize_system,
    random_patch_vector,
    patch_equation,
    lift_to_patch,
    repatch

const _HOMOGENIZE_ZERO_NUM = Num(0)

_homogenize_as_num(x) = x isa Num ? x : Num(x)

function _homogenize_symbolic_const_value(x)
    x isa Number && !(x isa Num) && return x
    isempty(Symbolics.get_variables(x)) || return nothing
    u = Symbolics.unwrap(x)
    return Symbolics.SymbolicUtils.isconst(u) ? Symbolics.SymbolicUtils.unwrap_const(u) : nothing
end

function _homogenize_is_symbolic_zero(expr)
    return isequal(Symbolics.simplify(expr), _HOMOGENIZE_ZERO_NUM)
end

function _rational_numerator(expr::Complex)
    real_num, real_den = _rational_num_den(Symbolics.simplify_fractions(real(expr)))
    imag_num, imag_den = _rational_num_den(Symbolics.simplify_fractions(imag(expr)))
    return Symbolics.expand(real_num * imag_den + im * imag_num * real_den)
end

function _rational_numerator(expr)
    num, _ = _rational_num_den(Symbolics.simplify_fractions(expr))
    return Symbolics.expand(num)
end

function _rational_num_den(expr)
    expr_num = _homogenize_as_num(expr)
    unwrapped = Symbolics.unwrap(expr_num)
    op = try
        Symbolics.operation(unwrapped)
    catch
        return expr_num, Num(1)
    end

    args = Symbolics.arguments(unwrapped)

    if op === (+)
        num = Num(0)
        den = Num(1)
        for arg in args
            arg_num, arg_den = _rational_num_den(_homogenize_as_num(arg))
            num = Symbolics.expand(num * arg_den + arg_num * den)
            den = Symbolics.expand(den * arg_den)
        end
        return num, den
    elseif op === (-)
        if length(args) == 1
            arg_num, arg_den = _rational_num_den(_homogenize_as_num(args[1]))
            return Symbolics.expand(-arg_num), arg_den
        elseif length(args) == 2
            left_num, left_den = _rational_num_den(_homogenize_as_num(args[1]))
            right_num, right_den = _rational_num_den(_homogenize_as_num(args[2]))
            return Symbolics.expand(left_num * right_den - right_num * left_den), Symbolics.expand(left_den * right_den)
        end
    elseif op === (*)
        num = Num(1)
        den = Num(1)
        for arg in args
            arg_num, arg_den = _rational_num_den(_homogenize_as_num(arg))
            num = Symbolics.expand(num * arg_num)
            den = Symbolics.expand(den * arg_den)
        end
        return num, den
    elseif op === (/)
        numerator_num, numerator_den = _rational_num_den(_homogenize_as_num(args[1]))
        denominator_num, denominator_den = _rational_num_den(_homogenize_as_num(args[2]))
        return Symbolics.expand(numerator_num * denominator_den), Symbolics.expand(numerator_den * denominator_num)
    elseif op === (^)
        exponent = _homogenize_symbolic_const_value(_homogenize_as_num(args[2]))
        if !(exponent isa Integer)
            throw(ArgumentError("Expected a polynomial or rational expression with integer exponents; found exponent $(args[2])."))
        end
        base_num, base_den = _rational_num_den(_homogenize_as_num(args[1]))
        if exponent >= 0
            return Symbolics.expand(base_num^exponent), Symbolics.expand(base_den^exponent)
        else
            return Symbolics.expand(base_den^(-exponent)), Symbolics.expand(base_num^(-exponent))
        end
    end

    return expr_num, Num(1)
end

function _assert_polynomial_and_total_degree(expr::Complex, vars)
    real_degree = _assert_polynomial_and_total_degree(real(expr), vars)
    imag_degree = _assert_polynomial_and_total_degree(imag(expr), vars)
    return max(real_degree, imag_degree)
end

function _assert_polynomial_and_total_degree(expr, vars)
    coeffs, remainder = Symbolics.polynomial_coeffs(Symbolics.expand(expr), vars)
    if !_homogenize_is_symbolic_zero(remainder)
        throw(ArgumentError("Expected a polynomial expression in the provided variables; non-polynomial remainder: $(remainder)."))
    end

    isempty(coeffs) && return 0

    max_degree = 0
    for monomial in keys(coeffs)
        degree = sum(Symbolics.degree(monomial, var) for var in vars)
        if !(degree isa Integer) || degree < 0
            throw(ArgumentError("Expected a polynomial expression with nonnegative integer exponents; found monomial $(monomial) with degree $(degree)."))
        end
        max_degree = max(max_degree, degree)
    end
    return max_degree
end

function _homogenize_with_degree(expr, vars, h_var, degree::Integer)
    substitutions = Dict(var => var / h_var for var in vars)
    expr_frac = Symbolics.substitute(expr, substitutions)
    return Symbolics.expand(Symbolics.simplify(expr_frac * h_var^degree))
end

function _homogenize_with_degree(expr::Complex, vars, h_var, degree::Integer)
    real_part = _homogenize_with_degree(real(expr), vars, h_var, degree)
    imag_part = _homogenize_with_degree(imag(expr), vars, h_var, degree)
    return Symbolics.expand(real_part + im * imag_part)
end

function homogenize_poly(expr, vars, u0_var)
    numerator = _rational_numerator(expr)
    degree = _assert_polynomial_and_total_degree(numerator, vars)
    return _homogenize_with_degree(numerator, vars, u0_var, degree)
end

function homogenize_expr(expr, vars, h_var)
    numerator = _rational_numerator(expr)
    degree = _assert_polynomial_and_total_degree(numerator, vars)
    return _homogenize_with_degree(numerator, vars, h_var, degree)
end

function homogenize_system(exprs, vars, h_var)
    return [homogenize_expr(expr, vars, h_var) for expr in exprs]
end

function random_patch_vector(n::Integer; rng=nothing, normalize=true)
    n > 0 || throw(ArgumentError("Patch vector length must be positive."))
    real_part = rng === nothing ? randn(n) : randn(rng, n)
    imag_part = rng === nothing ? randn(n) : randn(rng, n)
    a = complex.(real_part, imag_part)
    if normalize
        a_norm = norm(a)
        a_norm > 0 || throw(ArgumentError("Cannot normalize a zero patch vector."))
        return a ./ a_norm
    end
    return a
end

function patch_equation(projective_vars, a)
    length(projective_vars) == length(a) || throw(DimensionMismatch("projective_vars and a must have the same length."))
    !isempty(projective_vars) || throw(ArgumentError("projective_vars must be nonempty."))

    result = zero(first(a) * first(projective_vars))
    for i in 1:length(a)
        result += a[i] * projective_vars[i]
    end
    return result - 1
end

function _bilinear_dot(a, X)
    result = zero(first(a) * first(X))
    for i in 1:length(a)
        result += a[i] * X[i]
    end
    return result
end

function _assert_not_near_zero_denominator(denominator; atol=1e-12)
    if abs(denominator) <= atol
        throw(ArgumentError("Cannot scale to affine patch because a⋅X is zero or numerically too small."))
    end
    return denominator
end

function lift_to_patch(x_affine, a)
    length(a) == length(x_affine) + 1 || throw(DimensionMismatch("a must have length length(x_affine) + 1."))
    sample = isempty(x_affine) ? first(a) : first(a) * first(x_affine)
    T = typeof(sample)
    X = Vector{T}(undef, length(a))
    X[1] = one(sample)
    for i in eachindex(x_affine)
        X[i + 1] = x_affine[i]
    end

    denominator = _assert_not_near_zero_denominator(_bilinear_dot(a, X))
    return X ./ denominator
end

function repatch(X, a)
    length(X) == length(a) || throw(DimensionMismatch("X and a must have the same length."))
    !isempty(X) || throw(ArgumentError("X must be nonempty."))

    sample = first(a) * first(X)
    X_copy = Vector{typeof(sample)}(undef, length(X))
    for i in eachindex(X)
        X_copy[i] = X[i]
    end
    denominator = _assert_not_near_zero_denominator(_bilinear_dot(a, X_copy))
    return X_copy ./ denominator
end
