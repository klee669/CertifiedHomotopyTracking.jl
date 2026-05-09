const _complex_coeff_counter = Ref(0)

function _next_complex_coeff_variable()
    _complex_coeff_counter[] += 1
    return Symbolics.variable(Symbol("__complex_coeff_", _complex_coeff_counter[], "__"))
end

function _symbolic_const_value(x)
    x isa Number && !(x isa Num) && return x
    isempty(Symbolics.get_variables(x)) || return nothing
    u = Symbolics.unwrap(x)
    return Symbolics.SymbolicUtils.isconst(u) ? Symbolics.SymbolicUtils.unwrap_const(u) : nothing
end

function _is_symbolic_zero(x)
    value = _symbolic_const_value(x)
    return value !== nothing && iszero(value)
end

function _complex_coeff_value(re_coeff, im_coeff)
    re_value = _symbolic_const_value(re_coeff)
    im_value = _symbolic_const_value(im_coeff)
    if re_value === nothing || im_value === nothing
        error("Complex symbolic coefficients must be numeric constants to be converted into fixed parameters.")
    end
    return complex(re_value, im_value)
end

function _fixed_complex_coeff_parameter(coeff, coeff_vars, coeff_values, coeff_map)
    if haskey(coeff_map, coeff)
        return coeff_map[coeff]
    end

    var = _next_complex_coeff_variable()
    coeff_map[coeff] = var
    push!(coeff_vars, var)
    push!(coeff_values, coeff)
    return var
end

_as_num(x) = x isa Num ? x : Num(x)

function _split_numeric_factor(term)
    value = _symbolic_const_value(term)
    value !== nothing && return value, Num(1)

    term_num = _as_num(term)
    unwrapped = Symbolics.unwrap(term_num)
    op = try
        Symbolics.operation(unwrapped)
    catch
        return 1, term_num
    end
    args = Symbolics.arguments(unwrapped)

    if op === (*)
        coeff = 1
        factor = Num(1)
        for arg in args
            arg_coeff, arg_factor = _split_numeric_factor(_as_num(arg))
            coeff *= arg_coeff
            factor *= arg_factor
        end
        return coeff, Symbolics.expand(factor)
    elseif op === (/)
        num_coeff, num_factor = _split_numeric_factor(_as_num(args[1]))
        den_coeff, den_factor = _split_numeric_factor(_as_num(args[2]))
        return num_coeff / den_coeff, Symbolics.expand(num_factor / den_factor)
    else
        return 1, term_num
    end
end

function _numeric_coeff_value(coeff)
    value = _symbolic_const_value(coeff)
    value !== nothing && return value
    error("Complex symbolic coefficients must be numeric constants; got coefficient $coeff.")
end

function _term_coefficients(expr)
    terms = Dict{String, Tuple{Num, Any}}()
    for term in Symbolics.terms(Symbolics.expand(expr))
        coeff, factor = _split_numeric_factor(term)
        coeff = _numeric_coeff_value(coeff)
        key = string(factor)
        if haskey(terms, key)
            old_factor, old_coeff = terms[key]
            terms[key] = (old_factor, old_coeff + coeff)
        else
            terms[key] = (factor, coeff)
        end
    end
    return terms
end

function _parameterize_complex_coefficients(expr, x_vars, coeff_vars, coeff_values, coeff_map)
    expr isa Num && return expr
    expr isa Complex{Num} || error("Expected a Symbolics Num or Complex{Num} expression, got $(typeof(expr)).")

    re_terms = _term_coefficients(real(expr))
    im_terms = _term_coefficients(imag(expr))
    factor_keys = collect(union(keys(re_terms), keys(im_terms)))
    sort!(factor_keys)

    out = Num(0)
    for key in factor_keys
        factor = haskey(re_terms, key) ? re_terms[key][1] : im_terms[key][1]
        re_coeff = haskey(re_terms, key) ? re_terms[key][2] : 0
        im_coeff = haskey(im_terms, key) ? im_terms[key][2] : 0

        if iszero(im_coeff)
            out += re_coeff * factor
        else
            coeff = complex(re_coeff, im_coeff)
            coeff_param = _fixed_complex_coeff_parameter(coeff, coeff_vars, coeff_values, coeff_map)
            out += coeff_param * factor
        end
    end

    return Symbolics.expand(out)
end

function _parameterize_complex_coefficients(exprs::AbstractVector, x_vars)
    coeff_vars = Num[]
    coeff_values = Any[]
    coeff_map = Dict{Any, Num}()

    new_exprs = [
        _parameterize_complex_coefficients(expr, x_vars, coeff_vars, coeff_values, coeff_map)
        for expr in exprs
    ]

    return new_exprs, coeff_vars, coeff_values, coeff_map
end

function _parameterize_complex_coefficients!(exprs::AbstractVector, x_vars, coeff_vars, coeff_values, coeff_map)
    return [
        _parameterize_complex_coefficients(expr, x_vars, coeff_vars, coeff_values, coeff_map)
        for expr in exprs
    ]
end

function _coefficient_value(CCRing, coeff)
    coeff isa Complex ? CCRing(real(coeff), imag(coeff)) : CCRing(coeff)
end

