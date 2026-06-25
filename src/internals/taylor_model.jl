function _tm_real_interval(CC::AcbField, h::ArbFieldElem)
    RR = parent(h)
    upper = Nemo.midpoint(h) + Nemo.radius(h)
    half_width = RR(upper / 2)
    real_interval = Nemo.ball(half_width, half_width)
    return CC(real_interval, RR(0))
end

"""
    TaylorModel3

Third-order Taylor model with an interval remainder over a step interval `h`.

Fields `c0`, `c1`, `c2`, `c3` store polynomial coefficients and `rem` stores an
ACB remainder bound. This type is primarily used by certified Hermite predictor
validation in [`track_path`](@ref).
"""
struct TaylorModel3{C,R}
    c0::C
    c1::C
    c2::C
    c3::C
    rem::C
    h::R
end

_tm3(c0::AcbFieldElem, c1::AcbFieldElem, c2::AcbFieldElem, c3::AcbFieldElem, rem::AcbFieldElem, h::ArbFieldElem) =
    TaylorModel3{AcbFieldElem,ArbFieldElem}(c0, c1, c2, c3, rem, h)

TaylorModel3(c0::AcbFieldElem, c1::AcbFieldElem, c2::AcbFieldElem, c3::AcbFieldElem, rem::AcbFieldElem, h::ArbFieldElem) =
    _tm3(c0, c1, c2, c3, rem, h)

function TaylorModel3(c0, c1, c2, c3, rem::AcbFieldElem, h::ArbFieldElem)
    CC = parent(rem)
    TaylorModel3(CC(c0), CC(c1), CC(c2), CC(c3), rem, h)
end

function TaylorModel3(c::NTuple{4}, rem::AcbFieldElem, h::ArbFieldElem)
    TaylorModel3(c[1], c[2], c[3], c[4], rem, h)
end

function TaylorModel3(c::AbstractVector, rem::AcbFieldElem, h::ArbFieldElem)
    length(c) == 4 || throw(DimensionMismatch("TaylorModel3 requires exactly four coefficients."))
    TaylorModel3(c[1], c[2], c[3], c[4], rem, h)
end

Base.zero(tm::TaylorModel3) = TaylorModel3(parent(tm.rem)(0), parent(tm.rem)(0), parent(tm.rem)(0), parent(tm.rem)(0), parent(tm.rem)(0), parent(tm.h)(0))
Base.one(tm::TaylorModel3) = TaylorModel3(parent(tm.rem)(1), parent(tm.rem)(0), parent(tm.rem)(0), parent(tm.rem)(0), parent(tm.rem)(0), tm.h)
Base.Broadcast.broadcastable(x::TaylorModel3) = Ref(x)
Base.getproperty(tm::TaylorModel3, s::Symbol) = s === :c ? (getfield(tm, :c0), getfield(tm, :c1), getfield(tm, :c2), getfield(tm, :c3)) : getfield(tm, s)

function evaluate_taylor!(res::AcbFieldElem, tm::TaylorModel3, cache::TMCache)
    CC = parent(res)
    RR = parent(tm.h)
    
    t_int = _tm_real_interval(CC, tm.h)
    
    Nemo.add!(res, cache.zero_cc, tm.c3)
    
    Nemo.mul!(cache.temp1, t_int, res)
    Nemo.add!(res, tm.c2, cache.temp1)
    
    Nemo.mul!(cache.temp1, t_int, res)
    Nemo.add!(res, tm.c1, cache.temp1)
    
    Nemo.mul!(cache.temp1, t_int, res)
    Nemo.add!(res, tm.c0, cache.temp1)
    
    Nemo.add!(res, res, tm.rem)
    
    return res
end

"""
    evaluate_taylor(tm::TaylorModel3)

Evaluate a third-order Taylor model on its stored interval `h`, including the
remainder bound.
"""
function evaluate_taylor(tm::TaylorModel3)
    CC = parent(tm.rem); RR = parent(tm.h)
    t_int = _tm_real_interval(CC, tm.h)
    return _evaluate_poly(tm.c0, tm.c1, tm.c2, tm.c3, t_int) + tm.rem
end

_evaluate_poly(c0, c1, c2, c3, t) = c0 + t * (c1 + t * (c2 + t * c3))

function Base.:+(a::TaylorModel3, b::TaylorModel3)
    _tm3(a.c0 + b.c0, a.c1 + b.c1, a.c2 + b.c2, a.c3 + b.c3, a.rem + b.rem, a.h)
end
function Base.:+(a::TaylorModel3, b::AcbFieldElem)
    _tm3(a.c0 + b, a.c1, a.c2, a.c3, a.rem, a.h)
end
function Base.:+(a::TaylorModel3, b::Number)
    CC = parent(a.rem)
    _tm3(a.c0 + CC(b), a.c1, a.c2, a.c3, a.rem, a.h)
end
Base.:+(b::Union{Number, AcbFieldElem}, a::TaylorModel3) = a + b

function Base.:-(a::TaylorModel3, b::TaylorModel3)
    _tm3(a.c0 - b.c0, a.c1 - b.c1, a.c2 - b.c2, a.c3 - b.c3, a.rem - b.rem, a.h)
end
Base.:-(a::TaylorModel3, b::Union{Number, AcbFieldElem}) = a + (-b)
function Base.:-(b::AcbFieldElem, a::TaylorModel3)
    _tm3(b - a.c0, -a.c1, -a.c2, -a.c3, -a.rem, a.h)
end
function Base.:-(b::Number, a::TaylorModel3)
    CC = parent(a.rem)
    _tm3(CC(b) - a.c0, -a.c1, -a.c2, -a.c3, -a.rem, a.h)
end

function Base.:*(a::TaylorModel3, b::TaylorModel3)
    CC = parent(a.rem); RR = parent(a.h)
    h = a.h
    t_int = _tm_real_interval(CC, h)
    
    C0 = a.c0*b.c0
    C1 = a.c0*b.c1 + a.c1*b.c0
    C2 = a.c0*b.c2 + a.c1*b.c1 + a.c2*b.c0
    C3 = a.c0*b.c3 + a.c1*b.c2 + a.c2*b.c1 + a.c3*b.c0
    
    term4 = a.c1*b.c3 + a.c2*b.c2 + a.c3*b.c1
    term5 = a.c2*b.c3 + a.c3*b.c2
    term6 = a.c3*b.c3
    
    t2 = t_int * t_int
    t4 = t2 * t2
    t5 = t4 * t_int
    t6 = t4 * t2
    trunc_error = term4 * t4 + term5 * t5 + term6 * t6
    
    pA = _evaluate_poly(a.c0, a.c1, a.c2, a.c3, t_int)
    pB = _evaluate_poly(b.c0, b.c1, b.c2, b.c3, t_int)
    
    prop_error = pA * b.rem + pB * a.rem + a.rem * b.rem
    new_rem = trunc_error + prop_error
    
    _tm3(C0, C1, C2, C3, new_rem, h)
end

function Base.:*(a::TaylorModel3, b::AcbFieldElem)
    CC = parent(a.rem)
    m = get_mid(b)
    dev = b - m
    t_int = _tm_real_interval(CC, a.h)
    
    pA = _evaluate_poly(a.c0, a.c1, a.c2, a.c3, t_int)
    new_rem = a.rem * b + pA * dev
    
    _tm3(a.c0 * m, a.c1 * m, a.c2 * m, a.c3 * m, new_rem, a.h)
end

function Base.:*(a::TaylorModel3, b::Number)
    CC = parent(a.rem)
    val = CC(b)
    m = val
    new_rem = a.rem * val
    
    _tm3(a.c0 * m, a.c1 * m, a.c2 * m, a.c3 * m, new_rem, a.h)
end
Base.:*(b::Union{Number, AcbFieldElem}, a::TaylorModel3) = a * b

function Base.:^(a::TaylorModel3, n::Integer)
    if n == 0 return one(a) end
    if n == 1 return a end
    res = a
    for i in 2:n
        res = res * a
    end
    return res
end
Base.literal_pow(::typeof(^), a::TaylorModel3, ::Val{0}) = one(a)
Base.literal_pow(::typeof(^), a::TaylorModel3, ::Val{1}) = a
Base.literal_pow(::typeof(^), a::TaylorModel3, ::Val{2}) = a * a
Base.literal_pow(::typeof(^), a::TaylorModel3, ::Val{3}) = (a * a) * a

function Base.:/(a::TaylorModel3, b::TaylorModel3)
    CC = parent(a.rem)
    b0 = b.c0
    if contains(abs(b0), 0)
        error("Division by zero in Taylor Model")
    end
    inv_b0 = 1 / b0
    
    C0 = a.c0 * inv_b0
    val2 = a.c1 - C0*b.c1; C1 = val2 * inv_b0
    val3 = a.c2 - (C0*b.c2 + C1*b.c1); C2 = val3 * inv_b0
    val4 = a.c3 - (C0*b.c3 + C1*b.c2 + C2*b.c1); C3 = val4 * inv_b0
    
    tm_H_poly = _tm3(C0, C1, C2, C3, CC(0), a.h)
    tm_residue = a - (tm_H_poly * b)
    range_B = evaluate_taylor(b)
    
    if contains(range_B, 0)
        error("Division by zero in interval evaluation")
    end
    
    new_rem = evaluate_taylor(tm_residue) / range_B
    return _tm3(C0, C1, C2, C3, new_rem, a.h)
end

function Base.:/(a::TaylorModel3, b::Union{Number, AcbFieldElem})
    CC = parent(a.rem)
    inv_b = 1 / CC(b)
    return a * inv_b
end

function Base.:/(a::Union{Number, AcbFieldElem}, b::TaylorModel3)
    CC = parent(b.rem)
    tm_a = TaylorModel3(CC(a), CC(0), CC(0), CC(0), CC(0), b.h)
    return tm_a / b
end
