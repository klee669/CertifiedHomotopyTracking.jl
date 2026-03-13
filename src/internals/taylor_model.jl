export TaylorModel3, evaluate_taylor

struct TaylorModel3
    c::Vector{AcbFieldElem}
    rem::AcbFieldElem       
    h::ArbFieldElem         
end

function TaylorModel3(c0, c1, c2, c3, rem::AcbFieldElem, h::ArbFieldElem)
    CC = parent(rem)
    TaylorModel3([CC(c0), CC(c1), CC(c2), CC(c3)], rem, h)
end

Base.zero(tm::TaylorModel3) = TaylorModel3(parent(tm.rem)(0), parent(tm.rem)(0), parent(tm.rem)(0), parent(tm.rem)(0), parent(tm.rem)(0), parent(tm.h)(0))
Base.one(tm::TaylorModel3) = TaylorModel3(parent(tm.rem)(1), parent(tm.rem)(0), parent(tm.rem)(0), parent(tm.rem)(0), parent(tm.rem)(0), tm.h)
Base.Broadcast.broadcastable(x::TaylorModel3) = Ref(x)

function evaluate_taylor!(res::AcbFieldElem, tm::TaylorModel3, cache::TMCache)
    CC = parent(res)
    RR = parent(tm.h)
    
    t_int = CC(RR(0), tm.h) 
    
    Nemo.add!(res, cache.zero_cc, tm.c[4])
    
    Nemo.mul!(cache.temp1, t_int, res)
    Nemo.add!(res, tm.c[3], cache.temp1)
    
    Nemo.mul!(cache.temp1, t_int, res)
    Nemo.add!(res, tm.c[2], cache.temp1)
    
    Nemo.mul!(cache.temp1, t_int, res)
    Nemo.add!(res, tm.c[1], cache.temp1)
    
    Nemo.add!(res, res, tm.rem)
    
    return res
end

function evaluate_taylor(tm::TaylorModel3)
    CC = parent(tm.rem); RR = parent(tm.h)
    t_int = CC(RR(0), tm.h) 
    val = tm.c[4]
    val = tm.c[3] + t_int * val
    val = tm.c[2] + t_int * val
    val = tm.c[1] + t_int * val
    return val + tm.rem
end

function Base.:+(a::TaylorModel3, b::TaylorModel3)
    TaylorModel3(a.c .+ b.c, a.rem + b.rem, a.h)
end
function Base.:+(a::TaylorModel3, b::Union{Number, AcbFieldElem})
    CC = parent(a.rem)
    new_c = copy(a.c)
    new_c[1] += CC(b)
    TaylorModel3(new_c, a.rem, a.h)
end
Base.:+(b::Union{Number, AcbFieldElem}, a::TaylorModel3) = a + b

function Base.:-(a::TaylorModel3, b::TaylorModel3)
    TaylorModel3(a.c .- b.c, a.rem - b.rem, a.h)
end
Base.:-(a::TaylorModel3, b::Union{Number, AcbFieldElem}) = a + (-b)
function Base.:-(b::Union{Number, AcbFieldElem}, a::TaylorModel3)
    CC = parent(a.rem)
    new_c = -a.c
    new_c[1] += CC(b)
    TaylorModel3(new_c, -a.rem, a.h)
end

function Base.:*(a::TaylorModel3, b::TaylorModel3)
    CC = parent(a.rem); RR = parent(a.h)
    h = a.h
    t_int = CC(RR(0), h)
    
    A = a.c; B = b.c
    C = [CC(0) for _ in 1:4]
    
    C[1] = A[1]*B[1]
    C[2] = A[1]*B[2] + A[2]*B[1]
    C[3] = A[1]*B[3] + A[2]*B[2] + A[3]*B[1]
    C[4] = A[1]*B[4] + A[2]*B[3] + A[3]*B[2] + A[4]*B[1]
    
    term4 = A[2]*B[4] + A[3]*B[3] + A[4]*B[2]
    term5 = A[3]*B[4] + A[4]*B[3]
    term6 = A[4]*B[4]
    
    trunc_error = term4 * (t_int^4) + term5 * (t_int^5) + term6 * (t_int^6)
    
    pA = evaluate_taylor(TaylorModel3(A[1], A[2], A[3], A[4], CC(0), h))
    pB = evaluate_taylor(TaylorModel3(B[1], B[2], B[3], B[4], CC(0), h))
    
    prop_error = pA * b.rem + pB * a.rem + a.rem * b.rem
    new_rem = trunc_error + prop_error
    
    TaylorModel3(C, new_rem, h)
end

function Base.:*(a::TaylorModel3, b::Union{Number, AcbFieldElem})
    CC = parent(a.rem)
    val = CC(b)
    m = get_mid(val)
    dev = val - m
    
    new_c = a.c .* m
    pA = evaluate_taylor(TaylorModel3(a.c[1], a.c[2], a.c[3], a.c[4], CC(0), a.h))
    new_rem = a.rem * val + pA * dev
    
    TaylorModel3(new_c, new_rem, a.h)
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

function Base.:/(a::TaylorModel3, b::TaylorModel3)
    CC = parent(a.rem)
    b0 = b.c[1]
    if contains(abs(b0), 0)
        error("Division by zero in Taylor Model")
    end
    inv_b0 = 1 / b0
    
    A = a.c; B = b.c
    C = Vector{AcbFieldElem}(undef, 4)
    
    C[1] = A[1] * inv_b0
    val2 = A[2] - C[1]*B[2]; C[2] = val2 * inv_b0
    val3 = A[3] - (C[1]*B[3] + C[2]*B[2]); C[3] = val3 * inv_b0
    val4 = A[4] - (C[1]*B[4] + C[2]*B[3] + C[3]*B[2]); C[4] = val4 * inv_b0
    
    tm_H_poly = TaylorModel3(C[1], C[2], C[3], C[4], CC(0), a.h)
    tm_residue = a - (tm_H_poly * b)
    range_B = evaluate_taylor(b)
    
    if contains(range_B, 0)
        error("Division by zero in interval evaluation")
    end
    
    new_rem = evaluate_taylor(tm_residue) / range_B
    return TaylorModel3(C, new_rem, a.h)
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