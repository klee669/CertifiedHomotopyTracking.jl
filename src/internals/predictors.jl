export speed_vector,
       linear_predictor,
       hermite_predictor,
       compute_velocity,
       construct_hermite_predictor_tm


# --------------------------------------------------------------------------
# Compute the speed vector for the linear predictor
function speed_vector(
    H::Union{Matrix, Vector}, 
    x::Vector{AcbFieldElem}, 
    tval::Number, 
    A::AcbMatrix,
)
    ring = base_ring(H[1])
    n    = length(H)

    # Evaluate dH/dt symbolically at t
    result = evaluate_matrix(Matrix(transpose(hcat(H))),tval);
    for i = 1:n
        result[i]=Nemo.evaluate(derivative(H[i]),tval);
    end
    result = evaluate_matrix(result, x);

    midpoint_complex_box(-A*transpose(hcat(result)))
end


# --------------------------------------------------------------------------
# Linear predictor step: x + η·v
function linear_predictor(
    H::Union{Matrix, Vector}, 
    v::Vector{AcbFieldElem}, 
    x::Vector{AcbFieldElem},
)
    eR    = base_ring(H[1])
    η     = gens(eR)[end]
    delta = [vi * η for vi in v]

    return x + delta
end


# --------------------------------------------------------------------------
# Hermite predictor (uses previous step data)
function hermite_predictor(
    H::Union{Matrix, Vector}, 
    x::Vector, 
    xprev::Vector, 
    v::Vector, 
    vprev::Vector, 
    hprev::Number,
)
    eR = base_ring(H[1])
    η  = gens(eR)[end]
    n  = length(v)

    result = Vector{typeof(η)}(undef, n)

    for i in 1:n
        dη   = v[i] * η
        η²   = η^2
        η³   = η^3

        dx   = x[i] - xprev[i]
        dv   = v[i] - vprev[i]
        h²   = hprev^2
        h³   = hprev^3

        t2 = (3 * v[i] / hprev - dv / hprev - 3 * dx / h²) * η²
        t3 = (2 * v[i] / h²   - dv / h²   - 2 * dx / h³) * η³

        result[i] = dη + t2 + t3
    end

    return x + result
end


function compute_velocity(sys::HCSystem, x::AbstractVector{AcbFieldElem}, t, A::AbstractMatrix{AcbFieldElem})
    CC = sys.CC
    x_mid = get_mid_vec(x)
    t_mid = t isa AcbFieldElem ? get_mid(t) : CC(t)
    Ht = evaluate_dt(sys, x_mid, t_mid)
    return -A * Ht
end

function construct_hermite_predictor_tm(sys::HCSystem, x, x_prev, v, v_prev, h_prev, h_curr)
    CC = sys.CC; RR = sys.RR
    dx = x - x_prev
    dv = v - v_prev
    
    h_p = CC(h_prev)
    h_sq = h_p^2
    h_cb = h_p^3
    
    c0 = x
    c1 = v
    c2 = (3 .* v ./ h_p) .- (dv ./ h_p) .- (3 .* dx ./ h_sq)
    c3 = (2 .* v ./ h_sq) .- (dv ./ h_sq) .- (2 .* dx ./ h_cb)
    
    n = length(x)
    h_arb = RR(h_curr) 
    
    tm_vec = Vector{TaylorModel3}(undef, n)
    for i in 1:n
        tm_vec[i] = TaylorModel3(c0[i], c1[i], c2[i], c3[i], CC(0), h_arb)
    end
    return tm_vec
end