export straight_line_homotopy, linear_path, specified_system, max_degree

# ------------------------------------------------------------------------------
# Straight Line Homotopy
# ------------------------------------------------------------------------------

_straight_line_field_value(CC, value) = value isa Complex ? CC(real(value), imag(value)) : CC(value)

function straight_line_homotopy(F, G, t; gamma=nothing)
    n = length(F)
    
    R = parent(F[1])
    CC = base_ring(R) 
    
    gamma_val = gamma === nothing ? CC(rand(Float64), rand(Float64)) : _straight_line_field_value(CC, gamma)
    
    HR = parent(t)
    H = zeros(HR, 0)
    for i in 1:n
        H = push!(H, (1 - t) * gamma_val * G[i] + t * F[i])
    end
    H
end
function straight_line_homotopy(F_exprs::AbstractVector{<:Union{Num, Complex{Num}}}, 
    G_exprs::AbstractVector{<:Union{Num, Complex{Num}}}, 
    x_vars::AbstractVector{Num}; 
    CCRing=AcbField(256),
    homogeneous=false,
    projective=false,
    patch_vector=nothing,
    patch_rng=nothing,
    gamma=nothing,
)
    @variables __gamma_trick_internal_param__
    @variables t_var

    coeff_vars = Num[]
    coeff_values = Any[]
    coeff_map = Dict{Any, Num}()
    F_param = _parameterize_complex_coefficients!(F_exprs, x_vars, coeff_vars, coeff_values, coeff_map)
    G_param = _parameterize_complex_coefficients!(G_exprs, x_vars, coeff_vars, coeff_values, coeff_map)
    
    H = [(1 - t_var) * __gamma_trick_internal_param__ * G_param[i] + t_var * F_param[i] for i in 1:length(F_param)]
    compiled_H = compile_edge_homotopy(
        H,
        x_vars,
        [__gamma_trick_internal_param__];
        homogeneous=homogeneous,
        projective=projective,
        patch_vector=patch_vector,
        patch_rng=patch_rng,
        const_vars=coeff_vars,
    )
    
    gamma_val = gamma === nothing ? CCRing(complex(randn(), randn())) : _coefficient_value(CCRing, gamma)
    coeff_vals = [_coefficient_value(CCRing, coeff) for coeff in coeff_values]
    sys = make_edge_system(compiled_H, [gamma_val], [gamma_val], coeff_vals)
    
    return sys
end
# ------------------------------------------------------------------------------
# Max Degree Calculation
# ------------------------------------------------------------------------------
function max_degree(H::Matrix)
    HR = parent(H[1])       # R[t]
    R  = base_ring(HR)      # R
    CC = base_ring(R)       # CC (AcbField)
    
    rand_t = CC(rand(ComplexF64))
    generic_F = evaluate_matrix(H, rand_t)

    return map(i -> maximum(map(j -> sum(degrees(j)), AbstractAlgebra.monomials(generic_F[i]))), 1:length(generic_F))
end

# ------------------------------------------------------------------------------
# Linear Path Construction
# ------------------------------------------------------------------------------
function linear_path(p0, p1, t)
    n = length(p0)
    HR = parent(t)
    p = zeros(HR, n)
    for i = 1:n
        p[i] = p0[i] * (1 - t) + p1[i] * t
    end
    p
end

# ------------------------------------------------------------------------------
# Specified System Construction
# ------------------------------------------------------------------------------
function specified_system(p0, p1, F)
    
    HR = base_ring(F[1]) 
    
    n = length(F)
    t = gens(HR)[1] 
    
    p = linear_path(p0, p1, t)

    Fp = zeros(HR, n)
    for i in 1:n
        Fp[i] = AbstractAlgebra.evaluate(F[i], p)
    end
    hcat(Fp)
end
