export homogenize_poly,
    homogenize_expr

function homogenize_poly(expr, vars, u0_var)
    @variables _t_deg
    sub_t = Dict(v => _t_deg * v for v in vars)
    try
        expr_t = Symbolics.substitute(expr, sub_t)
        expr_poly = Symbolics.expand(expr_t) 
        d = Symbolics.degree(expr_poly, _t_deg)
        
        sub_hom = Dict(v => v / u0_var for v in vars)
        expr_frac = Symbolics.substitute(expr, sub_hom)
        
        return expr_frac * (u0_var^d)
    catch
        sub_hom = Dict(v => v / u0_var for v in vars)
        return Symbolics.substitute(expr, sub_hom)
    end
end

function homogenize_expr(expr, vars, h_var)
    @variables _t_deg_check
    sub_deg = Dict(v => _t_deg_check * v for v in vars)
    expr_t = substitute(expr, sub_deg)
    expr_poly = Symbolics.expand(expr_t)
    
    d = Symbolics.degree(expr_poly, _t_deg_check)
    
    sub_hom = Dict(v => v / h_var for v in vars)
    expr_frac = substitute(expr, sub_hom)
    return Symbolics.expand(expr_frac * (h_var^d))
end
