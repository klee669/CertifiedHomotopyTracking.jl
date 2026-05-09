export compile_edge_homotopy, compile_homotopy

function compile_homotopy(H_eqs, x_vars, t_var; homogeneous=false)
    println("Compiling Direct Homotopy System...")
    
    target_vars = x_vars
    
    if homogeneous
        @variables u0
        target_vars = [u0; x_vars]
        H_eqs = [homogenize_poly(eq, x_vars, u0) for eq in H_eqs]
        println("-> System Homogenized. Vars: $target_vars")
    end
    
    Jx_sub = Symbolics.jacobian(H_eqs, target_vars)
    dt_sub = Symbolics.derivative.(H_eqs, t_var)
    
    compile_args = [target_vars, t_var]
    
    func_H_raw = build_function(H_eqs, compile_args...; expression=Val{false})[1]
    func_Jx_raw = build_function(Jx_sub, compile_args...; expression=Val{false})[1]
    func_dt_raw = build_function(dt_sub, compile_args...; expression=Val{false})[1]

    println("Compilation Done.")
    return CompiledHomotopy(func_H_raw, func_Jx_raw, func_dt_raw, homogeneous)
end

function compile_edge_homotopy(F_eqs, x_vars, p_vars; homogeneous=false, const_vars=Num[])
    @variables t_var
    n_params = length(p_vars)
    n_consts = length(const_vars)
    
    @variables p_starts[1:n_params] p_ends[1:n_params]
    @variables p_consts[1:n_consts]
    path_exprs = [ (1-t_var)*p_starts[i] + t_var*p_ends[i] for i in 1:n_params ]
    
    subs_dict = Dict(p_vars[i] => path_exprs[i] for i in 1:n_params)
    for i in 1:n_consts
        subs_dict[const_vars[i]] = p_consts[i]
    end
    H_sub = [Symbolics.substitute(eq, subs_dict) for eq in F_eqs]
    
    target_vars = x_vars
    
    if homogeneous
        @variables u0
        target_vars = [u0; x_vars]
        H_sub = [homogenize_poly(eq, x_vars, u0) for eq in H_sub]
        println("-> System Homogenized. Vars: $target_vars")
    end
    
    Jx_sub = Symbolics.jacobian(H_sub, target_vars)
    dt_sub = Symbolics.derivative.(H_sub, t_var)
    
    compile_args = [target_vars, t_var, p_starts..., p_ends..., p_consts...]
    
    func_H_raw = build_function(H_sub, compile_args...; expression=Val{false})[1]
    func_Jx_raw = build_function(Jx_sub, compile_args...; expression=Val{false})[1]
    func_dt_raw = build_function(dt_sub, compile_args...; expression=Val{false})[1]

    println("Compilation Done.")
    return CompiledHomotopy(func_H_raw, func_Jx_raw, func_dt_raw, homogeneous)
end
