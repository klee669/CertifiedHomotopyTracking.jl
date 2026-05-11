export compile_edge_homotopy, compile_homotopy

function _compiled_patch_template(n_vars::Int, patch_vector, patch_rng)
    raw_patch = patch_vector === nothing ?
        random_patch_vector(n_vars; rng=patch_rng) :
        patch_vector
    length(raw_patch) == n_vars || throw(DimensionMismatch("patch_vector length must match the number of projective variables."))
    return ComplexF64[ComplexF64(a) for a in raw_patch]
end

function compile_homotopy(H_eqs, x_vars, t_var; homogeneous=false, projective=false, patch_vector=nothing, patch_rng=nothing)
    println("Compiling Direct Homotopy System...")
    
    target_vars = x_vars
    use_projective_coords = homogeneous || projective
    
    if use_projective_coords
        @variables u0
        target_vars = [u0; x_vars]
        H_eqs = homogenize_system(H_eqs, x_vars, u0)
        println("-> System Homogenized. Vars: $target_vars")
    end
    
    Jx_sub = Symbolics.jacobian(H_eqs, target_vars)
    dt_sub = Symbolics.derivative.(H_eqs, t_var)
    
    compile_args = [target_vars, t_var]
    
    func_H_raw = build_function(H_eqs, compile_args...; expression=Val{false})[1]
    func_Jx_raw = build_function(Jx_sub, compile_args...; expression=Val{false})[1]
    func_dt_raw = build_function(dt_sub, compile_args...; expression=Val{false})[1]

    println("Compilation Done.")
    patch_template = projective ? _compiled_patch_template(length(target_vars), patch_vector, patch_rng) : ComplexF64[]
    return CompiledHomotopy(func_H_raw, func_Jx_raw, func_dt_raw, use_projective_coords, length(target_vars), projective, patch_template)
end

function compile_edge_homotopy(F_eqs, x_vars, p_vars; homogeneous=false, projective=false, patch_vector=nothing, patch_rng=nothing, const_vars=Num[])
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
    use_projective_coords = homogeneous || projective
    
    if use_projective_coords
        @variables u0
        target_vars = [u0; x_vars]
        H_sub = homogenize_system(H_sub, x_vars, u0)
        println("-> System Homogenized. Vars: $target_vars")
    end
    
    Jx_sub = Symbolics.jacobian(H_sub, target_vars)
    dt_sub = Symbolics.derivative.(H_sub, t_var)
    
    compile_args = [target_vars, t_var, p_starts..., p_ends..., p_consts...]
    
    func_H_raw = build_function(H_sub, compile_args...; expression=Val{false})[1]
    func_Jx_raw = build_function(Jx_sub, compile_args...; expression=Val{false})[1]
    func_dt_raw = build_function(dt_sub, compile_args...; expression=Val{false})[1]

    println("Compilation Done.")
    patch_template = projective ? _compiled_patch_template(length(target_vars), patch_vector, patch_rng) : ComplexF64[]
    return CompiledHomotopy(func_H_raw, func_Jx_raw, func_dt_raw, use_projective_coords, length(target_vars), projective, patch_template)
end
