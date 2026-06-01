export compile_edge_homotopy, compile_homotopy

function _projective_chart_vars(n::Int)
    return only(@variables z[1:n])
end

function _chart_substitution(projective_vars, chart_vars, chart_idx::Int)
    subs = Dict{Any,Any}()
    coord_idx = 1
    for i in eachindex(projective_vars)
        if i == chart_idx
            subs[projective_vars[i]] = 1
        else
            subs[projective_vars[i]] = chart_vars[coord_idx]
            coord_idx += 1
        end
    end
    return subs
end

function _compile_projective_charts(H_eqs, projective_vars, t_var, extra_args)
    chart_H = Function[]
    chart_Jx = Function[]
    chart_dt = Function[]
    dim = length(projective_vars) - 1

    for chart_idx in eachindex(projective_vars)
        chart_vars = _projective_chart_vars(dim)
        subs = _chart_substitution(projective_vars, chart_vars, chart_idx)
        H_chart = [Symbolics.substitute(eq, subs) for eq in H_eqs]
        Jx_chart = Symbolics.jacobian(H_chart, chart_vars)
        dt_chart = Symbolics.derivative.(H_chart, t_var)
        compile_args = [chart_vars, t_var, extra_args...]

        push!(chart_H, build_function(H_chart, compile_args...; expression=Val{false})[1])
        push!(chart_Jx, build_function(Jx_chart, compile_args...; expression=Val{false})[1])
        push!(chart_dt, build_function(dt_chart, compile_args...; expression=Val{false})[1])
    end

    return chart_H, chart_Jx, chart_dt
end

function compile_homotopy(H_eqs, x_vars, t_var; homogeneous=false, projective=false, patch_vector=nothing, patch_rng=nothing)
    println("Compiling Direct Homotopy System...")
    patch_vector === nothing || throw(ArgumentError("patch_vector is no longer supported; projective=true uses coordinate charts."))
    patch_rng === nothing || throw(ArgumentError("patch_rng is no longer supported; projective=true uses coordinate charts."))
    
    target_vars = x_vars
    use_projective_coords = homogeneous || projective
    
    if use_projective_coords
        @variables u0
        target_vars = [u0; x_vars]
        H_eqs = homogenize_system(H_eqs, x_vars, u0)
        println("-> System Homogenized. Vars: $target_vars")
    end
    
    if projective
        chart_H, chart_Jx, chart_dt = _compile_projective_charts(H_eqs, target_vars, t_var, Any[])
        println("Compilation Done.")
        return CompiledHomotopy(chart_H[1], chart_Jx[1], chart_dt[1], chart_H, chart_Jx, chart_dt, true, length(target_vars), true, ComplexF64[])
    end
    
    Jx_sub = Symbolics.jacobian(H_eqs, target_vars)
    dt_sub = Symbolics.derivative.(H_eqs, t_var)
    compile_args = [target_vars, t_var]
    
    func_H_raw = build_function(H_eqs, compile_args...; expression=Val{false})[1]
    func_Jx_raw = build_function(Jx_sub, compile_args...; expression=Val{false})[1]
    func_dt_raw = build_function(dt_sub, compile_args...; expression=Val{false})[1]

    println("Compilation Done.")
    return CompiledHomotopy(func_H_raw, func_Jx_raw, func_dt_raw, Function[], Function[], Function[], use_projective_coords, length(target_vars), false, ComplexF64[])
end

function compile_edge_homotopy(F_eqs, x_vars, p_vars; homogeneous=false, projective=false, patch_vector=nothing, patch_rng=nothing, const_vars=Num[])
    patch_vector === nothing || throw(ArgumentError("patch_vector is no longer supported; projective=true uses coordinate charts."))
    patch_rng === nothing || throw(ArgumentError("patch_rng is no longer supported; projective=true uses coordinate charts."))
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
    
    compile_args = [target_vars, t_var, p_starts..., p_ends..., p_consts...]

    if projective
        chart_H, chart_Jx, chart_dt = _compile_projective_charts(H_sub, target_vars, t_var, [p_starts..., p_ends..., p_consts...])
        println("Compilation Done.")
        return CompiledHomotopy(chart_H[1], chart_Jx[1], chart_dt[1], chart_H, chart_Jx, chart_dt, true, length(target_vars), true, ComplexF64[])
    end

    Jx_sub = Symbolics.jacobian(H_sub, target_vars)
    dt_sub = Symbolics.derivative.(H_sub, t_var)
    
    func_H_raw = build_function(H_sub, compile_args...; expression=Val{false})[1]
    func_Jx_raw = build_function(Jx_sub, compile_args...; expression=Val{false})[1]
    func_dt_raw = build_function(dt_sub, compile_args...; expression=Val{false})[1]

    println("Compilation Done.")
    return CompiledHomotopy(func_H_raw, func_Jx_raw, func_dt_raw, Function[], Function[], Function[], use_projective_coords, length(target_vars), false, ComplexF64[])
end
