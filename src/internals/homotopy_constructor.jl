export compile_edge_homotopy, compile_homotopy, compile_system

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

"""
    compile_homotopy(H_eqs, x_vars, t_var; projective=false) -> CompiledHomotopy

Compile an explicitly given symbolic homotopy `H_eqs(x, t)`.

Use this when you have already written the homotopy equations yourself. The
returned [`CompiledHomotopy`](@ref) can be passed directly to
[`track_path`](@ref). The tracking precision is inferred from the ACB field of
the starting point.

Complex numeric coefficients in `H_eqs` are internally parameterized and stored
as fixed constants, so expressions such as `(1 + im) * (x^3 - 1)` work.

# Options

- `projective=false`: homogenize in `x_vars` and compile coordinate-chart
  evaluators for projective tracking.

# Example

```julia
using CertifiedHomotopyTracking;

@variables x y t;
gamma = 1 + im;
compiled = compile_homotopy([(1 - t) * gamma * (x^3 - 1) + t * (x^3 + 2y - 2),
                             (1 - t) * gamma * (y^2 - 1) + t * (y^2 + x - 2)],
                            [x, y], t);
CC = AcbField(256);
res = track_path(compiled, [CC(1), CC(1)])
success(res)
solution(res)
```
"""
function compile_homotopy(H_eqs, x_vars, t_var; projective=false, patch_vector=nothing, patch_rng=nothing)
    println("Compiling Direct Homotopy System...")
    patch_vector === nothing || throw(ArgumentError("patch_vector is no longer supported; projective=true uses coordinate charts."))
    patch_rng === nothing || throw(ArgumentError("patch_rng is no longer supported; projective=true uses coordinate charts."))
    coeff_vars = Num[]
    coeff_values = Any[]
    coeff_map = Dict{Any, Num}()
    H_param = _parameterize_complex_coefficients!(collect(H_eqs), x_vars, coeff_vars, coeff_values, coeff_map)
    source = HomotopySourceData(:direct, collect(H_param), collect(x_vars), Any[], t_var, coeff_vars, projective)
    n_consts = length(coeff_vars)
    @variables p_consts[1:n_consts]
    subs_dict = Dict{Any, Any}()
    for i in 1:n_consts
        subs_dict[coeff_vars[i]] = p_consts[i]
    end
    H_sub = isempty(subs_dict) ? H_param : [Symbolics.substitute(eq, subs_dict) for eq in H_param]
    
    target_vars = x_vars
    use_projective_coords = projective
    
    if use_projective_coords
        @variables u0
        target_vars = [u0; x_vars]
        H_sub = homogenize_system(H_sub, x_vars, u0)
        println("-> System Homogenized. Vars: $target_vars")
    end
    
    compile_args = [target_vars, t_var, p_consts...]

    if projective
        chart_H, chart_Jx, chart_dt = _compile_projective_charts(H_sub, target_vars, t_var, p_consts)
        println("Compilation Done.")
        return CompiledHomotopy(chart_H[1], chart_Jx[1], chart_dt[1], chart_H, chart_Jx, chart_dt, true, length(target_vars), true, ComplexF64[], source, coeff_values)
    end
    
    Jx_sub = Symbolics.jacobian(H_sub, target_vars)
    dt_sub = Symbolics.derivative.(H_sub, t_var)
    
    func_H_raw = build_function(H_sub, compile_args...; expression=Val{false})[1]
    func_Jx_raw = build_function(Jx_sub, compile_args...; expression=Val{false})[1]
    func_dt_raw = build_function(dt_sub, compile_args...; expression=Val{false})[1]

    println("Compilation Done.")
    return CompiledHomotopy(func_H_raw, func_Jx_raw, func_dt_raw, Function[], Function[], Function[], use_projective_coords, length(target_vars), false, ComplexF64[], source, coeff_values)
end

"""
    compile_system(F_eqs, x_vars; projective=false, const_vars=Num[]) -> CompiledHomotopy

Compile a static symbolic polynomial system `F_eqs(x)`.

The result is used by [`variety_system`](@ref) and can also be evaluated through
a [`SpecializedHomotopy`](@ref) with `t = 0`.

# Options

- `projective=false`: homogenize in `x_vars` and compile projective chart evaluators.
- `const_vars=Num[]`: symbolic constants that should be supplied as fixed ACB
  values when constructing the `SpecializedHomotopy`.

# Example

```julia
using CertifiedHomotopyTracking;

@variables x y;
CC = AcbField(128);
compiled = compile_system([x^2 + y^2 - 1], [x, y]);
sys = SpecializedHomotopy(compiled, CC);
evaluate_H(sys, [CC(1), CC(0)], CC(0))
```
"""
function compile_system(F_eqs, x_vars; projective=false, patch_vector=nothing, patch_rng=nothing, const_vars=Num[])
    println("Compiling Static Polynomial System...")
    patch_vector === nothing || throw(ArgumentError("patch_vector is no longer supported; projective=true uses coordinate charts."))
    patch_rng === nothing || throw(ArgumentError("patch_rng is no longer supported; projective=true uses coordinate charts."))
    @variables t_var
    n_consts = length(const_vars)
    @variables p_consts[1:n_consts]

    subs_dict = Dict{Any, Any}()
    for i in 1:n_consts
        subs_dict[const_vars[i]] = p_consts[i]
    end
    F_sub = isempty(subs_dict) ? collect(F_eqs) : [Symbolics.substitute(eq, subs_dict) for eq in F_eqs]
    source = HomotopySourceData(:system, collect(F_eqs), collect(x_vars), Any[], nothing, collect(const_vars), projective)

    target_vars = x_vars
    use_projective_coords = projective

    if use_projective_coords
        @variables u0
        target_vars = [u0; x_vars]
        F_sub = homogenize_system(F_sub, x_vars, u0)
        println("-> System Homogenized. Vars: $target_vars")
    end

    compile_args = [target_vars, t_var, p_consts...]

    if projective
        chart_H, chart_Jx, chart_dt = _compile_projective_charts(F_sub, target_vars, t_var, p_consts)
        println("Compilation Done.")
        return CompiledHomotopy(chart_H[1], chart_Jx[1], chart_dt[1], chart_H, chart_Jx, chart_dt, true, length(target_vars), true, ComplexF64[], source)
    end

    Jx_sub = Symbolics.jacobian(F_sub, target_vars)
    dt_sub = Symbolics.derivative.(F_sub, t_var)

    func_H_raw = build_function(F_sub, compile_args...; expression=Val{false})[1]
    func_Jx_raw = build_function(Jx_sub, compile_args...; expression=Val{false})[1]
    func_dt_raw = build_function(dt_sub, compile_args...; expression=Val{false})[1]

    println("Compilation Done.")
    return CompiledHomotopy(func_H_raw, func_Jx_raw, func_dt_raw, Function[], Function[], Function[], use_projective_coords, length(target_vars), false, ComplexF64[], source)
end

"""
    compile_edge_homotopy(F_eqs, x_vars, p_vars; projective=false, const_vars=Num[]) -> CompiledHomotopy

Compile a parameter homotopy for a system `F_eqs(x, p)`.

The compiled homotopy represents the linear parameter path from `p_start` to
`p_end`. Specialize it with [`make_edge_system`](@ref), or pass it directly to
[`solve_monodromy`](@ref).

Complex numeric coefficients in `F_eqs` are internally parameterized and stored
as fixed constants, so expressions such as `(1 + 2im) * x^2` work.

# Options

- `projective=false`: homogenize in `x_vars` and use projective coordinate charts.
- `const_vars=Num[]`: symbolic constants that are not part of the moving
  parameter vector `p_vars`.

# Example

```julia
using CertifiedHomotopyTracking;

@variables x y p q;
CC = AcbField(256);
F = [p*x^2 + 3y - 4, y^2 + q];
compiled = compile_edge_homotopy(F, [x, y], [p, q]);

p0 = [CC(1), CC(-1)];
p1 = [CC(cis(0.3)), CC(cis(0.7))];
sys = make_edge_system(compiled, p0, p1);
res = track_path(sys, [CC(1), CC(1)])
success(res)
```
"""
function compile_edge_homotopy(F_eqs, x_vars, p_vars; projective=false, patch_vector=nothing, patch_rng=nothing, const_vars=Num[])
    patch_vector === nothing || throw(ArgumentError("patch_vector is no longer supported; projective=true uses coordinate charts."))
    patch_rng === nothing || throw(ArgumentError("patch_rng is no longer supported; projective=true uses coordinate charts."))
    coeff_vars = Num[]
    coeff_values = Any[]
    coeff_map = Dict{Any, Num}()
    F_param = _parameterize_complex_coefficients!(collect(F_eqs), x_vars, coeff_vars, coeff_values, coeff_map)
    all_const_vars = Num[collect(const_vars); coeff_vars]
    source = HomotopySourceData(:edge, collect(F_param), collect(x_vars), collect(p_vars), nothing, all_const_vars, projective)
    @variables t_var
    n_params = length(p_vars)
    n_consts = length(all_const_vars)
    
    @variables p_starts[1:n_params] p_ends[1:n_params]
    @variables p_consts[1:n_consts]
    path_exprs = [ (1-t_var)*p_starts[i] + t_var*p_ends[i] for i in 1:n_params ]
    
    subs_dict = Dict(p_vars[i] => path_exprs[i] for i in 1:n_params)
    for i in 1:n_consts
        subs_dict[all_const_vars[i]] = p_consts[i]
    end
    H_sub = [Symbolics.substitute(eq, subs_dict) for eq in F_param]
    
    target_vars = x_vars
    use_projective_coords = projective
    
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
        return CompiledHomotopy(chart_H[1], chart_Jx[1], chart_dt[1], chart_H, chart_Jx, chart_dt, true, length(target_vars), true, ComplexF64[], source, coeff_values)
    end

    Jx_sub = Symbolics.jacobian(H_sub, target_vars)
    dt_sub = Symbolics.derivative.(H_sub, t_var)
    
    func_H_raw = build_function(H_sub, compile_args...; expression=Val{false})[1]
    func_Jx_raw = build_function(Jx_sub, compile_args...; expression=Val{false})[1]
    func_dt_raw = build_function(dt_sub, compile_args...; expression=Val{false})[1]

    println("Compilation Done.")
    return CompiledHomotopy(func_H_raw, func_Jx_raw, func_dt_raw, Function[], Function[], Function[], use_projective_coords, length(target_vars), false, ComplexF64[], source, coeff_values)
end
