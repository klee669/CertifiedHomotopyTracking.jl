export CompiledHomotopy, compile_edge_homotopy, HCSystem, evaluate_H, evaluate_Jac, evaluate_dt

struct CompiledHomotopy
    func_H::Function
    func_Jx::Function
    func_dt::Function
    homogeneous::Bool
end

mutable struct HCSystem
    compiled::CompiledHomotopy
    p_start::Tuple{Vararg{AcbFieldElem}}
    p_end::Tuple{Vararg{AcbFieldElem}}
    homogeneous::Bool
    patch_idx::Int
    
    CC::AcbField 
    RR::ArbField
    
    function HCSystem(compiled::CompiledHomotopy, p_start::Vector{AcbFieldElem}, p_end::Vector{AcbFieldElem})
        CC = parent(p_start[1]) 
        RR = ArbField(precision(CC))
        new(compiled, Tuple(p_start), Tuple(p_end), compiled.homogeneous, 1, CC, RR)
    end
end

struct TMCache
    temp1::AcbFieldElem
    temp2::AcbFieldElem
    temp3::AcbFieldElem
    zero_cc::AcbFieldElem 
    
    function TMCache(CC::AcbField)
        new(CC(0), CC(0), CC(0), CC(0))
    end
end


function compile_edge_homotopy(F_eqs, x_vars, p_vars, t_var; homogeneous=false)
    println("Compiling Homotopy System (Only Once!)...")
    n_params = length(p_vars)
    
    @variables p_starts[1:n_params] p_ends[1:n_params]
    path_exprs = [ (1-t_var)*p_starts[i] + t_var*p_ends[i] for i in 1:n_params ]
    
    subs_dict = Dict(p_vars[i] => path_exprs[i] for i in 1:n_params)
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
    
    compile_args = [target_vars, t_var, p_starts..., p_ends...]
    
    func_H_raw = build_function(H_sub, compile_args...; expression=Val{false})[1]
    func_Jx_raw = build_function(Jx_sub, compile_args...; expression=Val{false})[1]
    func_dt_raw = build_function(dt_sub, compile_args...; expression=Val{false})[1]

    println("Compilation Done.")
    return CompiledHomotopy(func_H_raw, func_Jx_raw, func_dt_raw, homogeneous)
end

function evaluate_H_augmented(sys::HCSystem, x, t)
    val_sys = sys.compiled.func_H(x, t, sys.p_start..., sys.p_end...)
    if !sys.homogeneous
        return val_sys
    else
        val_patch = x[sys.patch_idx] - 1
        return [val_sys; val_patch]
    end
end
evaluate_H(sys::HCSystem, x, t) = evaluate_H_augmented(sys, x, t)

function evaluate_Jac(sys::HCSystem, x, t)
    CC = sys.CC 
    J_sys = sys.compiled.func_Jx(x, t, sys.p_start..., sys.p_end...)
    if !sys.homogeneous
        return J_sys
    else
        n_rows, n_cols = size(J_sys)
        J_aug = Matrix{AcbFieldElem}(undef, n_rows + 1, n_cols)
        for i in 1:n_rows, j in 1:n_cols
            J_aug[i,j] = CC(J_sys[i,j])
        end
        for j in 1:n_cols
            J_aug[n_rows+1, j] = (j == sys.patch_idx) ? CC(1) : CC(0)
        end
        return J_aug
    end
end

function evaluate_dt_augmented(sys::HCSystem, x, t)
    CC = sys.CC 
    val_dt = sys.compiled.func_dt(x, t, sys.p_start..., sys.p_end...)
    if !sys.homogeneous
        return val_dt
    else
        return [val_dt; CC(0)]
    end
end
evaluate_dt(sys::HCSystem, x, t) = evaluate_dt_augmented(sys, x, t)
