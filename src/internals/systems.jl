export CompiledHomotopy, HCSystem, TMCache, evaluate_H, evaluate_Jac, evaluate_dt

struct CompiledHomotopy
    func_H::Function
    func_Jx::Function
    func_dt::Function
    homogeneous::Bool
end

function Base.show(io::IO, compiled::CompiledHomotopy)
    print(io, "CompiledHomotopy(Homogeneous: ", compiled.homogeneous, ")")
end

mutable struct HCSystem
    compiled::CompiledHomotopy
    p_start::Tuple{Vararg{AcbFieldElem}}
    p_end::Tuple{Vararg{AcbFieldElem}}
    p_const::Tuple{Vararg{AcbFieldElem}}
    homogeneous::Bool
    patch_idx::Int
    
    CC::AcbField 
    RR::ArbField
    
    function HCSystem(compiled::CompiledHomotopy, p_start::Vector{AcbFieldElem}, p_end::Vector{AcbFieldElem}, p_const::Vector{AcbFieldElem}=AcbFieldElem[])
        CC = !isempty(p_start) ? parent(p_start[1]) : parent(p_const[1])
        RR = ArbField(precision(CC))
        new(compiled, Tuple(p_start), Tuple(p_end), Tuple(p_const), compiled.homogeneous, 1, CC, RR)
    end

    function HCSystem(compiled::CompiledHomotopy, CC::AcbField)
        RR = ArbField(precision(CC))
        new(compiled, (), (), (), compiled.homogeneous, 1, CC, RR)
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

function evaluate_H_augmented(sys::HCSystem, x, t)
    val_sys = sys.compiled.func_H(x, t, sys.p_start..., sys.p_end..., sys.p_const...)
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
    J_sys = sys.compiled.func_Jx(x, t, sys.p_start..., sys.p_end..., sys.p_const...)
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
    val_dt = sys.compiled.func_dt(x, t, sys.p_start..., sys.p_end..., sys.p_const...)
    if !sys.homogeneous
        return val_dt
    else
        return [val_dt; CC(0)]
    end
end
evaluate_dt(sys::HCSystem, x, t) = evaluate_dt_augmented(sys, x, t)
