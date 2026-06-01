export CompiledHomotopy, HCSystem, TMCache, evaluate_H, evaluate_Jac, evaluate_dt,
       system_with_precision

struct CompiledHomotopy
    func_H::Function
    func_Jx::Function
    func_dt::Function
    chart_func_H::Vector{Function}
    chart_func_Jx::Vector{Function}
    chart_func_dt::Vector{Function}
    homogeneous::Bool
    n_vars::Int
    projective_patch::Bool
    patch_template::Vector{ComplexF64}
end
CompiledHomotopy(func_H::Function, func_Jx::Function, func_dt::Function, homogeneous::Bool) =
    CompiledHomotopy(func_H, func_Jx, func_dt, Function[], Function[], Function[], homogeneous, 0, false, ComplexF64[])
CompiledHomotopy(func_H::Function, func_Jx::Function, func_dt::Function, homogeneous::Bool, n_vars::Int) =
    CompiledHomotopy(func_H, func_Jx, func_dt, Function[], Function[], Function[], homogeneous, n_vars, false, ComplexF64[])

function Base.show(io::IO, compiled::CompiledHomotopy)
    print(io, "CompiledHomotopy(Homogeneous: ", compiled.homogeneous, ", n_vars=", compiled.n_vars, ", projective_patch=", compiled.projective_patch, ")")
end

mutable struct HCSystem
    compiled::CompiledHomotopy
    p_start::Tuple{Vararg{AcbFieldElem}}
    p_end::Tuple{Vararg{AcbFieldElem}}
    p_const::Tuple{Vararg{AcbFieldElem}}
    homogeneous::Bool
    patch_idx::Int
    patch_vector::Tuple{Vararg{AcbFieldElem}}
    
    CC::AcbField 
    RR::ArbField
    
    function HCSystem(
        compiled::CompiledHomotopy,
        p_start::Vector{AcbFieldElem},
        p_end::Vector{AcbFieldElem},
        p_const::Vector{AcbFieldElem}=AcbFieldElem[];
        patch_vector::Vector{AcbFieldElem}=AcbFieldElem[],
    )
        CC = !isempty(p_start) ? parent(p_start[1]) : parent(p_const[1])
        RR = ArbField(precision(CC))
        if compiled.projective_patch && !isempty(patch_vector)
            throw(ArgumentError("patch_vector is no longer used for projective tracking; coordinate charts are used instead."))
        end
        patch_values = patch_vector
        patch_tuple = Tuple(patch_values)
        if !isempty(patch_tuple) && compiled.n_vars != 0 && length(patch_tuple) != compiled.n_vars
            throw(DimensionMismatch("patch_vector length must match the number of compiled variables."))
        end
        new(compiled, Tuple(p_start), Tuple(p_end), Tuple(p_const), compiled.homogeneous, 1, patch_tuple, CC, RR)
    end

    function HCSystem(compiled::CompiledHomotopy, CC::AcbField; patch_vector::Vector{AcbFieldElem}=AcbFieldElem[])
        RR = ArbField(precision(CC))
        if compiled.projective_patch && !isempty(patch_vector)
            throw(ArgumentError("patch_vector is no longer used for projective tracking; coordinate charts are used instead."))
        end
        patch_values = patch_vector
        patch_tuple = Tuple(patch_values)
        if !isempty(patch_tuple) && compiled.n_vars != 0 && length(patch_tuple) != compiled.n_vars
            throw(DimensionMismatch("patch_vector length must match the number of compiled variables."))
        end
        new(compiled, (), (), (), compiled.homogeneous, 1, patch_tuple, CC, RR)
    end
end

function Base.show(io::IO, sys::HCSystem)
    print(
        io,
        "HCSystem(",
        sys.compiled,
        ", params=", length(sys.p_start),
        ", consts=", length(sys.p_const),
        ", precision=", precision(sys.CC),
        ")",
    )
end

has_projective_patch(sys::HCSystem) = !isempty(sys.patch_vector)
uses_projective_charts(sys::HCSystem) = sys.compiled.projective_patch

_convert_acb_vector(CC::AcbField, values) = AcbFieldElem[CC(value) for value in values]

function system_with_precision(sys::HCSystem, precision_bits::Integer)
    precision_bits > 0 || throw(ArgumentError("precision_bits must be positive."))
    precision(sys.CC) == precision_bits && return sys

    CC = AcbField(precision_bits)
    patch_vector = _convert_acb_vector(CC, sys.patch_vector)
    rebuilt = if isempty(sys.p_start)
        HCSystem(sys.compiled, CC; patch_vector=patch_vector)
    else
        HCSystem(
            sys.compiled,
            _convert_acb_vector(CC, sys.p_start),
            _convert_acb_vector(CC, sys.p_end),
            _convert_acb_vector(CC, sys.p_const);
            patch_vector=patch_vector,
        )
    end
    rebuilt.patch_idx = sys.patch_idx
    return rebuilt
end

function _patch_value(sys::HCSystem, x)
    result = zero(sys.patch_vector[1] * x[1])
    for i in 1:length(sys.patch_vector)
        result += sys.patch_vector[i] * x[i]
    end
    return result - 1
end

_coerce_eval_values(sys::HCSystem, values, x::AbstractVector{<:AcbFieldElem}) = sys.CC.(values)
_coerce_eval_values(sys::HCSystem, values, x) = values

function _chart_input_for_evaluation(sys::HCSystem, x)
    !uses_projective_charts(sys) && return x, sys.patch_idx
    length(x) == sys.compiled.n_vars - 1 && return x, sys.patch_idx
    length(x) == sys.compiled.n_vars || throw(DimensionMismatch("Projective chart evaluation expects chart length n or homogeneous length n + 1."))

    chart_idx = sys.patch_idx
    if mag_complex(x[chart_idx]) <= 1e-30
        _, chart_idx = findmax([mag_complex(xi) for xi in x])
    end
    scale = x[chart_idx]
    result = Vector{eltype(x)}(undef, length(x) - 1)
    out_idx = 1
    for i in eachindex(x)
        if i != chart_idx
            result[out_idx] = x[i] / scale
            out_idx += 1
        end
    end
    return result, chart_idx
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
    if uses_projective_charts(sys)
        x_chart, chart_idx = _chart_input_for_evaluation(sys, x)
        return _coerce_eval_values(sys, sys.compiled.chart_func_H[chart_idx](x_chart, t, sys.p_start..., sys.p_end..., sys.p_const...), x_chart)
    end
    val_sys = sys.compiled.func_H(x, t, sys.p_start..., sys.p_end..., sys.p_const...)
    if !sys.homogeneous
        return val_sys
    elseif has_projective_patch(sys)
        return [val_sys; _patch_value(sys, x)]
    else
        val_patch = x[sys.patch_idx] - 1
        return [val_sys; val_patch]
    end
end
evaluate_H(sys::HCSystem, x, t) = evaluate_H_augmented(sys, x, t)

function evaluate_Jac(sys::HCSystem, x, t)
    CC = sys.CC 
    if uses_projective_charts(sys)
        x_chart, chart_idx = _chart_input_for_evaluation(sys, x)
        return _coerce_eval_values(sys, sys.compiled.chart_func_Jx[chart_idx](x_chart, t, sys.p_start..., sys.p_end..., sys.p_const...), x_chart)
    end
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
            if has_projective_patch(sys)
                J_aug[n_rows+1, j] = sys.patch_vector[j]
            else
                J_aug[n_rows+1, j] = (j == sys.patch_idx) ? CC(1) : CC(0)
            end
        end
        return J_aug
    end
end

function evaluate_dt_augmented(sys::HCSystem, x, t)
    CC = sys.CC 
    if uses_projective_charts(sys)
        x_chart, chart_idx = _chart_input_for_evaluation(sys, x)
        return _coerce_eval_values(sys, sys.compiled.chart_func_dt[chart_idx](x_chart, t, sys.p_start..., sys.p_end..., sys.p_const...), x_chart)
    end
    val_dt = sys.compiled.func_dt(x, t, sys.p_start..., sys.p_end..., sys.p_const...)
    if !sys.homogeneous
        return val_dt
    else
        return [val_dt; CC(0)]
    end
end
evaluate_dt(sys::HCSystem, x, t) = evaluate_dt_augmented(sys, x, t)
