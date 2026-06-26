using CertifiedHomotopyTracking
using Statistics

const AA = AbstractAlgebra

function _as_column_system(H)
    return Matrix(transpose(hcat(H)))
end

function _system_at(H, t)
    return evaluate_matrix(_as_column_system(H), t)
end

function _constant_dt_bound(H1, x, r, A)
    CC = parent(x[1])
    evalH1 = evaluate_matrix(H1, x)
    rB = CC("$r +/- $r")
    r_point = map(z -> z + rB, x)
    evalJH1 = evaluate_matrix(jac(H1), r_point)

    norm_evalH1 = max_norm(A * transpose(evalH1))
    norm_evalJH1 = max_norm(A * transpose(evalJH1))
    return (3 / 4) * r / (norm_evalH1 + norm_evalJH1 * r)
end

function _path_geometry_diagnostics(Ft, H1, x, r, A, rho)
    CC = parent(x[1])
    rB = CC("$r +/- $r")
    r_point = map(z -> z + rB, x)

    beta = max_norm(A * transpose(evaluate_matrix(Ft, x)))
    norm_H1 = max_norm(A * transpose(evaluate_matrix(H1, x)))
    norm_JH1 = max_norm(A * transpose(evaluate_matrix(jac(H1), r_point)))

    r_theory = beta / (1 + rho)
    radius_gap = r_theory == 0 ? Inf : r / r_theory
    eta = norm_H1 == 0 ? Inf : r * norm_JH1 / norm_H1

    return (
        beta=beta,
        r_theory=r_theory,
        radius_gap=radius_gap,
        eta=eta,
        norm_H1=norm_H1,
        norm_JH1=norm_JH1,
    )
end

function _jacobian_inverse_from_jac(J, vec)
    Jv = evaluate_matrix(J, vec)
    return inv(Jv)
end

function _krawczyk_test_cached_jac(system, J, point, r, A, rho)
    n = length(system)
    CC = base_ring(system[1])
    B = CC("+/- 1", "+/-1")
    I = identity_matrix(CC, n)
    box = fill(B, n, 1)

    eval_sys = evaluate_matrix(system, point)
    eval_jac = evaluate_matrix(J, vec(point + r * box))
    K = (-1 / r) * A * transpose(eval_sys) + (I - A * eval_jac) * matrix(box)
    return max_norm(K) < rho
end

function _refine_moore_box_cached_jac(system, J, point, r, A, rho)
    y = point
    while _krawczyk_test_cached_jac(system, J, y, r, A, rho) == false
        d = A * transpose(evaluate_matrix(system, y))
        if max_norm(d) <= (1 / 64) * rho * r
            r = (1 / 2) * r
        else
            y = midpoint_complex_box(y - d[:, 1])
        end
        A = _jacobian_inverse_from_jac(J, y)
    end

    while 2 * r <= 1 && _krawczyk_test_cached_jac(system, J, point, 2 * r, A, rho)
        r = 2 * r
    end

    return y, r, A
end

function _finite_mean(values)
    finite_values = filter(isfinite, values)
    return isempty(finite_values) ? Inf : mean(finite_values)
end

function _reaches_float_endpoint(t, dt)
    return t + dt >= 1 || t + dt == t || 1 - t <= 8 * eps(Float64)
end

function _acb_column(v)
    CC = parent(v[1])
    col = zero_matrix(CC, length(v), 1)
    for i in eachindex(v)
        col[i, 1] = v[i]
    end
    return col
end

function _filled_x_values(R, xvals, S)
    vals = Vector{elem_type(typeof(S))}(undef, AA.ngens(R))
    for i in 1:AA.ngens(R)
        vals[i] = i <= length(xvals) ? xvals[i] : zero(S)
    end
    return vals
end

function _compose_homotopy_polynomial(h, tpoly, xvals)
    S = parent(tpoly)
    result = zero(S)
    for k in 0:AA.degree(h)
        result += AA.evaluate(AA.coeff(h, k), xvals) * tpoly^k
    end
    return result
end

function _compose_homotopy(H, tpoly, xvals)
    return [_compose_homotopy_polynomial(h, tpoly, xvals) for h in H]
end

function _jac_homotopy(H, n)
    HR = parent(H[1])
    R = base_ring(HR)
    t = AA.gen(HR)
    JH = Matrix{typeof(H[1])}(undef, n, n)
    for row in 1:n, col in 1:n
        entry = zero(HR)
        for k in 0:AA.degree(H[row])
            entry += HR(derivative(AA.coeff(H[row], k), AA.gen(R, col))) * t^k
        end
        JH[row, col] = entry
    end
    return JH
end

function _jet_zero(CC, q)
    return [zero(CC) for _ in 0:q]
end

function _jet_const(c, CC, q)
    out = _jet_zero(CC, q)
    out[1] = CC(c)
    return out
end

function _jet_add(a, b, q)
    return [a[i] + b[i] for i in 1:q + 1]
end

function _jet_mul(a, b, q)
    CC = parent(a[1])
    out = _jet_zero(CC, q)
    for i in 0:q
        ai = a[i + 1]
        iszero(ai) && continue
        for j in 0:q - i
            bj = b[j + 1]
            iszero(bj) && continue
            out[i + j + 1] += ai * bj
        end
    end
    return out
end

function _jet_scale(c, a, q)
    return [c * a[i] for i in 1:q + 1]
end

function _jet_pow(a, exponent::Int, q)
    CC = parent(a[1])
    result = _jet_const(one(CC), CC, q)
    exponent == 0 && return result
    base = a
    e = exponent
    while e > 0
        if isodd(e)
            result = _jet_mul(result, base, q)
        end
        e ÷= 2
        e > 0 && (base = _jet_mul(base, base, q))
    end
    return result
end

function _eval_mpoly_jet(p, vals, q, CC)
    total = _jet_zero(CC, q)
    for term in AA.terms(p)
        c = AA.coeff(term, 1)
        exps = AA.exponent_vector(term, 1)
        nz = findall(!iszero, exps)
        if isempty(nz)
            term_jet = _jet_const(c, CC, q)
        elseif length(nz) == 1 && exps[nz[1]] == 1
            term_jet = _jet_scale(c, vals[nz[1]], q)
        elseif length(nz) == 1 && exps[nz[1]] == 2
            term_jet = _jet_scale(c, _jet_mul(vals[nz[1]], vals[nz[1]], q), q)
        elseif length(nz) == 2 && exps[nz[1]] == 1 && exps[nz[2]] == 1
            term_jet = _jet_scale(c, _jet_mul(vals[nz[1]], vals[nz[2]], q), q)
        else
            term_jet = _jet_const(c, CC, q)
            for var_idx in nz
                term_jet = _jet_mul(term_jet, _jet_pow(vals[var_idx], exps[var_idx], q), q)
            end
        end
        total = _jet_add(total, term_jet, q)
    end
    return total
end

function _compose_homotopy_jet(H, tjet, xjets, q, CC)
    out = Vector{Vector{AcbFieldElem}}(undef, length(H))
    for i in eachindex(H)
        hi = H[i]
        total = _jet_zero(CC, q)
        for k in 0:AA.degree(hi)
            coeff_jet = _eval_mpoly_jet(AA.coeff(hi, k), xjets, q, CC)
            total = _jet_add(total, _jet_mul(coeff_jet, _jet_pow(tjet, k, q), q), q)
        end
        out[i] = total
    end
    return out
end

function _jet_zero_over(R, q)
    return [zero(R) for _ in 0:q]
end

function _jet_const_over(c, R, q)
    out = _jet_zero_over(R, q)
    out[1] = R(c)
    return out
end

function _jet_add_generic(a, b, R, q)
    return [a[i] + b[i] for i in 1:q + 1]
end

function _jet_mul_generic(a, b, R, q)
    out = _jet_zero_over(R, q)
    for i in 0:q
        ai = a[i + 1]
        iszero(ai) && continue
        for j in 0:q - i
            bj = b[j + 1]
            iszero(bj) && continue
            out[i + j + 1] += ai * bj
        end
    end
    return out
end

function _jet_pow_generic(a, exponent::Int, R, q)
    result = _jet_const_over(one(R), R, q)
    exponent == 0 && return result
    base = a
    e = exponent
    while e > 0
        if isodd(e)
            result = _jet_mul_generic(result, base, R, q)
        end
        e ÷= 2
        e > 0 && (base = _jet_mul_generic(base, base, R, q))
    end
    return result
end

function _eval_mpoly_jet_generic(p, vals, R, q)
    total = _jet_zero_over(R, q)
    for term in AA.terms(p)
        term_jet = _jet_const_over(AA.coeff(term, 1), R, q)
        exps = AA.exponent_vector(term, 1)
        for (var_idx, exponent) in enumerate(exps)
            exponent == 0 && continue
            term_jet = _jet_mul_generic(term_jet, _jet_pow_generic(vals[var_idx], exponent, R, q), R, q)
        end
        total = _jet_add_generic(total, term_jet, R, q)
    end
    return total
end

function _compose_homotopy_polynomial_jet_generic(h, tjet, xjets, R, q)
    total = _jet_zero_over(R, q)
    for k in 0:AA.degree(h)
        coeff_jet = _eval_mpoly_jet_generic(AA.coeff(h, k), xjets, R, q)
        total = _jet_add_generic(total, _jet_mul_generic(coeff_jet, _jet_pow_generic(tjet, k, R, q), R, q), R, q)
    end
    return total
end

function _predictor_coefficients(H, t0, x, A, q)
    CC = parent(x[1])
    R = base_ring(parent(H[1]))
    coeffs = [zero(CC) for _ in 1:length(x), _ in 1:q]

    for m in 1:q
        xseries = [
            begin
                jet = _jet_zero(CC, q)
                jet[1] = x[j]
                for ell in 1:m - 1
                    jet[ell + 1] = coeffs[j, ell]
                end
                jet
            end
            for j in 1:length(x)
        ]
        vals = Vector{Vector{AcbFieldElem}}(undef, AA.ngens(R))
        for i in 1:AA.ngens(R)
            vals[i] = i <= length(xseries) ? xseries[i] : _jet_zero(CC, q)
        end
        tjet = _jet_zero(CC, q)
        tjet[1] = CC(t0)
        tjet[2] = one(CC)
        Fseries = _compose_homotopy_jet(H, tjet, vals, q, CC)
        b = [Fseries[i][m + 1] for i in 1:length(x)]
        delta = A * _acb_column(b)
        for j in 1:length(x)
            coeffs[j, m] = -delta[j, 1]
        end
    end

    return coeffs
end

function _predictor_value(x, coeffs, h)
    CC = parent(x[1])
    q = size(coeffs, 2)
    hh = CC(h)
    return [
        x[j] + sum(coeffs[j, m] * hh^m for m in 1:q; init=zero(CC))
        for j in 1:length(x)
    ]
end

function _homotopy_s_degree_bound(H, q)
    maxdeg = 0
    for h in H
        for k in 0:AA.degree(h)
            coeff = AA.coeff(h, k)
            for term in AA.terms(coeff)
                degree_x = sum(AA.exponent_vector(term, 1))
                maxdeg = max(maxdeg, k + q * degree_x)
            end
        end
    end
    return maxdeg
end

function _c_norms(H, t0, x, A, coeffs)
    CC = parent(x[1])
    R = base_ring(parent(H[1]))
    q = size(coeffs, 2)
    maxdeg = _homotopy_s_degree_bound(H, q)
    norms = zeros(Float64, max(0, maxdeg))
    maxdeg <= q && return norms

    xjets = [
        begin
            jet = _jet_zero(CC, maxdeg)
            jet[1] = x[j]
            for m in 1:q
                jet[m + 1] = coeffs[j, m]
            end
            jet
        end
        for j in 1:length(x)
    ]
    vals = Vector{Vector{AcbFieldElem}}(undef, AA.ngens(R))
    for i in 1:AA.ngens(R)
        vals[i] = i <= length(xjets) ? xjets[i] : _jet_zero(CC, maxdeg)
    end
    tjet = _jet_zero(CC, maxdeg)
    tjet[1] = CC(t0)
    tjet[2] = one(CC)
    Fseries = _compose_homotopy_jet(H, tjet, vals, maxdeg, CC)

    for k in q + 1:maxdeg
        ck = A * _acb_column([Fseries[i][k + 1] for i in 1:length(x)])
        norms[k] = max_norm(ck)
    end
    return norms
end

function _unit_complex_box(CC)
    return CC("0 +/- 1") + onei(CC) * CC("0 +/- 1")
end

function _evaluate_u_box(p, ubox)
    return length(ubox) == 1 ? AA.evaluate(p, ubox[1]) : AA.evaluate(p, ubox)
end

function _M_norms(H, t0, x, r, A, coeffs)
    CC = parent(x[1])
    n = length(x)
    if n == 1
        U, (u) = CC["u_1"]
        uvars = [u]
    else
        U, _ = CC[["u_$i" for i in 1:n]...]
        uvars = [AA.gen(U, i) for i in 1:n]
    end
    R = base_ring(parent(H[1]))
    q = size(coeffs, 2)
    s_order = q + 1

    xq = [
        begin
            jet = _jet_zero_over(U, s_order)
            jet[1] = U(x[j])
            for m in 1:q
                jet[m + 1] = U(coeffs[j, m])
            end
            jet
        end
        for j in 1:n
    ]
    moving_vals = [
        begin
            jet = copy(xq[j])
            jet[1] += U(CC(r)) * uvars[j]
            jet
        end
        for j in 1:n
    ]
    base_vals = [
        begin
            jet = _jet_zero_over(U, s_order)
            jet[1] = U(x[j]) + U(CC(r)) * uvars[j]
            jet
        end
        for j in 1:n
    ]
    while length(moving_vals) < AA.ngens(R)
        push!(moving_vals, _jet_zero_over(U, s_order))
        push!(base_vals, _jet_zero_over(U, s_order))
    end

    moving_t = _jet_zero_over(U, s_order)
    moving_t[1] = U(CC(t0))
    moving_t[2] = one(U)
    base_t = _jet_zero_over(U, s_order)
    base_t[1] = U(CC(t0))

    JH = _jac_homotopy(H, n)
    moving = Matrix{Vector{elem_type(typeof(U))}}(undef, n, n)
    base = Matrix{Vector{elem_type(typeof(U))}}(undef, n, n)
    for row in 1:n, col in 1:n
        moving[row, col] = _compose_homotopy_polynomial_jet_generic(
            JH[row, col],
            moving_t,
            moving_vals,
            U,
            s_order,
        )
        base[row, col] = _compose_homotopy_polynomial_jet_generic(
            JH[row, col],
            base_t,
            base_vals,
            U,
            s_order,
        )
    end

    maxdeg = s_order
    norms = zeros(Float64, max(0, maxdeg))
    ubox = [_unit_complex_box(CC) for _ in 1:n]

    for k in 1:maxdeg
        interval_matrix = Matrix{AcbFieldElem}(undef, n, n)
        for row in 1:n, col in 1:n
            entry = zero(CC)
            for ell in 1:n
                delta = _evaluate_u_box(moving[ell, col][k + 1], ubox) -
                    _evaluate_u_box(base[ell, col][k + 1], ubox)
                entry += A[row, ell] * delta
            end
            interval_matrix[row, col] = entry
        end
        norms[k] = max_norm(interval_matrix)
    end
    return norms
end

function _solve_positive_polynomial_bound(coeffs, rhs, upper)
    rhs <= 0 && return 0.0
    all(iszero, coeffs) && return upper

    value_at(h) = sum(coeffs[k] * h^k for k in eachindex(coeffs))
    hi = upper
    lo = 0.0
    if value_at(hi) <= rhs
        return hi
    end
    for _ in 1:80
        mid = (lo + hi) / 2
        if value_at(mid) <= rhs
            lo = mid
        else
            hi = mid
        end
    end
    return lo
end

function higher_order_dt_bound(H, t0, x, r, A; order=2, rho=1 / 8, tau=7 / 8, max_step=1 - t0)
    coeffs = _predictor_coefficients(H, t0, x, A, order)
    c_norms = _c_norms(H, t0, x, A, coeffs)
    M_norms = _M_norms(H, t0, x, r, A, coeffs)
    degree_bound = max(length(c_norms), length(M_norms))
    bound_coeffs = zeros(Float64, degree_bound)

    for k in 1:degree_bound
        bound_coeffs[k] =
            (k <= length(c_norms) ? c_norms[k] : 0.0) +
            r * (k <= length(M_norms) ? M_norms[k] : 0.0)
    end

    dt = _solve_positive_polynomial_bound(bound_coeffs, r * (tau - rho), max_step)
    return dt, coeffs, c_norms, M_norms, bound_coeffs
end

function tracking_higher_order_apriori(
    H,
    x,
    r;
    order=2,
    rho=1 / 8,
    tau=7 / 8,
    show_display=true,
    iterations_count=false,
    final_refine=true,
    max_iterations=typemax(Int),
    compare_constant=true,
)
    t = 0.0
    iter = 0
    dt_list = Float64[]
    r_list = Float64[]
    constant_dt_list = Float64[]
    beta_list = Float64[]
    r_theory_list = Float64[]
    radius_gap_list = Float64[]
    eta_list = Float64[]
    H1 = evaluate_matrix(_as_column_system(map(derivative, H)), 0)
    status = :success

    G = _system_at(H, 0)
    A = jacobian_inverse(G, x)

    while t < 1
        if iter >= max_iterations
            status = :max_iterations
            break
        end
        Ft = _system_at(H, t)
        JFt = jac(Ft)
        x, r, A = _refine_moore_box_cached_jac(Ft, JFt, x, r, A, rho)
        diag = _path_geometry_diagnostics(Ft, H1, x, r, A, rho)

        constant_dt = compare_constant ? min(1 - t, _constant_dt_bound(H1, x, r, A)) : NaN
        dt, coeffs = higher_order_dt_bound(H, t, x, r, A; order, rho, tau, max_step=1 - t)
        dt = min(dt, 1 - t)
        dt <= 0 && error("Computed a non-positive stepsize at t=$t.")

        push!(constant_dt_list, constant_dt)
        push!(dt_list, dt)
        push!(r_list, r)
        push!(beta_list, diag.beta)
        push!(r_theory_list, diag.r_theory)
        push!(radius_gap_list, diag.radius_gap)
        push!(eta_list, diag.eta)

        reaches_endpoint = _reaches_float_endpoint(t, dt)
        show_display && print("\r(order=$order, dt=$(round(dt; sigdigits=5)), constant=$(round(constant_dt; sigdigits=5)), t=$(round(t; sigdigits=8)), r=$r)")
        x = _predictor_value(x, coeffs, dt)
        t = reaches_endpoint ? 1.0 : t + dt
        iter += 1
    end

    if final_refine && status == :success
        Ft = _system_at(H, 1)
        JFt = jac(Ft)
        x, r, A = _refine_moore_box_cached_jac(Ft, JFt, x, r, A, rho)
    end
    show_display && println()

    if iterations_count
        return x, (
            status=status,
            iterations=iter,
            final_t=t,
            min_dt=minimum(dt_list),
            median_dt=median(dt_list),
            max_dt=maximum(dt_list),
            mean_dt=mean(dt_list),
            mean_constant_dt=mean(constant_dt_list),
            min_ratio=minimum(dt_list ./ constant_dt_list),
            median_ratio=median(dt_list ./ constant_dt_list),
            mean_ratio=mean(dt_list ./ constant_dt_list),
            min_radius=minimum(r_list),
            median_radius=median(r_list),
            mean_radius_theory_gap=_finite_mean(radius_gap_list),
            max_eta=maximum(eta_list),
            mean_beta=mean(beta_list),
            min_r_theory=minimum(r_theory_list),
            dt_list=dt_list,
            radius_list=r_list,
            r_theory_list=r_theory_list,
            radius_gap_list=radius_gap_list,
            eta_list=eta_list,
        )
    end
    return x
end

function tracking_constant_apriori(
    H,
    x,
    r;
    rho=1 / 8,
    show_display=true,
    iterations_count=false,
    final_refine=true,
)
    t = 0.0
    iter = 0
    dt_list = Float64[]
    r_list = Float64[]
    beta_list = Float64[]
    r_theory_list = Float64[]
    radius_gap_list = Float64[]
    eta_list = Float64[]
    H1 = evaluate_matrix(_as_column_system(map(derivative, H)), 0)

    G = _system_at(H, 0)
    A = jacobian_inverse(G, x)

    while t < 1
        Ft = _system_at(H, t)
        JFt = jac(Ft)
        x, r, A = _refine_moore_box_cached_jac(Ft, JFt, x, r, A, rho)
        diag = _path_geometry_diagnostics(Ft, H1, x, r, A, rho)

        dt = min(1 - t, _constant_dt_bound(H1, x, r, A))
        dt <= 0 && error("Computed a non-positive constant stepsize at t=$t.")

        push!(dt_list, dt)
        push!(r_list, r)
        push!(beta_list, diag.beta)
        push!(r_theory_list, diag.r_theory)
        push!(radius_gap_list, diag.radius_gap)
        push!(eta_list, diag.eta)

        reaches_endpoint = _reaches_float_endpoint(t, dt)
        show_display && print("\r(constant, dt=$(round(dt; sigdigits=5)), t=$(round(t; sigdigits=8)), r=$r)")
        t = reaches_endpoint ? 1.0 : t + dt
        iter += 1
    end

    if final_refine
        Ft = _system_at(H, 1)
        JFt = jac(Ft)
        x, r, A = _refine_moore_box_cached_jac(Ft, JFt, x, r, A, rho)
    end
    show_display && println()

    if iterations_count
        return x, (
            iterations=iter,
            min_dt=minimum(dt_list),
            median_dt=median(dt_list),
            max_dt=maximum(dt_list),
            mean_dt=mean(dt_list),
            min_radius=minimum(r_list),
            median_radius=median(r_list),
            mean_radius_theory_gap=_finite_mean(radius_gap_list),
            max_eta=maximum(eta_list),
            mean_beta=mean(beta_list),
            min_r_theory=minimum(r_theory_list),
            dt_list=dt_list,
            radius_list=r_list,
            r_theory_list=r_theory_list,
            radius_gap_list=radius_gap_list,
            eta_list=eta_list,
        )
    end
    return x
end

function demo_univariate(; m=10_000, order=2, r=0.1)
    CC = AcbField()
    R, (a, eta) = CC["a", "eta"]
    HR, (t) = R["t"]

    F = [a^2 - m]
    G = [a^2 - 1]
    H = (((1 - t) .* transpose(hcat(G))) + (t .* transpose(hcat(F))))[1, :]

    x, stats = tracking_higher_order_apriori(
        H,
        [CC(1.0)],
        r;
        order,
        iterations_count=true,
    )
    println("final x = ", x)
    println("stats = ", stats)
    return x, stats
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo_univariate(order=2)
end
