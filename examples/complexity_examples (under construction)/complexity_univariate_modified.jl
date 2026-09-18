using Printf

# Exact scalar experiment for
#   F_t(x) = x^2 - (1 + (m - 1)t),  x(t) = sqrt(1 + (m - 1)t).
#
# The complex box is B = [-1,1] + i[-1,1], equipped with the max norm on
# real and imaginary parts.  For this norm, ||uv|| <= 2||u||||v||, so
#
#   gamma(F_t, x) = ||JF_t(x)^(-1) J^2F_t(x) / 2!|| = 1 / |x|.
#
# We cap the algorithmic radius at one to match the definition
# r_alpha(t) = sup_{0 < r <= 1}{...} used by the complexity theorem.

const RHO = BigFloat(1) / 8
const TAU = BigFloat(7) / 8
const U_RHO = BigFloat("0.05")
const ALPHA = BigFloat("0.5")
const LAMBDA = BigFloat(154)

solution_value(m, t) = sqrt(1 + (m - 1) * t)
gamma_exact(x) = inv(abs(x))
tracking_radius(x) = min(BigFloat(1), U_RHO / (2 * gamma_exact(x)))

function intrinsic_radius(x)
    return min(BigFloat(1), ALPHA * x)
end

function intrinsic_length(m)
    endpoint = sqrt(BigFloat(m))
    breakpoint = inv(ALPHA)
    if endpoint <= breakpoint
        return log(endpoint) / ALPHA
    end
    return log(breakpoint) / ALPHA + endpoint - breakpoint
end

# Rigorous analytic upper bound for the normalized Krawczyk operator.
# For Y = 1/(2y) and radius R,
#
#   ||K(F_t,y,R,Y)|| / R
#       <= |F_t(y)|/(2|y|R) + 2R/|y|.
#
# The factor two in the second term comes from multiplying two complex
# square boxes in the real/imaginary max norm.
function normalized_krawczyk_bound(m, t, y, radius)
    a = 1 + (m - 1) * t
    beta_over_radius = abs(y^2 - a) / (2 * abs(y) * radius)
    jacobian_term = 2 * radius / abs(y)
    return beta_over_radius + jacobian_term
end

function refine_for_modified_track(m, t, y; max_newton=100)
    for newton_steps in 0:max_newton
        radius = tracking_radius(y)
        krawczyk_r = normalized_krawczyk_bound(m, t, y, radius)
        krawczyk_2r = normalized_krawczyk_bound(m, t, y, 2 * radius)
        if krawczyk_r <= RHO && krawczyk_2r <= RHO
            return y, radius, newton_steps, krawczyk_r, krawczyk_2r
        end

        a = 1 + (m - 1) * t
        y = (y + a / y) / 2
    end
    error("Newton refinement did not certify m=$m at t=$t")
end

function modified_apriori_track(m::Integer; max_steps=1_000_000)
    m > 1 || throw(ArgumentError("m must be greater than one"))

    m_big = BigFloat(m)
    t = BigFloat(0)
    y = BigFloat(1)
    steps = 0
    total_newton_steps = 0
    minimum_dt = BigFloat(Inf)
    maximum_dt = BigFloat(0)
    maximum_krawczyk = BigFloat(0)

    while t < 1
        steps < max_steps || error("maximum step count reached for m=$m")
        y, radius, newton_steps, kr, k2r = refine_for_modified_track(m_big, t, y)
        total_newton_steps += newton_steps
        maximum_krawczyk = max(maximum_krawczyk, kr, k2r)

        # JF_t(y)^(-1)F^(1)(y) = -(m-1)/(2y), and JF^(1)=0.
        speed_bound = (m_big - 1) / (2 * abs(y))
        dt = (TAU - RHO) * radius / speed_bound
        dt = min(dt, 1 - t)

        minimum_dt = min(minimum_dt, dt)
        maximum_dt = max(maximum_dt, dt)
        t += dt
        steps += 1
    end

    # Endpoint refinement is not counted as a tracking step.
    y, _, endpoint_newton_steps, kr, k2r = refine_for_modified_track(m_big, BigFloat(1), y)
    total_newton_steps += endpoint_newton_steps
    maximum_krawczyk = max(maximum_krawczyk, kr, k2r)

    exact_endpoint = sqrt(m_big)
    length = intrinsic_length(m_big)
    return (
        m=m,
        iterations=steps,
        length=length,
        ratio=BigFloat(steps) / length,
        newton_steps=total_newton_steps,
        endpoint_error=abs(y - exact_endpoint),
        minimum_dt=minimum_dt,
        maximum_dt=maximum_dt,
        maximum_krawczyk=maximum_krawczyk,
    )
end

function check_constants()
    u_limit = 1 - inv(sqrt(1 + RHO))
    lambda_limit = 2 * sqrt(BigFloat(2)) * exp(BigFloat(1)) / U_RHO
    lower_alpha_limit = 2 * RHO / (1 - RHO)
    upper_margin = 1 - RHO - ALPHA * (1 + RHO)

    @assert 0 < U_RHO < u_limit
    @assert LAMBDA > lambda_limit
    @assert lower_alpha_limit < ALPHA < 1
    @assert upper_margin > 0

    return (
        u_limit=u_limit,
        lambda_limit=lambda_limit,
        lower_alpha_limit=lower_alpha_limit,
        upper_margin=upper_margin,
    )
end

function write_results(path, results)
    open(path, "w") do io
        println(io, "m,iterations,L,iterations_per_L,newton_steps,endpoint_error,min_dt,max_dt,max_krawczyk")
        for result in results
            @printf(
                io,
                "%d,%d,%.12g,%.12g,%d,%.6e,%.6e,%.6e,%.6e\n",
                result.m,
                result.iterations,
                Float64(result.length),
                Float64(result.ratio),
                result.newton_steps,
                Float64(result.endpoint_error),
                Float64(result.minimum_dt),
                Float64(result.maximum_dt),
                Float64(result.maximum_krawczyk),
            )
        end
    end
end

function main()
    setprecision(BigFloat, 256) do
        constants = check_constants()
        @printf("rho=%.6f, tau=%.6f, u_rho=%.6f, alpha=%.6f\n",
            Float64(RHO), Float64(TAU), Float64(U_RHO), Float64(ALPHA))
        @printf("u_rho upper limit: %.9f\n", Float64(constants.u_limit))
        @printf("Lambda lower limit: %.9f; chosen Lambda: %.1f\n",
            Float64(constants.lambda_limit), Float64(LAMBDA))
        @printf("alpha lower limit: %.9f; upper-condition margin: %.9f\n\n",
            Float64(constants.lower_alpha_limit), Float64(constants.upper_margin))

        m_values = [10, 100, 1000, 10000, 20000, 30000]
        results = [modified_apriori_track(m) for m in m_values]
        @assert all(result.maximum_krawczyk <= RHO for result in results)

        println("       m    iterations              L          iters/L    Newton")
        for result in results
            @printf("%8d %13d %14.6f %16.6f %9d\n",
                result.m,
                result.iterations,
                Float64(result.length),
                Float64(result.ratio),
                result.newton_steps)
        end

        output_path = joinpath(@__DIR__, "results_modified_univariate.csv")
        write_results(output_path, results)
        println("\nWrote $output_path")
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
