using CertifiedHomotopyTracking
import HomotopyContinuation as HC

# Alt's problem in the 8-variable formulation of Hauenstein and Helmer,
# "Probabilistic Saturations and Alt's Problem", Section 2.1.
#
# Execute this file section by section.  HC.find_start_pair constructs one
# generic parameter/start-point pair, and CHT then certifies that single path
# to the requested Alt system.

# -----------------------------------------------------------------------------
# 1. Eight design variables
# -----------------------------------------------------------------------------

@variables a abar b bbar x xbar y ybar
VARS = [a, abar, b, bbar, x, xbar, y, ybar]

PREC_BITS = 256
CC = AcbField(PREC_BITS)

# Algebraic conjugation in isotropic coordinates.  This swaps every variable
# with its barred partner; it does not conjugate numerical coefficients.
bar_swap = Dict(a => abar, abar => a, b => bbar, bbar => b,
                x => xbar, xbar => x, y => ybar, ybar => y)
algconj(f) = Symbolics.substitute(f, bar_swap)

# -----------------------------------------------------------------------------
# 2. Coupler-curve basis f_1,...,f_15
# -----------------------------------------------------------------------------

# These expressions are transcribed from the authors' public Magma code.  The
# intermediate names m1,...,m15 follow that code; PAPER_BASIS reorders them to
# the f_1,...,f_15 order used in the paper.
m1 = (xbar - ybar) * (x - y)
m2 = (x - y) * (abar*xbar - 2*abar*ybar + 2*bbar*xbar - bbar*ybar)
m3 = (abar^2*ybar - 2*abar*bbar*xbar + 2*abar*bbar*ybar - bbar^2*xbar) * (x - y)
m4 = abar*bbar*(x - y)*(abar*ybar - bbar*xbar)

m6 = -a*abar*x*xbar + a*abar*x*ybar + a*abar*xbar*y - 2*a*abar*y*ybar -
     2*a*bbar*x*xbar + a*bbar*x*ybar + 4*a*bbar*xbar*y - 2*a*bbar*y*ybar +
     a*x*xbar*ybar + a*xbar^2*y - 2*a*xbar*y*ybar - 2*abar*b*x*xbar +
     4*abar*b*x*ybar + abar*b*xbar*y - 2*abar*b*y*ybar + abar*x^2*ybar +
     abar*x*xbar*y - 2*abar*x*y*ybar - 2*b*bbar*x*xbar + b*bbar*x*ybar +
     b*bbar*xbar*y - b*bbar*y*ybar - 2*b*x*xbar*ybar + b*x*ybar^2 +
     b*xbar*y*ybar - 2*bbar*x*xbar*y + bbar*x*y*ybar + bbar*xbar*y^2 -
     x^2*ybar^2 + 2*x*xbar*y*ybar - xbar^2*y^2

m7 = 2*a*abar*bbar*x*xbar - a*abar*bbar*x*ybar - 2*a*abar*bbar*xbar*y +
     2*a*abar*bbar*y*ybar - a*abar*x*xbar*ybar + 2*a*abar*xbar*y*ybar +
     a*bbar^2*x*xbar - 2*a*bbar^2*xbar*y - a*bbar*x*xbar*ybar -
     2*a*bbar*xbar^2*y + 2*a*bbar*xbar*y*ybar - 2*abar^2*b*x*ybar +
     abar^2*b*y*ybar - abar^2*x^2*ybar + 2*abar^2*x*y*ybar +
     2*abar*b*bbar*x*xbar - 2*abar*b*bbar*x*ybar - abar*b*bbar*xbar*y +
     2*abar*b*bbar*y*ybar + 2*abar*b*x*xbar*ybar - 2*abar*b*x*ybar^2 -
     abar*b*xbar*y*ybar - abar*bbar*x^2*ybar - abar*bbar*xbar*y^2 +
     2*abar*x^2*ybar^2 - 2*abar*x*xbar*y*ybar + 2*b*bbar*x*xbar*ybar -
     b*bbar*xbar*y*ybar + 2*bbar^2*x*xbar*y - bbar^2*xbar*y^2 -
     2*bbar*x*xbar*y*ybar + 2*bbar*xbar^2*y^2

m8 = -a*abar*bbar^2*x*xbar + a*abar*bbar^2*xbar*y +
     a*abar*bbar*x*xbar*ybar - 2*a*abar*bbar*xbar*y*ybar +
     a*bbar^2*xbar^2*y + abar^2*b*bbar*x*ybar - abar^2*b*bbar*y*ybar +
     abar^2*b*x*ybar^2 + abar^2*bbar*x^2*ybar - abar^2*bbar*x*y*ybar -
     abar^2*x^2*ybar^2 - 2*abar*b*bbar*x*xbar*ybar +
     abar*b*bbar*xbar*y*ybar - abar*bbar^2*x*xbar*y +
     abar*bbar^2*xbar*y^2 + 2*abar*bbar*x*xbar*y*ybar - bbar^2*xbar^2*y^2

m11 = a^2*bbar^2*xbar*y + 2*a^2*bbar*xbar^2*y -
      2*a^2*bbar*xbar*y*ybar - a^2*xbar^2*y*ybar -
      2*a*abar*b*bbar*x*xbar + a*abar*b*bbar*x*ybar +
      a*abar*b*bbar*xbar*y - 2*a*abar*b*bbar*y*ybar +
      a*abar*b*x*ybar^2 - a*abar*b*xbar*y*ybar -
      a*abar*bbar*x*y*ybar + a*abar*bbar*xbar*y^2 -
      a*b*bbar*x*xbar*ybar + a*b*bbar*xbar^2*y +
      a*b*x*xbar*ybar^2 + a*b*xbar^2*y*ybar -
      2*a*bbar^2*x*xbar*y + 2*a*bbar^2*xbar*y^2 +
      3*a*bbar*x*xbar*y*ybar - 3*a*bbar*xbar^2*y^2 +
      abar^2*b^2*x*ybar + 2*abar^2*b*x^2*ybar -
      2*abar^2*b*x*y*ybar - abar^2*x^2*y*ybar -
      2*abar*b^2*x*xbar*ybar + 2*abar*b^2*x*ybar^2 +
      abar*b*bbar*x^2*ybar - abar*b*bbar*x*xbar*y -
      3*abar*b*x^2*ybar^2 + 3*abar*b*x*xbar*y*ybar +
      abar*bbar*x^2*y*ybar + abar*bbar*x*xbar*y^2 -
      b^2*x*xbar*ybar^2 - bbar^2*x*xbar*y^2

m12 = (a*bbar*xbar*y - abar*b*x*ybar) *
      (a*bbar*xbar - a*xbar*ybar - abar*b*ybar - abar*bbar*x +
       abar*bbar*y + abar*x*ybar + b*xbar*ybar - bbar*xbar*y)

m5 = algconj(m2)
m9 = algconj(m3)
m13 = algconj(m4)
m10 = algconj(m7)
m14 = algconj(m8)
m15 = algconj(m12)

# Paper order: f2/f3, f4/f5, f6/f7, f9/f10, f11/f12, and f14/f15
# are algebraic-conjugate pairs.
PAPER_BASIS = [m1, m2, m5, m3, m9, m4, m13, m6,
               m7, m10, m8, m14, m11, m12, m15]

# -----------------------------------------------------------------------------
# 3. The coefficient-parametrized 8-equation system
# -----------------------------------------------------------------------------

@variables coeff[1:120]
COEFFS = collect(coeff)

G_param = [
    sum(COEFFS[(i - 1)*15 + j] * PAPER_BASIS[j] for j in 1:15)
    for i in 1:8
]

# The coefficient vector c(p,pbar) from equation (Section 2.1) of the paper.
function coupler_coefficients(p, pbar=conj(p))
    return ComplexF64[
        p^3*pbar^3,
        p^3*pbar^2,
        p^2*pbar^3,
        p^3*pbar,
        p*pbar^3,
        p^3,
        pbar^3,
        p^2*pbar^2,
        p^2*pbar,
        p*pbar^2,
        p^2,
        pbar^2,
        p*pbar,
        p,
        pbar,
    ]
end

# -----------------------------------------------------------------------------
# 4. Target Alt instance: nine prescribed points
# -----------------------------------------------------------------------------

coupler_points = ComplexF64[
    0.8961867 - 0.09802917im,
    1.2156535 - 1.18749100im,
    1.5151435 - 0.85449808im,
    1.6754775 - 0.48768058im,
    1.7138690 - 0.30099232im,
    1.7215236 + 0.03269953im,
    1.6642029 + 0.33241088im,
    1.4984171 + 0.74435576im,
    1.3011834 + 0.92153806im,
]

# The paper places the first point at the origin.
precision_points = coupler_points[2:9] .- coupler_points[1]
target_coefficients = reduce(vcat, coupler_coefficients.(precision_points))

# -----------------------------------------------------------------------------
# 5. Compile the coefficient-parametrized system
# -----------------------------------------------------------------------------

compiled = compile_edge_homotopy(G_param, VARS, COEFFS)

# This temporary specialization is used only to expose the corresponding
# parameterized HC.jl System.  Its endpoints are not tracked or certified.
H_template = make_edge_system(
    compiled,
    CC.(target_coefficients),
    CC.(target_coefficients),
)

# Reuse CHT's source metadata to obtain the exactly corresponding HC.jl system.
posteriori_tracker = prepare_posteriori_tracker(H_template)
hc_system = posteriori_hc_system(posteriori_tracker)

# -----------------------------------------------------------------------------
# 6. Ask HomotopyContinuation.jl for one generic start pair
# -----------------------------------------------------------------------------

# `p` is a generic 120-entry coefficient vector and `start_point` satisfies
# hc_system(start_point; parameters=p) approximately.
start_point, p = HC.find_start_pair(
    hc_system;
    max_tries=100,
    atol=0.0,
    rtol=1e-12,
)

# Build the actual CHT homotopy from the discovered generic parameter `p` to
# the coefficients determined by the eight prescribed precision points.
H = make_edge_system(
    compiled,
    CC.(p),
    CC.(target_coefficients),
)

start_point = CC.(start_point)

# -----------------------------------------------------------------------------
# 7. A posteriori certification with CHT
# -----------------------------------------------------------------------------

# Kept only as a reference; do not execute for the intended workflow.
# res = track_path(H, start_point; show_progress=true)

# Execute the path from the start pair found above.
cert = certify_posteriori(H, start_point; show_progress=true, visualize=(axes=(:t, (1, :real), (2, :imag)), show_trace=true))

cert.success
cert.hc_trace.status
cert.total_boxes
cert.max_depth
cert.failed_segments
export_path_tikz(
                  cert,
                  "posteriori_visualization.tex";
                  axes=(:t, (1, :real), (1, :imag)),
                  show_trace=true,
              )