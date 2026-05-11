#Pkg.instantiate()

using CertifiedHomotopyTracking

# ------------------------------------------------------------------------------
# 1. Setup Polynomial Ring & System
# ------------------------------------------------------------------------------
# @setupfield to define variables and fields
@variables x y z λ 
@variables u1 u2 u3
PREC_BITS = 256
CC = AcbField(PREC_BITS)


F = [2*(x-u1)-6*λ*(x^2+y^2)^2*x, 2*(y-u2)-6*λ*(x^2+y^2)^2*y, 2*(z-u3)+4*λ*z^3, 0*u1+z^4-(x^2+y^2)^3]
vars = [x, y, z, λ]
pars = [u1, u2, u3]

# ------------------------------------------------------------------------------
# 2. Setup Vertices & Initial Points
# ------------------------------------------------------------------------------
# define base_point 
bp = [CC(.09868,.675389), CC(.423238,.713082), CC(.592351,.144969)]
x0 = [CC(1.23836,-.422501), CC(1.19574,-1.0474), CC(2.08916,1.85256), CC(-.0126777,.0505892)]


v1 = vertex(bp,[x0])
vertices = [v1]

for i in 1:3
    rand_u = [CC(cis(rand()*2*pi)) for _ in 1:length(pars)]
    push!(vertices, vertex(rand_u))
end

compiled_homotopy = compile_edge_homotopy(F, vars, pars)

# ------------------------------------------------------------------------------
# 3. Solve Monodromy (Tracking)
# ------------------------------------------------------------------------------
edges = solve_monodromy(compiled_homotopy, vertices; max_roots=8, show_progress=true)

# ------------------------------------------------------------------------------
# 4. GAP Group Construction
# ------------------------------------------------------------------------------
G = build_gap_group(8, edges) # Find a group of size 8 from edge correspondences

if G !== nothing
    println("Structure Description:")
    println(GAP.Globals.StructureDescription(G)) 
    println("Galois Width:")
    gw = galois_width(G)
    println(gw) # 3
end
