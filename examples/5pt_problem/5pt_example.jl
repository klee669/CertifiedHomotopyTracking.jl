#Pkg.instantiate()

using CertifiedHomotopyTracking

println("=== Running 5 Point Example ===")

# ------------------------------------------------------------------------------
# 1. Setup Polynomial Ring & System
# ------------------------------------------------------------------------------
# @setupfield to define variables and fields
@monodromy_setup begin
    vars = (X,Y,Z)
    params = (x11,x21,x31,x41,x51,x12,x22,x32,x42,x52,y11,y21,y31,y41,y51,y12,y22,y32,y42,y52)
end
CCi = _CCi

include("eqs.jl")
include("pts.jl")

# ------------------------------------------------------------------------------
# 2. Setup Vertices & Initial Points
# ------------------------------------------------------------------------------
# define base_point 
bp = [CCi(.831187,-.555993), CCi(-.077487,.996993), CCi(-.031994,-.999488), CCi(.019587,.999808), CCi(.949892,-.312579), CCi(-.881637,.471928), CCi(.203765,.97902), CCi(-.819274,.573402), CCi(-.973777,-.227505), CCi(-.972427,-.233208), CCi(.192661,-.981265), CCi(-.856378,-.516349), CCi(.296842,-.954927), CCi(.525122,-.851027), CCi(-.999761,-.021852), CCi(-.679705,-.733486), CCi(-.583596,.812044), CCi(.875249,.483673), CCi(.43059,-.902548), CCi(.438152,.898901)]

x0 = p_list[1] 

# initial vertex
v1 = vertex(bp, [x0])

# construct parameter points (total 4 vertices, 20 points)
vs = parameter_points(v1, 20, 4) 


# ------------------------------------------------------------------------------
# 3. Solve Monodromy (Tracking)
# ------------------------------------------------------------------------------
println("Starting Monodromy Tracking...")

edges = solve_monodromy(F, vs; max_roots=20)

# ------------------------------------------------------------------------------
# 4. GAP Group Construction
# ------------------------------------------------------------------------------
println("Building GAP Group...")

G = build_gap_group(20, edges)

# ------------------------------------------------------------------------------
# 5. Analysis (GAP)
# ------------------------------------------------------------------------------

@gap("Read(\"~/Documents/GitHub/certified_monodromy_comp/examples/5pt_problem/5pt_perm_list_20.txt\");")
@gap("G;")

if G !== nothing
    println("Structure Description:")
    println(GAP.Globals.StructureDescription(G)) # (C2 x C2 x C2 x C2 x C2 x C2 x C2 x C2 x C2) : S10
    println("Galois Width:")
    gw = galois_width(G)
    println(gw) # 10
end