

using CertifiedHomotopyTracking
@variables x y p q
PREC_BITS = 256
CC = AcbField(PREC_BITS)

F_exprs = [p*x^2 + 3*y - 4, y^2 + q]
x_vars = [x, y]
p_vars = [p, q]

bp = [CC(1), CC(-1)] 
x0 = [CC(1) , CC(1)] 


compiled_homotopy = compile_edge_homotopy(F_exprs, x_vars, p_vars; homogeneous=false)

# graph construction
edge_list = [
    (1, 2), (2, 3), (3, 1), 
    (3, 4), (4, 5), (5, 3), 
    (5, 6), (6, 1)          
]

v1 = vertex(bp, [x0])
vertices = [v1]
for i in 1:5 
    rand_u = [CC(cis(rand()*2*pi)) for _ in 1:2]
    push!(vertices, vertex(rand_u))
end
# graph making helper function
custom_edges = build_edges(vertices, edge_list)

# solve monodromy
edges = solve_monodromy(compiled_homotopy, vertices, custom_edges; max_roots=4)

# GAP analysis
G = build_gap_group(4, edges) 
if G !== nothing
    println("Structure Description:")
    println(GAP.Globals.StructureDescription(G)) 
    println("Galois Width:")
    gw = galois_width(G)
    println(gw) 
end