

using Nemo
using AbstractAlgebra
using Symbolics
using GAP
using CertifiedMonodromyComputation

@variables x t
@variables c0,c1,c2,c3,c4,c5,c6,c7,c8

PREC_BITS = 256 
CC = AcbField(PREC_BITS) # Complex Field (acb)
RR = ArbField(PREC_BITS)    # Real Field (arb)
CCi = CC

F_exprs = [
    c8*x^8+c7*x^7+c6*x^6+c5*x^5+c4*x^4+c3*x^3+c2*x^2+c1*x+c0
]
bp = [CCi(1/2), CCi(9/7), CCi(9/7), CCi(45/7), CCi(45/56), CCi(-2/3), CCi(1/3), CCi(3/4), CCi(3/7)]
x0 = [CCi(.0720838,-.481488)]
v1 = vertex(bp,[x0])
n_nodes = 6;


# 나머지 Vertex 생성 (랜덤 파라미터)
# 총 3개의 지점으로 루프를 구성하여 Monodromy 수행
vertices = [v1]

# Reproducibility를 위해 시드 고정 (선택사항)
using Random
Random.seed!() 
x_vars = [x]
p_vars = [c0,c1,c2,c3,c4,c5,c6,c7,c8]

for i in 1:4
    rand_u = [CC(cis(rand()*2*pi)) for _ in 1:9]
    push!(vertices, vertex(rand_u))
end

compiled_homotopy = compile_edge_homotopy(F_exprs, x_vars, p_vars, t; homogeneous=false)

edges = solve_monodromy(compiled_homotopy, vertices; max_roots=8)

G = build_gap_group(8, edges) # Find a group of size 4 from edge correspondences

if G !== nothing
    println("Structure Description:")
    println(GAP.Globals.StructureDescription(G)) 
    println("Galois Width:")
    gw = galois_width(G)
    println(gw) #
end