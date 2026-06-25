using CertifiedHomotopyTracking
const CHT = CertifiedHomotopyTracking

HomotopyContinuation.@var a_1, a_2, b_1, b_2, c_1, c_2, d_1, d_2

a=0

f = HomotopyContinuation.System([
    a_1^3+b_1^3+c_1^3+d_1^3,
    3*a_1^2*a_2+3*b_1^2*b_2+3*c_1^2*c_2+3*d_1^2*d_2,
    3*a_1*a_2^2+3*b_1*b_2^2+3*c_1*c_2^2+3*d_1*d_2^2,
    a_2^3+b_2^3+c_2^3+d_2^3,
    .111589*a_1+1.48565*a_2+2.523357*b_1-1.646118*b_2+.30993*c_1+.36067*c_2+.0452546*d_1+.429238*d_2+.0146597,
    .896976*a_1+.138818*a_2+.707287*b_1-.0657292*b_2+.26806*c_1+.866177*c_2+.18616*d_1+.896809*d_2+.417192,
    .473248*a_1+.643963*a_2+.12494*b_1+.520899*b_2+.576117*c_1+.486441*c_2+.384317*d_1+.552473*d_2+.18095,
    .779979*a_1+.443068*a_2+.901413*b_1+2.349995*b_2+.283387*c_1+.646477*c_2+.549256*d_1+.135106*d_2+.569817
])

sols = HomotopyContinuation.solve(f)
nonsingular_sols = HomotopyContinuation.solutions(sols, only_nonsingular=true)

# ------------------------------------------------------------------------------
# 1. Setup System
# ------------------------------------------------------------------------------
@variables a_1 a_2 b_1 b_2 c_1 c_2 d_1 d_2
@variables alpha
const PREC_BITS = 256 
const CC = AcbField(PREC_BITS) # Complex Field (acb)

# convert solutions to CC format
p_list = map(i-> map(j -> CC(j), i), nonsingular_sols)

v1 = vertex([CC(0)], [p_list[1]])



# System Definition
f = [a_1^3+b_1^3+c_1^3+6*alpha*b_1*c_1*d_1+d_1^3,
    3*a_1^2*a_2+3*b_1^2*b_2+3*c_1^2*c_2+6*alpha*b_2*c_1*d_1+6*alpha*b_1*c_2*d_1+6*alpha*b_1*c_1*d_2+3*d_1^2*d_2,
     3*a_1*a_2^2+3*b_1*b_2^2+3*c_1*c_2^2+6*alpha*b_2*c_2*d_1+6*alpha*b_2*c_1*d_2+6*alpha*b_1*c_2*d_2+3*d_1*d_2^2,
     a_2^3+b_2^3+c_2^3+6*alpha*b_2*c_2*d_2+d_2^3,
    #.911589*a_1-.48565*a_2+2.523357*b_1-.646118*b_2+.30993*c_1+.36067*c_2+.0452546*d_1+.429238*d_2+.0146597,
    .111589*a_1+1.48565*a_2+2.523357*b_1-1.646118*b_2+.30993*c_1+.36067*c_2+.0452546*d_1+.429238*d_2+.0146597,
    .896976*a_1+.138818*a_2+.707287*b_1-.0657292*b_2+.26806*c_1+.866177*c_2+.18616*d_1+.896809*d_2+.417192,
    .473248*a_1+.643963*a_2+.12494*b_1+.520899*b_2+.576117*c_1+.486441*c_2+.384317*d_1+.552473*d_2+.18095,
    .779979*a_1+.443068*a_2+.901413*b_1+2.349995*b_2+.283387*c_1+.646477*c_2+.549256*d_1+.135106*d_2+.569817
    ];
x_vars = [a_1, a_2, b_1, b_2, c_1, c_2, d_1, d_2]
p_vars = [alpha]

# ------------------------------------------------------------------------------
# 2. Local Helper Functions
# ------------------------------------------------------------------------------
function _track_lines27_edge(F, p_start, p_target, x_start; posteriori = false, posteriori_options = (;), show_progress = true)
    sys = make_edge_system(F, p_start, p_target)
    if posteriori
        cert = certify_posteriori(
            sys,
            x_start;
            show_progress = show_progress,
            posteriori_options...,
        )
        return cert.success ? certified_region(cert) : nothing, cert
    end

    result = track_path(sys, x_start; t_end = 1.0, h_init = 0.1, show_progress = show_progress)
    return success(result) ? certified_region(result) : nothing, result
end

function track_loop(bp, a, b, x0, p_list, i, F; posteriori = false, posteriori_options = (;), show_progress = true)
    println("Root Number $i: Tracking the first edge")
    x1, res_x1 = _track_lines27_edge(
        F,
        bp,
        a,
        x0;
        posteriori = posteriori,
        posteriori_options = posteriori_options,
        show_progress = show_progress,
    )
    x1 === nothing && return nothing, nothing

    println("Root Number $i: Tracking the second edge")
    x2, res_x2 = _track_lines27_edge(
        F,
        a,
        b,
        x1;
        posteriori = posteriori,
        posteriori_options = posteriori_options,
        show_progress = show_progress,
    )
    x2 === nothing && return nothing, nothing

    println("Root Number $i: Tracking the third edge")
    x3, res_x3 = _track_lines27_edge(
        F,
        b,
        bp,
        x2;
        posteriori = posteriori,
        posteriori_options = posteriori_options,
        show_progress = show_progress,
    )
    x3 === nothing && return nothing, nothing

    F3 = make_edge_system(F, b, bp)
    ind = CHT.search_point_certified(F3, x3, p_list)
    println("Result: Mapped to $ind")
    return x3, ind
end

function generate_perm(F, bp, a, b, p_list; posteriori = false, posteriori_options = (;), show_progress = true)
    n = length(p_list)
    perm = []
    res_list = []
    for i = 1:n
        res, ind = track_loop(
            bp,
            a,
            b,
            p_list[i],
            p_list,
            i,
            F;
            posteriori = posteriori,
            posteriori_options = posteriori_options,
            show_progress = show_progress,
        )
        
        if res === nothing 
            println("Stopped by user. Returning partial permutation.")
            break 
        end
        
        push!(perm, ind)
        push!(res_list, res)
    end
    return res_list, perm
end

# ------------------------------------------------------------------------------
# 3. Manual Loop Definitions
# ------------------------------------------------------------------------------


starting_point = [CC(0)]
leftloop1 = [CC(-2,5)]
leftloop2 = [CC(-2,-5)]

USE_POSTERIORI = true
POSTERIORI_OPTIONS = (;
    max_depth = 12,
)

p1 = generate_perm(
    compiled_homotopy,
    starting_point,
    leftloop1,
    leftloop2,
    p_list;
    posteriori = USE_POSTERIORI,
    posteriori_options = POSTERIORI_OPTIONS,
);


urloop1 = [CC(-2,5)]
urloop2 = [CC(10,0)]


p2 = generate_perm(
    compiled_homotopy,
    starting_point,
    urloop1,
    urloop2,
    p_list;
    posteriori = USE_POSTERIORI,
    posteriori_options = POSTERIORI_OPTIONS,
); 



drloop1 = [CC(-2,-5)]
drloop2 = [CC(10,0)]


p3 = generate_perm(
    compiled_homotopy,
    starting_point,
    drloop1,
    drloop2,
    p_list;
    posteriori = USE_POSTERIORI,
    posteriori_options = POSTERIORI_OPTIONS,
);


# ------------------------------------------------------------------------------
# 4. GAP Analysis
# ------------------------------------------------------------------------------
println("\n=== GAP Analysis ===")

p1_gap = GAP.Globals.PermList(GAP.Obj(p1[2]))
p2_gap = GAP.Globals.PermList(GAP.Obj(p2[2]))
p3_gap = GAP.Globals.PermList(GAP.Obj(p3[2]))

G = GAP.Globals.Group(p1_gap, p2_gap,p3_gap)
println("Group G defined.")
println("Structure Description:")
println(GAP.Globals.StructureDescription(G)) # S3
println("Galois Width:")
gw = galois_width(G) # 3
