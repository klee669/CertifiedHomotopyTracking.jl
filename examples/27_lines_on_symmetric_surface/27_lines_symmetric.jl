
using CertifiedHomotopyTracking

println("=== Running 27 Lines Example (Manual Loops) ===")

# ------------------------------------------------------------------------------
# 1. Setup System
# ------------------------------------------------------------------------------
@variables a₁ a₂ b₁ b₂ c₁ c₂ d₁ d₂
@variables a2100 a1110
PREC_BITS = 256 
CC = AcbField(PREC_BITS) # Complex Field (acb)


include("27_lines_sol_list.txt") 

# System Definition
f1= a₁*d₁^2*a2100 + a₁^2*d₁*a2100 + b₁*a₁^2*a2100 + b₁*c₁^2*a2100 + b₁*d₁^2*a2100 + b₁^2*a₁*a2100 + b₁^2*c₁*a2100 + b₁^2*d₁*a2100 + c₁*a₁^2*a2100 + c₁*d₁^2*a2100 + c₁^2*a₁*a2100 + c₁^2*d₁*a2100 + b₁*a₁*d₁*a1110 + b₁*c₁*a₁*a1110 + b₁*c₁*d₁*a1110 + c₁*a₁*d₁*a1110 + 1.0*a₁^3 + 1.0*b₁^3 + 1.0*c₁^3 + 1.0*d₁^3
f2=3.0*a₁^2*a₂ + 3.0*b₂*b₁^2 + 3.0*c₁^2*c₂ + 3.0*d₂*d₁^2 + a₁^2*d₂*a2100 + a₂*d₁^2*a2100 + b₁^2*a₂*a2100 + b₁^2*c₂*a2100 + b₁^2*d₂*a2100 + b₂*a₁^2*a2100 + b₂*c₁^2*a2100 + b₂*d₁^2*a2100 + c₁^2*a₂*a2100 + c₁^2*d₂*a2100 + c₂*a₁^2*a2100 + c₂*d₁^2*a2100 + 2*a₁*a₂*d₁*a2100 + 2*a₁*d₂*d₁*a2100 + 2*b₁*a₁*a₂*a2100 + b₁*a₁*d₂*a1110 + b₁*a₂*d₁*a1110 + b₁*c₁*a₂*a1110 + 2*b₁*c₁*c₂*a2100 + b₁*c₁*d₂*a1110 + b₁*c₂*a₁*a1110 + b₁*c₂*d₁*a1110 + 2*b₁*d₂*d₁*a2100 + b₂*a₁*d₁*a1110 + 2*b₂*b₁*a₁*a2100 + 2*b₂*b₁*c₁*a2100 + 2*b₂*b₁*d₁*a2100 + b₂*c₁*a₁*a1110 + b₂*c₁*d₁*a1110 + 2*c₁*a₁*a₂*a2100 + c₁*a₁*d₂*a1110 + c₁*a₂*d₁*a1110 + 2*c₁*c₂*a₁*a2100 + 2*c₁*c₂*d₁*a2100 + 2*c₁*d₂*d₁*a2100 + c₂*a₁*d₁*a1110
f3= 3.0*a₁*a₂^2 + 3.0*b₂^2*b₁ + 3.0*c₁*c₂^2 + 3.0*d₂^2*d₁ + a₁*d₂^2*a2100 + a₂^2*d₁*a2100 + b₁*a₂^2*a2100 + b₁*c₂^2*a2100 + b₁*d₂^2*a2100 + b₂^2*a₁*a2100 + b₂^2*c₁*a2100 + b₂^2*d₁*a2100 + c₁*a₂^2*a2100 + c₁*d₂^2*a2100 + c₂^2*a₁*a2100 + c₂^2*d₁*a2100 + 2*a₁*a₂*d₂*a2100 + 2*a₂*d₂*d₁*a2100 + b₁*a₂*d₂*a1110 + b₁*c₂*a₂*a1110 + b₁*c₂*d₂*a1110 + 2*b₂*a₁*a₂*a2100 + b₂*a₁*d₂*a1110 + b₂*a₂*d₁*a1110 + 2*b₂*b₁*a₂*a2100 + 2*b₂*b₁*c₂*a2100 + 2*b₂*b₁*d₂*a2100 + b₂*c₁*a₂*a1110 + 2*b₂*c₁*c₂*a2100 + b₂*c₁*d₂*a1110 + b₂*c₂*a₁*a1110 + b₂*c₂*d₁*a1110 + 2*b₂*d₂*d₁*a2100 + c₁*a₂*d₂*a1110 + 2*c₁*c₂*a₂*a2100 + 2*c₁*c₂*d₂*a2100 + 2*c₂*a₁*a₂*a2100 + c₂*a₁*d₂*a1110 + c₂*a₂*d₁*a1110 + 2*c₂*d₂*d₁*a2100
f4=a₂*d₂^2*a2100 + a₂^2*d₂*a2100 + b₂*a₂^2*a2100 + b₂*c₂^2*a2100 + b₂*d₂^2*a2100 + b₂^2*a₂*a2100 + b₂^2*c₂*a2100 + b₂^2*d₂*a2100 + c₂*a₂^2*a2100 + c₂*d₂^2*a2100 + c₂^2*a₂*a2100 + c₂^2*d₂*a2100 + b₂*a₂*d₂*a1110 + b₂*c₂*a₂*a1110 + b₂*c₂*d₂*a1110 + c₂*a₂*d₂*a1110 + 1.0*a₂^3 + 1.0*b₂^3 + 1.0*c₂^3 + 1.0*d₂^3
f5= -1.0 - 0.506678639976439*a₁ - 1.88147657425506*a₂ + 2.98140090140213*b₁ - 0.0662841633453234*b₂ + 0.64589883691229*c₁ + 1.20555402987922*c₂ - 0.967371765008337*d₁ + 0.59080240629979*d₂
f6=-1.0 - 0.614671350751522*a₁ - 1.13214950346366*a₂ - 1.45701065536224*b₁ - 0.927555894895465*b₂ + 0.521704930167087*c₁ - 0.126225364275008*c₂ + 0.861780657701936*d₁ + 1.10700700597821*d₂
f7=-1.0 + 0.129606622277643*a₁ + 0.876163537518904*a₂ - 0.190549680783866*b₁ + 1.62684419040138*b₂ + 0.284493812409805*c₁ + 2.12502842074383*c₂ - 1.54300612629157*d₁ + 0.827955755122909*d₂
f8=-1.0 + 1.30252394968029*a₁ - 0.283220423701744*a₂ - 2.22448020204455*b₁ + 2.07503960766641*b₂ - 1.00077453580414*c₁ - 0.0136695870565606*c₂ - 1.33893635772738*d₁ - 1.1143731550125*d₂
f = [f1, f2, f3, f4, f5, f6, f7, f8];
x_vars = [a₁, a₂, b₁, b₂, c₁, c₂, d₁, d₂]
p_vars = [a2100, a1110]
# ------------------------------------------------------------------------------
# 2. Local Helper Functions
# ------------------------------------------------------------------------------
function track_loop(bp, a, b, x0, p_list, i, F)
    println("Root Number $i: Tracking the first edge")
    F1 = make_edge_system(F, bp, a)
    res_x1 = track_path(F1, x0; t_end=1.0, h_init=0.1, show_progress=true)
    x1 = certified_region(res_x1)
    if !succeeded(res_x1) return nothing, nothing end 

    println("Root Number $i: Tracking the second edge")
    F2 = make_edge_system(F, a, b)
    res_x2 = track_path(F2, x1; t_end=1.0, h_init=0.1, show_progress=true)
    x2 = certified_region(res_x2)
    if !succeeded(res_x2) return nothing, nothing end 

    println("Root Number $i: Tracking the third edge")
    F3 = make_edge_system(F, b, bp)
    res_x3 = track_path(F3, x2; t_end=1.0, h_init=0.1, show_progress=true)
    x3 = certified_region(res_x3)
    if !succeeded(res_x3) return nothing, nothing end 

    ind = search_point_certified(F3, x3, p_list)
    println("Result: Mapped to $ind")
    return x3, ind
end

function generate_perm(F, bp, a, b, p_list)
    n = length(p_list)
    perm = []
    res_list = []
    for i = 1:n
        res, ind = track_loop(bp, a, b, p_list[i], p_list, i, F)
        
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
# red loop
red1 = [CC(-500),CC(0,0)]
red2 = [CC(-500),CC(3200,-500)]
red3 = [CC(-500),CC(3200,500)]
v1 = vertex(red1, [p_list[1]])
v2 = vertex(red2)
v3 = vertex(red3)

vertices = [v1, v2, v3]

compiled_homotopy = compile_edge_homotopy(f, x_vars, p_vars)


println("\n[Loop 1] Red Loop")
p1 = generate_perm(compiled_homotopy, red1, red2, red3, p_list) #[8, 2, 3, 4, 24, 6, 7, 1, 27, 23, 11, 16, 13, 14, 25, 12, 17, 18, 19, 20, 21, 22, 10, 5, 15, 26, 9]


# green loop
green1 = [CC(-500),CC(0,0)]
green2 = [CC(-500),CC(-800,-500)]
green3 = [CC(-500),CC(-800,500)]
println("\n[Loop 2] Green Loop")
p2 = generate_perm(compiled_homotopy, green1, green2, green3, p_list) #[25, 4, 17, 2, 27, 11, 26, 15, 24, 12, 6, 10, 14, 13, 8, 23, 3, 18, 19, 22, 21, 20, 16, 9, 1, 7, 5]



# ------------------------------------------------------------------------------
# 4. GAP Analysis
# ------------------------------------------------------------------------------
println("\n=== GAP Analysis ===")

p1_gap = GAP.Globals.PermList(GAP.Obj(p1[2]))
p2_gap = GAP.Globals.PermList(GAP.Obj(p2[2]))

prod = p1_gap * p2_gap
println("p1 * p2 calculated in GAP.")

G = GAP.Globals.Group(p1_gap, p2_gap)
println("Group G defined.")

println("Structure Description:")
println(GAP.Globals.StructureDescription(G)) # C2 x C2 (= K4)

println("Galois Width:")
gw = galois_width(G) # 2
