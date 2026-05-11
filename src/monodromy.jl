# ------------------------------------------------------------------------------
# High-Level Interface
# ------------------------------------------------------------------------------
# Exporting Vertex/Edge allows users to inspect results and build custom graphs
export solve_monodromy, build_gap_group, vertex, edge, parameter_points, galois_width,
    Edge, Vertex, search_point, search_point_certified, same_root_krawczyk, build_edges

# ------------------------------------------------------------------------------
# Data Structures
# ------------------------------------------------------------------------------
const PointType = Vector{AcbFieldElem} 

mutable struct Vertex
    base_point::PointType
    sols::Vector{PointType}
    Edges::Vector{Any} # Type Any to avoid circular dependency
end

mutable struct Edge
    node1::Vertex
    node2::Vertex
    correspondence12::Vector{Tuple{Int, Int}}
    correspondence21::Vector{Tuple{Int, Int}}
end

# Constructors
vertex(p::Vector) = Vertex(p, PointType[], [])
vertex(p::Vector, x::Vector) = Vertex(p, x, [])
edge(p::Vertex, q::Vertex) = Edge(p, q, Tuple{Int,Int}[], Tuple{Int,Int}[])

Base.show(io::IO, v::Vertex) = print(io, "Vertex($(length(v.sols)) solutions)")
Base.show(io::IO, e::Edge) = print(io, "Edge($(length(e.correspondence12)) correspondences)")

function build_edges(vertices::Vector{Vertex}, edge_pairs::Vector{Tuple{Int, Int}})
    custom_edges = Edge[]
    for (i, j) in edge_pairs
        e = edge(vertices[i], vertices[j])
        push!(custom_edges, e)
        push!(vertices[i].Edges, e)
        push!(vertices[j].Edges, e)
    end
    return custom_edges
end


function make_edge_system(
    compiled_sys::CompiledHomotopy,
    p_start::Vector{AcbFieldElem},
    p_end::Vector{AcbFieldElem},
    p_const::Vector{AcbFieldElem}=AcbFieldElem[];
    patch_vector::Vector{AcbFieldElem}=AcbFieldElem[],
)
    return HCSystem(compiled_sys, p_start, p_end, p_const; patch_vector=patch_vector)
end


# ------------------------------------------------------------------------------
# Core Logic: Monodromy Solving
# ------------------------------------------------------------------------------

"""
    solve_monodromy(H, vertices; radius=0.1, max_roots=20, predictor=true, show_progress=true)

Tracks the complete graph to find the monodromy group.
Handles InterruptException (Ctrl+C) gracefully to return partial results.
Updates the `vertices` and edges in-place.
"""
function solve_monodromy(
    H::Union{Matrix,Vector},
    vertices::Vector{Vertex};
    radius::Number = 0.1,
    max_roots::Int = 20,
    predictor::Bool = true,
    show_progress::Bool = true,
)
    # 1. Initialize Edges (Connect all vertices)
    edges = Edge[]
    for i in 1:length(vertices)-1
        for j in i+1:length(vertices)
            e = edge(vertices[i], vertices[j])
            push!(edges, e)
            push!(vertices[i].Edges, e)
            push!(vertices[j].Edges, e)
        end
    end
    
    base_points = map(v -> v.base_point, vertices)
    
    # 2. Main Loop with Interrupt Handling
    try
        iter_stagnant = 0
        total_correspondences = 0
        
        while true
            current_correspondences = sum(map(e -> length(e.correspondence12), edges))
            
            # Check progress
            if current_correspondences == total_correspondences
                iter_stagnant += 1
            else
                iter_stagnant = 0
                total_correspondences = current_correspondences
            end

            # Stopping conditions
            if all(e -> length(e.correspondence12) == max_roots, edges)
                @info "Success: All edges reached $max_roots correspondences."
                break
            end

            if iter_stagnant > 10
                @warn "Stopping: No new solutions found after 10 iterations."
                break
            end

            # Track each edge
            for (idx, e) in enumerate(edges)
                n1_idx = search_point(e.node1.base_point, base_points)
                n2_idx = search_point(e.node2.base_point, base_points)
                
                # Track Forward & Backward
                track_edge!(H, e, true, radius; predictor=predictor, show_progress=show_progress, id1=n1_idx, id2=n2_idx)
                track_edge!(H, e, false, radius; predictor=predictor, show_progress=show_progress, id1=n2_idx, id2=n1_idx)
                
                # Log status in real-time
                counts = map(x -> length(x.correspondence12), edges)
                @info "Status (Edge $idx done): Correspondences => $counts"
            end
        end

    catch e
        # [Core Feature] Handle Ctrl+C to preserve data
        if isa(e, InterruptException)
            @warn "⚠️ Computation Interrupted by User!"
            @warn "Returning vertices and edges computed SO FAR."
            @warn "You can resume or analyze the partial data."
            return edges
        else
            # Re-throw other unexpected errors
            rethrow(e)
        end
    end

    return edges
end

"""
    track_edge!(H, edge, direction, radius; predictor, show_progress, id1, id2)

Tracks paths along an edge. Updates the edge and vertices IN-PLACE.
"""
function track_edge!(
    H::Union{Matrix, Vector}, 
    e::Edge, 
    from1to2::Bool, 
    r::Number;
    predictor = true,
    show_progress = true,
    id1 = "?", # For logging only
    id2 = "?"
)
    # Determine direction
    if from1to2
        source_v, target_v = e.node1, e.node2
        c_forward, c_backward = e.correspondence12, e.correspondence21
    else
        source_v, target_v = e.node2, e.node1
        c_forward, c_backward = e.correspondence21, e.correspondence12
    end

    source_sols = source_v.sols
    target_sols = target_v.sols
    
    # Identify which paths haven't been tracked yet
    tracked_indices = Set(x[1] for x in c_forward)
    untracked_indices = [i for i in 1:length(source_sols) if !(i in tracked_indices)]

    if isempty(untracked_indices)
        return
    end

    @info "Tracking Edge ($id1 -> $id2): $(length(untracked_indices)) new paths to track."

    # Create the homotopy system for this specific edge
    # Note: specified_system must be exported from Homotopy module
    Fab = specified_system(source_v.base_point, target_v.base_point, H)

    for src_idx in untracked_indices
        start_point = source_sols[src_idx]
        
        # Perform Tracking
        if predictor
            y = track(Fab, start_point; show_display = show_progress, refinement_threshold = 1/8)
        else
            y = tracking_without_predictor(Fab, start_point)
        end

        # Check if the result 'y' is a known point in target_v
        dest_idx = search_point(y, target_sols)
        
        if dest_idx === nothing
            # New solution found! Add it to the target vertex
            push!(target_sols, y)
            dest_idx = length(target_sols)
        end

        # Record the correspondence
        push!(c_forward, (src_idx, dest_idx))
        push!(c_backward, (dest_idx, src_idx))
    end
    
    # Sort correspondences for consistency
    sort!(c_forward)
    sort!(c_backward)
end


function track_edge!(
    compiled_sys::CompiledHomotopy,
    e::Edge,
    from1to2::Bool,
    edge_idx::Int,
    id_src::Int,
    id_dest::Int;
    show_progress::Bool=false,
    root_match::Symbol=:heuristic,
)
    if from1to2
        source_v, target_v = e.node1, e.node2
        c_forward, c_backward = e.correspondence12, e.correspondence21
    else
        source_v, target_v = e.node2, e.node1
        c_forward, c_backward = e.correspondence21, e.correspondence12
    end

    source_sols = source_v.sols
    target_sols = target_v.sols
    
    tracked_indices = Set(x[1] for x in c_forward)
    untracked_indices = [i for i in 1:length(source_sols) if !(i in tracked_indices)]

    if isempty(untracked_indices)
        return
    end

    @info "Tracking Edge ($id_src -> $id_dest): $(length(untracked_indices)) new paths to track."

    sys_edge = make_edge_system(compiled_sys, source_v.base_point, target_v.base_point)
    
    count_success = 0; count_new = 0; count_collision = 0
    
    for src_idx in untracked_indices
        start_point = source_sols[src_idx]
        
        y_end, success = track_path(sys_edge, start_point; t_end=1.0, h_init=0.1, show_progress=show_progress)
        
        if success
            dest_idx = _search_solution(sys_edge, y_end, target_sols; root_match=root_match)
            
            if dest_idx === nothing
                push!(target_sols, y_end)
                dest_idx = length(target_sols)
                
                push!(c_forward, (src_idx, dest_idx))
                push!(c_backward, (dest_idx, src_idx))
                
                printstyled("+", color=:green, bold=true)
                count_new += 1
                count_success += 1
            else
                already_mapped = any(p -> p[2] == dest_idx, c_forward)
                if !already_mapped
                    push!(c_forward, (src_idx, dest_idx))
                    push!(c_backward, (dest_idx, src_idx))
                    printstyled(".", color=:cyan)
                    count_success += 1
                else
                    printstyled("c", color=:yellow, bold=true)
                    count_collision += 1
                end
            end
        else
            printstyled("x", color=:red, bold=true)
        end
    end
    
    sort!(c_forward)
    sort!(c_backward)
    println(" Done. (Ok: $count_success, New: $count_new, Collision: $count_collision)")
end



function solve_monodromy(
    compiled_sys::CompiledHomotopy,
    vertices::Vector{Vertex},
    edges::Vector{Edge};
    max_roots=20,
    show_progress::Bool=false,
    root_match::Symbol=:heuristic,
)
    iter = 0
    iter_stagnant = 0
    total_correspondences = 0
    
    try
        while true
            iter += 1
            current_correspondences = sum(length(e.correspondence12) for e in edges)
            
            if current_correspondences == total_correspondences
                iter_stagnant += 1
            else
                iter_stagnant = 0
                total_correspondences = current_correspondences
            end
            
            max_sols_found = maximum(length(v.sols) for v in vertices)
            total_sols_stored = sum(length(v.sols) for v in vertices)
            
            println("\n" * "="^70)
            println(" Iteration $iter | Target: $max_roots | Max Found: $max_sols_found | Stored: $total_sols_stored | Stagnant: $iter_stagnant")
            println("-"^70)
            
            if all(e -> length(e.correspondence12) == max_roots, edges)
                println("🎉 Success! All edges reached $max_roots correspondences.")
                break
            end

            if iter_stagnant > 10
                println("⚠️ Stopping: No new solutions found after 10 iterations.")
                break
            end
            
            for (idx, e) in enumerate(edges)
                id1 = findfirst(==(e.node1), vertices)
                id2 = findfirst(==(e.node2), vertices)
                
                track_edge!(compiled_sys, e, true, idx, id1, id2; show_progress=show_progress, root_match=root_match)
                track_edge!(compiled_sys, e, false, idx, id2, id1; show_progress=show_progress, root_match=root_match)

                counts = map(x -> length(x.correspondence12), edges)
                @info "Status (Edge $idx done): Correspondences => $counts"
            end
        end
        
    catch e
        if isa(e, InterruptException)
            println("\n\n⚠️ Computation Interrupted by User!")
            println("Returning vertices and edges computed SO FAR.")
            return edges
        else
            rethrow(e)
        end
    end
    
    return edges
end
function solve_monodromy(
    compiled_sys::CompiledHomotopy,
    vertices::Vector{Vertex};
    max_roots=20,
    show_progress::Bool=false,
    root_match::Symbol=:heuristic,
)
    println("Building a complete graph for the given vertices...")
    edges = Edge[]
    for i in 1:length(vertices)-1
        for j in i+1:length(vertices)
            e = edge(vertices[i], vertices[j])
            push!(edges, e)
            push!(vertices[i].Edges, e)
            push!(vertices[j].Edges, e)
        end
    end
    
    return solve_monodromy(compiled_sys, vertices, edges; max_roots=max_roots, show_progress=show_progress, root_match=root_match)
end

# ------------------------------------------------------------------------------
# GAP Integration
# ------------------------------------------------------------------------------

"""
    build_gap_group(max_roots, edges)

Constructs a GAP permutation group directly from the tracking results.
Returns a GAP Group object.
"""
function build_gap_group(max_roots::Int, edges::Vector{Edge})
    perms_vec = get_permutations(max_roots, edges)
    
    if isempty(perms_vec)
        @warn "No valid permutations found to build a group."
        return nothing
    end

    # Use GAP.Obj to convert Julia Vector to GAP List
    gap_perms = [GAP.Globals.PermList(GAP.Obj(p)) for p in perms_vec]
    
    # Create the group in GAP
    G = GAP.Globals.Group(gap_perms...)
    
    return G
end

"""
    galois_width(G)

Calculates the Galois Width of the GAP group G using GAP.
Defines the function in GAP unconditionally to avoid scope issues.
"""
function galois_width(G::GAP.GapObj)
    # Define the function in GAP (unconditional definition)
    @gap("""
    GaloisWidth := function(G)
      local X, M, C, phi;
      if IsTrivial(G) then return 1;
      elif IsNaturalSymmetricGroup(G) or IsNaturalAlternatingGroup(G) then
        X := OrbitsDomain(G)[1];
        if Length(X) = 4 then return 3;
        else return Length(X);
        fi;
      elif IsCyclic(G) then return Maximum(Factors(Order(G)));
      elif not IsTransitive(G) then 
        return Maximum(List(Orbits(G), 
          O -> GaloisWidth(Image(ActionHomomorphism(G,O)))
        ));
      else
        X := OrbitsDomain(G)[1];
        if not IsPrimitive(G) then
          phi := ActionHomomorphism(G,Blocks(G, X),OnSets);
          return Maximum(GaloisWidth(Kernel(phi)), GaloisWidth(Image(phi)));
        elif IsSimple(G) then
          M := List(ConjugacyClassesMaximalSubgroups(G), H -> Representative(H));
          return Minimum(List(M, H -> Order(G)/Order(H)));
        else
          C := CompositionSeries(G);
          return Maximum(List([1..Length(C)-1], 
            i -> GaloisWidth(C[i]/C[i+1])
          ));
        fi;
      fi;
    end;
    """)

    return GAP.Globals.GaloisWidth(G)
end


# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------

function search_point(res::Vector{AcbFieldElem}, pool::Vector{Vector{AcbFieldElem}}; tol=1e-4)
    for (idx, sol) in enumerate(pool)
        dist = maximum(mag_complex.(res .- sol))
        if dist < tol
            return idx
        end
    end
    return nothing
end

function _tracking_coordinates(sys::HCSystem, x::Vector{AcbFieldElem})
    if has_projective_patch(sys) && length(x) == length(sys.patch_vector) - 1
        return lift_to_patch(x, collect(sys.patch_vector))
    elseif sys.homogeneous && !has_projective_patch(sys) && sys.compiled.n_vars != 0 && length(x) == sys.compiled.n_vars - 1
        return [sys.CC(1); x]
    end
    return copy(x)
end

function _enclosing_center_radius(x::Vector{AcbFieldElem}, y::Vector{AcbFieldElem}; inflate=1.1, min_radius=1e-12)
    length(x) == length(y) || throw(DimensionMismatch("Cannot compare roots with different dimensions."))
    c = get_mid_vec((get_mid_vec(x) .+ get_mid_vec(y)) ./ 2)
    return c, _enclosing_radius_about(c, x, y; inflate=inflate, min_radius=min_radius)
end

function _enclosing_radius_about(center::Vector{AcbFieldElem}, xs::Vector{AcbFieldElem}...; inflate=1.1, min_radius=1e-12)
    r = min_radius
    for x in xs
        length(center) == length(x) || throw(DimensionMismatch("Cannot compare roots with different dimensions."))
        for i in eachindex(center)
            r = max(r, inflate * max_int_norm(x[i] - center[i]))
        end
    end
    return r
end

function _krawczyk_same_on_hull(sys::HCSystem, x::Vector{AcbFieldElem}, y::Vector{AcbFieldElem}, t; rho, inflate, min_radius)
    center, radius = _enclosing_center_radius(x, y; inflate=inflate, min_radius=min_radius)
    A = compute_preconditioner(sys, center, t)
    passed, k_norm = krawczyk_test(sys, center, t, radius, A; rho=rho)
    return passed, k_norm, radius
end

function _krawczyk_same_at_center(sys::HCSystem, center::Vector{AcbFieldElem}, x::Vector{AcbFieldElem}, y::Vector{AcbFieldElem}, t; rho, inflate, min_radius)
    radius = _enclosing_radius_about(center, x, y; inflate=inflate, min_radius=min_radius)
    A = compute_preconditioner(sys, center, t)
    passed, k_norm = krawczyk_test(sys, center, t, radius, A; rho=rho)
    return passed, k_norm, radius
end

function _polish_for_root_comparison(sys::HCSystem, x::Vector{AcbFieldElem}, t; radius=1e-8)
    center = get_mid_vec(x)
    A = compute_preconditioner(sys, center, t)
    polished, _, _, success = refine_moore_box(sys, center, t, radius, A)
    return success ? get_mid_vec(polished) : center, success
end

function same_root_krawczyk(
    sys::HCSystem,
    x::Vector{AcbFieldElem},
    y::Vector{AcbFieldElem};
    t=1.0,
    rho=0.7,
    inflate=1.1,
    min_radius=1e-12,
)
    x_tracking = _tracking_coordinates(sys, x)
    y_tracking = _tracking_coordinates(sys, y)

    try
        k_norm = Inf
        radius = Inf
        radius_floors = (min_radius, max(min_radius, 1e-10), max(min_radius, 1e-8))
        for radius_floor in radius_floors
            passed, k_norm, radius = _krawczyk_same_on_hull(sys, x_tracking, y_tracking, t; rho=rho, inflate=inflate, min_radius=radius_floor)
            passed && return :same, k_norm, radius
        end

        x_polished, x_success = _polish_for_root_comparison(sys, x_tracking, t)
        y_polished, y_success = _polish_for_root_comparison(sys, y_tracking, t)
        if x_success || y_success
            centers = (
                get_mid_vec((x_polished .+ y_polished) ./ 2),
                x_polished,
                y_polished,
            )
            for radius_floor in radius_floors
                for center in centers
                    passed, k_norm, radius = _krawczyk_same_at_center(sys, center, x_polished, y_polished, t; rho=rho, inflate=inflate, min_radius=radius_floor)
                    passed && return :same, k_norm, radius
                end
            end
        end

        return :unknown, k_norm, radius
    catch
        return :unknown, Inf, Inf
    end
end

function search_point_certified(
    sys::HCSystem,
    res::Vector{AcbFieldElem},
    pool::Vector{Vector{AcbFieldElem}};
    t=1.0,
    rho=0.7,
    inflate=1.1,
    min_radius=1e-12,
)
    for (idx, sol) in enumerate(pool)
        status, _, _ = same_root_krawczyk(sys, res, sol; t=t, rho=rho, inflate=inflate, min_radius=min_radius)
        status == :same && return idx
    end
    return nothing
end

function _search_solution(sys::HCSystem, res::Vector{AcbFieldElem}, pool::Vector{Vector{AcbFieldElem}}; root_match::Symbol=:heuristic)
    if root_match == :heuristic
        return search_point(res, pool)
    elseif root_match == :certified
        return search_point_certified(sys, res, pool)
    elseif root_match == :certified_or_heuristic
        idx = search_point_certified(sys, res, pool)
        return idx === nothing ? search_point(res, pool) : idx
    else
        throw(ArgumentError("Unknown root_match=$root_match. Use :heuristic, :certified, or :certified_or_heuristic."))
    end
end

function parameter_points(v1::Vertex, sz_p::Int, n_vertices::Int)
    # Detect the correct ring from the input vertex
    CC = parent(v1.base_point[1]) 
    
    vertices = Vertex[v1]
    for i in 1:n_vertices-1
        v = AcbFieldElem[]
        for j in 1:sz_p
            # Random point on unit circle
            r_unit_circle = exp(rand(Int8) * im)
            push!(v, CC(real(r_unit_circle), imag(r_unit_circle)))
        end
        push!(vertices, vertex(v))
    end
    vertices
end

# ------------------------------------------------------------------------------
# Permutation Logic
# ------------------------------------------------------------------------------

function complete_correspondences(rc::Int, E::Vector{Edge})
    # Filter edges that are not fully tracked (must match max_roots 'rc')
    filter(e -> 
        length(unique(map(j -> j[2], e.correspondence12))) == rc && 
        length(map(j -> j[2], e.correspondence21)) == rc, E)
end

function neighbor(v::Vertex, e::Edge)
    if v == e.node1
        return e.node2
    elseif v == e.node2
        return e.node1
    else
        error("Edge is not incident at the given vertex.")
    end
end

function membership_test(v::Vertex, e::Edge, v_list::Vector{Vertex}, e_list::Vector{Edge})
    u = neighbor(v, e)
    !(u in v_list) && (e in e_list)
end

function p_compose(
    H1::Union{Vector{Tuple{Int64,Int64}}, Vector{Pair{Int64,Int64}}}, 
    H2::Vector{Pair{Int,Int}}
)
    l = length(H2)
    sorted_H1 = sort(H1) 
    
    map(k -> k => sorted_H1[H2[k][2]][2], 1:l)
end

function get_permutations(rc::Number, E::Vector{Edge})
    id_perm = [i => i for i in 1:rc]
    
    # 1. Get a subgraph of fully tracked edges
    EG = complete_correspondences(rc, E)
    if isempty(EG)
        return Vector{Int64}[]
    end

    VG = unique(vcat(map(i -> i.node1, EG), map(i -> i.node2, EG)))

    # Spanning tree construction logic
    uncovered_v = deleteat!(copy(VG), findall(i -> i == VG[1], VG))
    uncovered_e = copy(EG)

    T = [] # Spanning Tree

    while !isempty(uncovered_v)
        # Find a vertex in uncovered_v connected to the current tree
        v_candidates = filter(v -> any(e -> membership_test(v, e, uncovered_v, uncovered_e), v.Edges), uncovered_v)
        
        if !isempty(v_candidates)
            v = v_candidates[1]
            e_candidates = filter(e -> membership_test(v, e, uncovered_v, uncovered_e), v.Edges)
            
            if !isempty(e_candidates)
                e = e_candidates[1] # Take first valid edge
                push!(T, [v, e])
                deleteat!(uncovered_e, findall(x -> x == e, uncovered_e))
            end
        end
        deleteat!(uncovered_v, findall(x -> x == v, uncovered_v))
    end

    # Loop generation logic
    perms = Vector{Int64}[]
    for e in uncovered_e
        u = e.node1
        v = e.node2
        
        u_path = id_perm
        # Trace path for u
        ind = findall(i -> u == i[1], T)
        curr_u = u
        while !isempty(ind)
            eu = T[ind[1]][2]
            if curr_u == eu.node1
                u_path = p_compose(eu.correspondence12, u_path)
            else
                u_path = p_compose(eu.correspondence21, u_path) # sorted by default in tuple
            end
            curr_u = neighbor(curr_u, eu)
            ind = findall(i -> curr_u == i[1], T)
        end

        v_path = id_perm
        # Trace path for v
        ind = findall(i -> v == i[1], T)
        curr_v = v
        while !isempty(ind)
            ev = T[ind[1]][2]
            if curr_v == ev.node1
                v_path = p_compose(ev.correspondence12, v_path)
            else
                v_path = p_compose(ev.correspondence21, v_path)
            end
            curr_v = neighbor(curr_v, ev)
            ind = findall(i -> curr_v == i[1], T)
        end
        
        # Combine paths to form a cycle (permutation)
        # Perm = v_path * e_correspondence * inverse(u_path)
        final_perm_pairs = p_compose(v_path, p_compose(e.correspondence12, sort(map(i -> reverse(i), u_path))))
        push!(perms, map(i -> i[2], final_perm_pairs))
    end

    return perms
end
