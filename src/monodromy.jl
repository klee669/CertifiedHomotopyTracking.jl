# ------------------------------------------------------------------------------
# High-Level Interface
# ------------------------------------------------------------------------------
# Exporting Vertex/Edge allows users to inspect results and build custom graphs
export solve_monodromy, build_gap_group, vertex, edge, galois_width,
    Edge, Vertex, MonodromyResult, build_edges

# ------------------------------------------------------------------------------
# Data Structures
# ------------------------------------------------------------------------------
const PointType = Vector{AcbFieldElem} 

"""
    Vertex

Node of a monodromy graph.

Fields:

- `base_point`: parameter value attached to the vertex.
- `sols`: known solutions over `base_point`.
- `Edges`: incident graph edges.

Construct vertices with [`vertex`](@ref).
"""
mutable struct Vertex
    base_point::PointType
    sols::Vector{PointType}
    Edges::Vector{Any} # Type Any to avoid circular dependency
end

"""
    Edge

Undirected monodromy graph edge with tracked correspondences in both directions.

Fields:

- `node1`, `node2`: endpoint vertices.
- `correspondence12`: pairs `(i, j)` mapping solution `i` at `node1` to solution
  `j` at `node2`.
- `correspondence21`: reverse correspondences.

Construct edges with [`edge`](@ref) or [`build_edges`](@ref).
"""
mutable struct Edge
    node1::Vertex
    node2::Vertex
    correspondence12::Vector{Tuple{Int, Int}}
    correspondence21::Vector{Tuple{Int, Int}}
end

"""
    MonodromyResult

Return object from [`solve_monodromy`](@ref).

Fields:

- `vertices`, `edges`: mutated graph containing discovered solutions and correspondences.
- `success`: whether every edge reached `max_roots` correspondences.
- `status`: `:success`, `:stagnant`, `:interrupted`, or `:partial`.
- `iterations`, `stagnant_iterations`: solve-loop counters.
- `max_roots`: target number of roots per edge.
- `diagnostics`: optional a posteriori path diagnostics.

`MonodromyResult` behaves like its `edges` vector for indexing and iteration.
"""
struct MonodromyResult
    vertices::Vector{Vertex}
    edges::Vector{Edge}
    success::Bool
    status::Symbol
    iterations::Int
    stagnant_iterations::Int
    max_roots::Int
    diagnostics::NamedTuple
end

# Constructors
"""
    vertex(p)
    vertex(p, solutions)

Create a monodromy graph vertex at parameter point `p`.

Use `vertex(p, [x0])` when one or more start solutions are already known.
"""
vertex(p::Vector) = Vertex(p, PointType[], [])

vertex(p::Vector, x::Vector) = Vertex(p, x, [])

"""
    edge(v, w)

Create an empty monodromy graph edge between vertices `v` and `w`.
"""
edge(p::Vertex, q::Vertex) = Edge(p, q, Tuple{Int,Int}[], Tuple{Int,Int}[])

Base.show(io::IO, v::Vertex) = print(io, "Vertex($(length(v.sols)) solutions)")
Base.show(io::IO, e::Edge) = print(io, "Edge($(length(e.correspondence12)) correspondences)")
Base.show(io::IO, result::MonodromyResult) = print(
    io,
    "MonodromyResult(status=$(result.status), success=$(result.success), edges=$(length(result.edges)), iterations=$(result.iterations))",
)
Base.length(result::MonodromyResult) = length(result.edges)
Base.isempty(result::MonodromyResult) = isempty(result.edges)
Base.getindex(result::MonodromyResult, i::Int) = result.edges[i]
Base.iterate(result::MonodromyResult, state...) = iterate(result.edges, state...)

"""
    build_edges(vertices, edge_pairs) -> Vector{Edge}

Create graph edges from pairs of vertex indices and register each edge with its
endpoint vertices.

# Example

```julia
using CertifiedHomotopyTracking;

CC = AcbField(128);
vertices = [vertex([CC(i)]) for i in 1:3];
edges = build_edges(vertices, [(1, 2), (2, 3), (3, 1)])
length(edges)
```
"""
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


"""
    make_edge_system(compiled, p_start, p_end, p_const=AcbFieldElem[]) -> SpecializedHomotopy

Specialize a [`CompiledHomotopy`](@ref), usually from
[`compile_edge_homotopy`](@ref), to the linear parameter path from `p_start` to
`p_end`.

`p_const` supplies fixed constants declared through `const_vars` or introduced
internally for complex numeric coefficients.

# Example

```julia
using CertifiedHomotopyTracking;

@variables x y p q;
CC = AcbField(256);
F = [p*x^2 + 3y - 4, y^2 + q];
compiled = compile_edge_homotopy(F, [x, y], [p, q]);
sys = make_edge_system(compiled, [CC(1), CC(-1)], [CC(cis(0.3)), CC(cis(0.7))]);
res = track_path(sys, [CC(1), CC(1)])
success(res)
```
"""
function make_edge_system(
    compiled_sys::CompiledHomotopy,
    p_start::Vector{AcbFieldElem},
    p_end::Vector{AcbFieldElem},
    p_const::Vector{AcbFieldElem}=AcbFieldElem[];
    patch_vector::Vector{AcbFieldElem}=AcbFieldElem[],
    source::Union{Nothing,HomotopySourceData}=compiled_sys.source,
)
    return SpecializedHomotopy(compiled_sys, p_start, p_end, p_const; patch_vector=patch_vector, source=source)
end

function _posteriori_endpoint(sys_edge::SpecializedHomotopy, cert)
    if haskey(cert, :affine_endpoint)
        return sys_edge.CC.(cert.affine_endpoint), true
    end
    if haskey(cert, :segments) && !isempty(cert.segments)
        _, idx = findmax(seg -> seg.t_end, cert.segments)
        seg = cert.segments[idx]
        haskey(seg, :x_end) && return copy(seg.x_end), true
    end
    if haskey(cert, :hc_trace) && cert.hc_trace !== nothing && !isempty(cert.hc_trace.trace)
        return sys_edge.CC.(last(cert.hc_trace.trace).x), true
    end
    return AcbFieldElem[], false
end

function _empty_monodromy_diagnostics()
    return (
        total_attempted_paths = 0,
        total_successful_paths = 0,
        total_failed_paths = 0,
        total_skipped_failed_paths = 0,
        failed_paths = NamedTuple[],
        skipped_failed_paths = NamedTuple[],
    )
end

function _posteriori_failure_summary(cert)
    cert === nothing && return (;)
    summary = (
        success = haskey(cert, :success) ? cert.success : missing,
        status = haskey(cert, :status) ? cert.status : missing,
        method = haskey(cert, :method) ? cert.method : missing,
        certification_chart = haskey(cert, :certification_chart) ? cert.certification_chart : missing,
        total_boxes = haskey(cert, :total_boxes) ? cert.total_boxes : missing,
        original_segments = haskey(cert, :original_segments) ? cert.original_segments : missing,
        max_depth = haskey(cert, :max_depth) ? cert.max_depth : missing,
        max_krawczyk_norm = haskey(cert, :max_krawczyk_norm) ? cert.max_krawczyk_norm : missing,
        max_step_size = haskey(cert, :max_step_size) ? cert.max_step_size : missing,
        failed_segments = haskey(cert, :failed_segments) ? length(cert.failed_segments) : missing,
    )
    return summary
end

function _record_monodromy_failure!(diagnostics, failure)
    diagnostics === nothing && return nothing
    push!(diagnostics.failed_paths, failure)
    return nothing
end

function _record_monodromy_skip!(diagnostics, skipped)
    diagnostics === nothing && return nothing
    push!(diagnostics.skipped_failed_paths, skipped)
    return nothing
end

function _freeze_monodromy_diagnostics(diagnostics)
    diagnostics === nothing && return _empty_monodromy_diagnostics()
    return (
        total_attempted_paths = diagnostics.total_attempted_paths[],
        total_successful_paths = diagnostics.total_successful_paths[],
        total_failed_paths = diagnostics.total_failed_paths[],
        total_skipped_failed_paths = diagnostics.total_skipped_failed_paths[],
        failed_paths = copy(diagnostics.failed_paths),
        skipped_failed_paths = copy(diagnostics.skipped_failed_paths),
    )
end

function _monodromy_diagnostics_state(enabled::Bool)
    enabled || return nothing
    return (
        total_attempted_paths = Ref(0),
        total_successful_paths = Ref(0),
        total_failed_paths = Ref(0),
        total_skipped_failed_paths = Ref(0),
        failed_paths = NamedTuple[],
        skipped_failed_paths = NamedTuple[],
    )
end

function _posteriori_diagnostics_requested(posteriori::Bool, posteriori_options::NamedTuple)
    posteriori || return false
    mode = get(posteriori_options, :diagnostics, :off)
    return mode !== :off && mode != false && mode !== nothing
end

function _track_compiled_edge_path(
    sys_edge::SpecializedHomotopy,
    start_point;
    posteriori::Bool,
    posteriori_options::NamedTuple,
    show_progress::Bool,
    track_options::NamedTuple,
)
    if posteriori
        options = haskey(posteriori_options, :show_progress) ?
            posteriori_options :
            merge((; show_progress = show_progress), posteriori_options)
        cert = certify_posteriori(
            sys_edge,
            ComplexF64.(start_point);
            options...,
        )
        !cert.success && return AcbFieldElem[], false, cert
        y_end, success = _posteriori_endpoint(sys_edge, cert)
        return y_end, success, cert
    end

    y_end, success = track_path(
        sys_edge,
        start_point;
        t_end = 1.0,
        h_init = 0.1,
        show_progress = show_progress,
        track_options...,
    )
    return y_end, success, nothing
end


# ------------------------------------------------------------------------------
# Core Logic: Monodromy Solving
# ------------------------------------------------------------------------------

"""
    solve_monodromy(compiled, vertices[, edges]; kwargs...) -> MonodromyResult

Track a monodromy graph and discover solution correspondences.

The current Symbolics-based workflow uses a [`CompiledHomotopy`](@ref) from
[`compile_edge_homotopy`](@ref). If `edges` are omitted, a complete graph is
built on `vertices`. If `edges` are supplied, the custom graph is tracked.

# Options for the compiled workflow

- `max_roots=20`: target number of correspondences per edge.
- `max_attempts=10`: stop after this many stagnant iterations.
- `show_progress=false`: print path-tracking progress.
- `root_match=:certified_or_heuristic`: root matching mode used when deciding whether a
  tracked endpoint is already present in the target vertex. Supported values are
  `:heuristic`, `:certified`, and `:certified_or_heuristic`.
- `projective=false`: track in projective charts. Requires `compiled` to be
  built with `projective=true`.
- `track_options=(;)`: keyword options forwarded to [`track_path`](@ref), for
  example `(; adaptive_precision=false, h_init=0.05)`.
- `posteriori=false`: additionally certify paths with
  [`certify_posteriori`](@ref).
- `posteriori_options=(;)`: options forwarded to a posteriori certification.
- `return_result=true`: return a [`MonodromyResult`](@ref); set to `false` for
  the legacy `Vector{Edge}` return value.
- `threading=false`, `ntasks=Threads.nthreads()`: track independent paths in
  parallel when enabled.

The graph is mutated in-place. Interrupting with Ctrl-C returns partial data.

# Example

```julia
using CertifiedHomotopyTracking;

@variables x y p q;
CC = AcbField(256);
F = [p*x^2 + 3y - 4, y^2 + q];
compiled = compile_edge_homotopy(F, [x, y], [p, q]);

v1 = vertex([CC(1), CC(-1)], [[CC(1), CC(1)]]);
vertices = [v1; [vertex([CC(cis(0.2k)), CC(cis(0.3k))]) for k in 1:3]];

# Omitting edges uses the complete graph on vertices.
result = solve_monodromy(compiled, vertices; max_roots=4)
length(result.edges)

# Or provide a custom graph explicitly.
vertices = [v1; [vertex([CC(cis(0.2k)), CC(cis(0.3k))]) for k in 1:5]];
edges = build_edges(vertices, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)])
length(edges)
```
"""
function solve_monodromy(
    H::Union{Matrix,Vector},
    vertices::Vector{Vertex};
    radius::Number = 0.1,
    max_roots::Int = 20,
    max_attempts::Int = 10,
    predictor::Bool = true,
    show_progress::Bool = true,
    return_result::Bool = true,
    threading::Bool = false,
    ntasks::Integer = Base.Threads.nthreads(),
)
    _thread_count(threading, ntasks)
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
        iter = 0
        status = :partial
        iter_stagnant = 0
        total_correspondences = 0
        
        while true
            iter += 1
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
                status = :success
                break
            end

            if iter_stagnant > max_attempts
                @warn "Stopping: No new solutions found after $max_attempts iterations."
                status = :stagnant
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
            @warn "Computation interrupted by user."
            @warn "Returning vertices and edges computed SO FAR."
            @warn "You can resume or analyze the partial data."
            result = MonodromyResult(
                vertices,
                edges,
                false,
                :interrupted,
                @isdefined(iter) ? iter : 0,
                @isdefined(iter_stagnant) ? iter_stagnant : 0,
                max_roots,
                _empty_monodromy_diagnostics(),
            )
            return return_result ? result : edges
        else
            # Re-throw other unexpected errors
            rethrow(e)
        end
    end

    result = MonodromyResult(
        vertices,
        edges,
        (@isdefined(status) ? status : :partial) === :success,
        @isdefined(status) ? status : :partial,
        @isdefined(iter) ? iter : 0,
        @isdefined(iter_stagnant) ? iter_stagnant : 0,
        max_roots,
        _empty_monodromy_diagnostics(),
    )
    return return_result ? result : edges
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
    root_match::Symbol=:certified_or_heuristic,
    track_options::NamedTuple=(;),
    posteriori::Bool=false,
    posteriori_options::NamedTuple=(;),
    failed_paths=nothing,
    diagnostics=nothing,
    threading::Bool=false,
    ntasks::Integer=Base.Threads.nthreads(),
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
    untracked_indices = Int[]
    for i in 1:length(source_sols)
        i in tracked_indices && continue
        key = (edge_idx, from1to2, i)
        if failed_paths !== nothing && key in failed_paths
            diagnostics !== nothing && (diagnostics.total_skipped_failed_paths[] += 1)
            _record_monodromy_skip!(
                diagnostics,
                (
                    edge_index = edge_idx,
                    direction = from1to2 ? :forward : :backward,
                    source_vertex = id_src,
                    target_vertex = id_dest,
                    source_index = i,
                ),
            )
            continue
        end
        push!(untracked_indices, i)
    end

    if isempty(untracked_indices)
        return
    end

    @info "Tracking Edge ($id_src -> $id_dest): $(length(untracked_indices)) new paths to track."

    sys_edge = make_edge_system(compiled_sys, source_v.base_point, target_v.base_point)
    
    count_success = 0; count_new = 0; count_collision = 0

    if _threading_enabled(threading, ntasks) && length(untracked_indices) > 1
        count_success, count_new, count_collision = _track_edge_paths_threaded!(
            compiled_sys,
            sys_edge,
            source_sols,
            target_sols,
            c_forward,
            c_backward,
            untracked_indices;
            edge_idx = edge_idx,
            from1to2 = from1to2,
            id_src = id_src,
            id_dest = id_dest,
            show_progress = show_progress,
            root_match = root_match,
            track_options = track_options,
            posteriori = posteriori,
            posteriori_options = posteriori_options,
            failed_paths = failed_paths,
            diagnostics = diagnostics,
            ntasks = _thread_count(threading, ntasks),
        )
        sort!(c_forward)
        sort!(c_backward)
        println(" Done. (Ok: $count_success, New: $count_new, Collision: $count_collision)")
        return
    end

    for src_idx in untracked_indices
        start_point = source_sols[src_idx]
        
        diagnostics !== nothing && (diagnostics.total_attempted_paths[] += 1)
        y_end, success, cert = _track_compiled_edge_path(
            sys_edge,
            start_point;
            posteriori = posteriori,
            posteriori_options = posteriori_options,
            show_progress = show_progress,
            track_options = track_options,
        )
        
        if success
            diagnostics !== nothing && (diagnostics.total_successful_paths[] += 1)
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
            diagnostics !== nothing && (diagnostics.total_failed_paths[] += 1)
            failed_paths !== nothing && push!(failed_paths, (edge_idx, from1to2, src_idx))
            _record_monodromy_failure!(
                diagnostics,
                merge(
                    (
                        edge_index = edge_idx,
                        direction = from1to2 ? :forward : :backward,
                        source_vertex = id_src,
                        target_vertex = id_dest,
                        source_index = src_idx,
                        reason = posteriori ? :posteriori_failed : :tracking_failed,
                    ),
                    _posteriori_failure_summary(cert),
                ),
            )
            printstyled("x", color=:red, bold=true)
        end
    end
    
    sort!(c_forward)
    sort!(c_backward)
    println(" Done. (Ok: $count_success, New: $count_new, Collision: $count_collision)")
end

function _track_edge_paths_threaded!(
    compiled_sys::CompiledHomotopy,
    sys_edge::SpecializedHomotopy,
    source_sols,
    target_sols,
    c_forward,
    c_backward,
    untracked_indices;
    edge_idx,
    from1to2,
    id_src,
    id_dest,
    show_progress,
    root_match,
    track_options,
    posteriori,
    posteriori_options,
    failed_paths,
    diagnostics,
    ntasks,
)
    results = Vector{Any}(undef, length(untracked_indices))
    chunks = _thread_chunks(untracked_indices, ntasks)
    Base.Threads.@sync for chunk in chunks
        Base.Threads.@spawn begin
            local_sys_edge = make_edge_system(
                compiled_sys,
                AcbFieldElem[sys_edge.p_start...],
                AcbFieldElem[sys_edge.p_end...],
                AcbFieldElem[sys_edge.p_const...];
                patch_vector = AcbFieldElem[sys_edge.patch_vector...],
                source = sys_edge.source,
            )
            for pos in chunk
                src_idx = untracked_indices[pos]
                y_end, success, cert = _track_compiled_edge_path(
                    local_sys_edge,
                    source_sols[src_idx];
                    posteriori = posteriori,
                    posteriori_options = posteriori_options,
                    show_progress = show_progress,
                    track_options = track_options,
                )
                results[pos] = (src_idx = src_idx, y_end = y_end, success = success, cert = cert)
            end
        end
    end

    count_success = 0
    count_new = 0
    count_collision = 0
    for result in results
        diagnostics !== nothing && (diagnostics.total_attempted_paths[] += 1)
        src_idx = result.src_idx
        if result.success
            diagnostics !== nothing && (diagnostics.total_successful_paths[] += 1)
            dest_idx = _search_solution(sys_edge, result.y_end, target_sols; root_match=root_match)
            if dest_idx === nothing
                push!(target_sols, result.y_end)
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
            diagnostics !== nothing && (diagnostics.total_failed_paths[] += 1)
            failed_paths !== nothing && push!(failed_paths, (edge_idx, from1to2, src_idx))
            _record_monodromy_failure!(
                diagnostics,
                merge(
                    (
                        edge_index = edge_idx,
                        direction = from1to2 ? :forward : :backward,
                        source_vertex = id_src,
                        target_vertex = id_dest,
                        source_index = src_idx,
                        reason = posteriori ? :posteriori_failed : :tracking_failed,
                    ),
                    _posteriori_failure_summary(result.cert),
                ),
            )
            printstyled("x", color=:red, bold=true)
        end
    end
    return count_success, count_new, count_collision
end



"""
    solve_monodromy(compiled_sys, vertices, edges; kwargs...)

Track the given graph for a compiled homotopy and return a
`MonodromyResult`. The result stores the mutated `vertices`, `edges`,
success/status metadata, iteration counts, and optional diagnostics. Set
`return_result = false` for the legacy `Vector{Edge}` return value.

Within one solve run, paths that already failed are cached and skipped on later
iterations with the same options.
"""
function solve_monodromy(
    compiled_sys::CompiledHomotopy,
    vertices::Vector{Vertex},
    edges::Vector{Edge};
    max_roots=20,
    max_attempts::Int=10,
    show_progress::Bool=false,
    root_match::Symbol=:certified_or_heuristic,
    projective::Bool=false,
    track_options::NamedTuple=(;),
    posteriori::Bool=false,
    posteriori_options::NamedTuple=(;),
    return_result::Bool=true,
    threading::Bool=false,
    ntasks::Integer=Base.Threads.nthreads(),
)
    if projective && !compiled_sys.projective_patch
        throw(ArgumentError("projective=true requires compiled_sys to be built with projective=true."))
    end
    _thread_count(threading, ntasks)
    effective_track_options = merge((; projective=projective), track_options)
    iter = 0
    iter_stagnant = 0
    total_correspondences = 0
    status = :partial
    failed_paths = Set{Tuple{Int,Bool,Int}}()
    diagnostics_state = _monodromy_diagnostics_state(
        _posteriori_diagnostics_requested(posteriori, posteriori_options),
    )
    
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
                println("Success! All edges reached $max_roots correspondences.")
                status = :success
                break
            end

            if iter_stagnant > max_attempts
                println("Stopping: No new solutions found after $max_attempts iterations.")
                status = :stagnant
                break
            end
            
            for (idx, e) in enumerate(edges)
                id1 = findfirst(==(e.node1), vertices)
                id2 = findfirst(==(e.node2), vertices)
                
                track_edge!(
                    compiled_sys,
                    e,
                    true,
                    idx,
                    id1,
                    id2;
                    show_progress = show_progress,
                    root_match = root_match,
                    track_options = effective_track_options,
                    posteriori = posteriori,
                    posteriori_options = posteriori_options,
                    failed_paths = failed_paths,
                    diagnostics = diagnostics_state,
                    threading = threading,
                    ntasks = ntasks,
                )
                track_edge!(
                    compiled_sys,
                    e,
                    false,
                    idx,
                    id2,
                    id1;
                    show_progress = show_progress,
                    root_match = root_match,
                    track_options = effective_track_options,
                    posteriori = posteriori,
                    posteriori_options = posteriori_options,
                    failed_paths = failed_paths,
                    diagnostics = diagnostics_state,
                    threading = threading,
                    ntasks = ntasks,
                )

                counts = map(x -> length(x.correspondence12), edges)
                @info "Status (Edge $idx done): Correspondences => $counts"
            end
        end
        
    catch e
        if isa(e, InterruptException)
            println("\n\nComputation interrupted by user.")
            println("Returning vertices and edges computed SO FAR.")
            status = :interrupted
            result = MonodromyResult(
                vertices,
                edges,
                false,
                status,
                iter,
                iter_stagnant,
                max_roots,
                _freeze_monodromy_diagnostics(diagnostics_state),
            )
            return return_result ? result : edges
        else
            rethrow(e)
        end
    end
    
    success = status === :success
    result = MonodromyResult(
        vertices,
        edges,
        success,
        status,
        iter,
        iter_stagnant,
        max_roots,
        _freeze_monodromy_diagnostics(diagnostics_state),
    )
    return return_result ? result : edges
end
function solve_monodromy(
    compiled_sys::CompiledHomotopy,
    vertices::Vector{Vertex};
    max_roots=20,
    max_attempts::Int=10,
    show_progress::Bool=false,
    root_match::Symbol=:certified_or_heuristic,
    projective::Bool=false,
    track_options::NamedTuple=(;),
    posteriori::Bool=false,
    posteriori_options::NamedTuple=(;),
    return_result::Bool=true,
    threading::Bool=false,
    ntasks::Integer=Base.Threads.nthreads(),
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
    
    return solve_monodromy(
        compiled_sys,
        vertices,
        edges;
        max_roots = max_roots,
        max_attempts = max_attempts,
        show_progress = show_progress,
        root_match = root_match,
        projective = projective,
        track_options = track_options,
        posteriori = posteriori,
        posteriori_options = posteriori_options,
        return_result = return_result,
        threading = threading,
        ntasks = ntasks,
    )
end

# ------------------------------------------------------------------------------
# GAP Integration
# ------------------------------------------------------------------------------

"""
    build_gap_group(max_roots, edges_or_result)

Construct a GAP permutation group from monodromy edge correspondences.

Returns `nothing` if no valid permutations are available.
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

build_gap_group(max_roots::Int, result::MonodromyResult) =
    build_gap_group(max_roots, result.edges)

"""
    galois_width(G::GAP.GapObj)

Calculate the Galois width of a GAP group using a GAP helper function.
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

"""
    search_point(res, pool; tol=1e-4)

Heuristically find `res` in `pool` by comparing the max norm of interval-ball
differences to `tol`. Returns the matching index or `nothing`.
"""
function search_point(res::Vector{AcbFieldElem}, pool::Vector{Vector{AcbFieldElem}}; tol=1e-4)
    for (idx, sol) in enumerate(pool)
        dist = maximum(mag_complex.(res .- sol))
        if dist < tol
            return idx
        end
    end
    return nothing
end

function _tracking_coordinates(sys::SpecializedHomotopy, x::Vector{AcbFieldElem})
    if uses_projective_charts(sys)
        if length(x) == sys.compiled.n_vars - 1
            return copy(x)
        elseif length(x) == sys.compiled.n_vars
            return _projective_to_chart_coordinates(sys, x, sys.patch_idx)
        end
    elseif has_projective_patch(sys) && length(x) == length(sys.patch_vector) - 1
        return lift_to_patch(x, collect(sys.patch_vector))
    elseif sys.projective_coordinates && !has_projective_patch(sys) && sys.compiled.n_vars != 0 && length(x) == sys.compiled.n_vars - 1
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

function _krawczyk_same_on_hull(sys::SpecializedHomotopy, x::Vector{AcbFieldElem}, y::Vector{AcbFieldElem}, t; rho, inflate, min_radius)
    center, radius = _enclosing_center_radius(x, y; inflate=inflate, min_radius=min_radius)
    A = compute_preconditioner(sys, center, t)
    passed, k_norm = krawczyk_test(sys, center, t, radius, A; rho=rho)
    return passed, k_norm, radius
end

function _krawczyk_same_at_center(sys::SpecializedHomotopy, center::Vector{AcbFieldElem}, x::Vector{AcbFieldElem}, y::Vector{AcbFieldElem}, t; rho, inflate, min_radius)
    radius = _enclosing_radius_about(center, x, y; inflate=inflate, min_radius=min_radius)
    A = compute_preconditioner(sys, center, t)
    passed, k_norm = krawczyk_test(sys, center, t, radius, A; rho=rho)
    return passed, k_norm, radius
end

function _polish_for_root_comparison(sys::SpecializedHomotopy, x::Vector{AcbFieldElem}, t; radius=1e-8)
    center = get_mid_vec(x)
    A = compute_preconditioner(sys, center, t)
    polished, _, _, success = refine_moore_box(sys, center, t, radius, A)
    return success ? get_mid_vec(polished) : center, success
end

"""
    same_root_krawczyk(sys, x, y; t=1.0, rho=0.7, inflate=1.1, min_radius=1e-12)

Try to certify that `x` and `y` represent the same root of `sys` at parameter
value `t`.

Returns `(status, k_norm, radius)` where `status` is `:same` or `:unknown`.
"""
function same_root_krawczyk(
    sys::SpecializedHomotopy,
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

"""
    search_point_certified(sys, res, pool; t=1.0, rho=0.7, inflate=1.1, min_radius=1e-12)

Search for `res` in `pool` using [`same_root_krawczyk`](@ref). Returns the
matching index or `nothing`.

This is the pool-search wrapper around [`same_root_krawczyk`](@ref): it tests
`res` against each candidate and returns the first certified match.
"""
function search_point_certified(
    sys::SpecializedHomotopy,
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

function _search_solution(sys::SpecializedHomotopy, res::Vector{AcbFieldElem}, pool::Vector{Vector{AcbFieldElem}}; root_match::Symbol=:certified_or_heuristic)
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

complete_correspondences(rc::Int, result::MonodromyResult) =
    complete_correspondences(rc, result.edges)

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

get_permutations(rc::Number, result::MonodromyResult) =
    get_permutations(rc, result.edges)
