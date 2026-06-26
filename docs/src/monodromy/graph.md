# Homotopy Graph

Monodromy tracking uses a graph of parameter points. A [`Vertex`](@ref) stores a
parameter value and the solutions currently known over that value. An [`Edge`](@ref)
stores the correspondences found by tracking paths between two vertices.

The high-level solver can build a complete graph automatically, so most examples
do not need to construct edges by hand. Build explicit edges only when you want
to track a chosen graph.

## Vertices And Edges

```@setup monodromy_graph
using CertifiedHomotopyTracking
```

```@repl monodromy_graph
CC = AcbField(256);
v1 = vertex([CC(1), CC(-1)], [[CC(1), CC(1)]])
v2 = vertex([CC(cis(0.2)), CC(cis(0.3))])
e = edge(v1, v2)
e.correspondence12
```

## Custom Graph

Pass an explicit edge list to [`solve_monodromy`](@ref) when the complete graph
is not the graph you want to track. The following example builds a cyclic graph
on six vertices.

```@setup monodromy_graph_custom
using CertifiedHomotopyTracking
```

```@repl monodromy_graph_custom
@variables x y p q;
CC = AcbField(256);
F = [p*x^2 + 3y - 4, y^2 + q];
compiled = compile_edge_homotopy(F, [x, y], [p, q]);

v1 = vertex([CC(1), CC(-1)], [[CC(1), CC(1)]]);
vertices = [v1; [vertex([CC(cis(0.2k)), CC(cis(0.3k))]) for k in 1:5]];
edges = build_edges(vertices, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]);
length(edges)
```

```julia
result = solve_monodromy(compiled, vertices, edges; max_roots=4)
length(result.edges)
```

```@docs
Vertex
Edge
vertex
edge
build_edges
```
