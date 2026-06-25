# Solving

Build a compiled edge homotopy, seed one vertex with a known solution, and run
monodromy tracking.

```@repl monodromy_solving
using CertifiedHomotopyTracking;

@variables x y p q;
CC = AcbField(256);
F = [p*x^2 + 3y - 4, y^2 + q];
compiled = compile_edge_homotopy(F, [x, y], [p, q]);

v1 = vertex([CC(1), CC(-1)], [[CC(1), CC(1)]]);
vertices = [v1; [vertex([CC(cis(0.2k)), CC(cis(0.3k))]) for k in 1:3]];
length(vertices)
```

Omitting `edges` makes `solve_monodromy` build and track the complete graph on
the given vertices.

```julia
result = solve_monodromy(compiled, vertices; max_roots=4)
length(result.edges)
```

To track a custom graph instead of the complete graph, build the edges
explicitly.

```@repl monodromy_custom_graph
using CertifiedHomotopyTracking;

@variables x y p q;
CC = AcbField(256);
F = [p*x^2 + 3y - 4, y^2 + q];
compiled = compile_edge_homotopy(F, [x, y], [p, q]);

v1 = vertex([CC(1), CC(-1)], [[CC(1), CC(1)]]);
vertices = [v1; [vertex([CC(cis(0.2k)), CC(cis(0.3k))]) for k in 1:5]];
edges = build_edges(vertices, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]);
length(edges)
```

For threaded monodromy tracking, pass `threading=true` and choose the number of
tasks.

```julia
result = solve_monodromy(
    compiled,
    vertices,
    edges;
    threading=true,
    ntasks=4,
)
```

```@docs
solve_monodromy
vertex
edge
build_edges
```
