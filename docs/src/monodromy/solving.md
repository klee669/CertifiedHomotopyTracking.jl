# Solving Monodromy

Build a compiled parameter homotopy, seed one vertex with a known solution, and
run monodromy tracking. The basic workflow only needs vertices: if no explicit
edge list is supplied, [`solve_monodromy`](@ref) tracks the complete graph on
those vertices. The graph-based monodromy workflow follows the computational
monodromy framework of [duff2019solving](@cite), with certified homotopy-graph
ideas from [duff2026certifying](@cite).

```@setup monodromy_solving
using CertifiedHomotopyTracking
```

```@repl monodromy_solving
@variables x y p q;
CC = AcbField(256);
F = [p*x^2 + 3y - 4, y^2 + q];
compiled = compile_edge_homotopy(F, [x, y], [p, q]);

v1 = vertex([CC(1), CC(-1)], [[CC(1), CC(1)]]);
vertices = [v1; [vertex([CC(cis(0.2k)), CC(cis(0.3k))]) for k in 1:3]];
length(vertices)
```

Solve the complete monodromy graph by passing only `compiled` and `vertices`.

```julia
result = solve_monodromy(compiled, vertices; max_roots=4)
length(result.edges)
result.success
```

For threaded monodromy tracking, pass `threading=true` and choose the number of
tasks.

```julia
result = solve_monodromy(
    compiled,
    vertices;
    threading=true,
    ntasks=4,
)
```

Use [Homotopy Graph](graph.md) when you want to inspect vertices and edges as
data structures or supply a custom graph instead of the complete graph.

```@docs
solve_monodromy
```

## References

```@bibliography
Pages = [@__FILE__]
```
