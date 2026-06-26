# Results

`solve_monodromy` returns a `MonodromyResult` by default. It stores the mutated
vertices and edges, success/status metadata, iteration counts, and optional
diagnostics.

```julia
result = solve_monodromy(compiled, vertices; max_roots=4)

result.success
result.status
result.iterations
result.stagnant_iterations
length(result.vertices)
length(result.edges)
```

Edges store discovered correspondences in both directions.

```julia
e = result.edges[1]
e.correspondence12
e.correspondence21
```

The result also behaves like its `edges` vector for indexing and iteration.

```julia
length(result)
first(result)
```

```@docs
MonodromyResult
```
