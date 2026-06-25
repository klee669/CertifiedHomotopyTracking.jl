# GAP Computations

After solving enough correspondences, convert the monodromy data to a GAP
permutation group and ask GAP for group-theoretic information.

```julia
result = solve_monodromy(compiled, vertices; max_roots=4)

G = build_gap_group(4, result)
GAP.Globals.StructureDescription(G)
galois_width(G)
```

`build_gap_group` returns `nothing` if the available correspondences do not
define any valid permutations.

```julia
G === nothing
```

```@docs
build_gap_group
galois_width
```
