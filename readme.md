# CertifiedHomotopyTracking.jl

**CertifiedHomotopyTracking.jl** is a Julia package for **certified homotopy tracking** and **monodromy group computation**. It uses interval arithmetic (via [Nemo.jl](https://github.com/Nemocas/Nemo.jl) and Arb) to provide mathematically rigorous results and integrates with [GAP](https://www.gap-system.org/) for group-theoretic analysis.

## Features

* **Certified tracking:** Uses interval arithmetic and the Krawczyk test to certify solution paths.
* **Certified monodromy computation:** Automatically tracks the complete homotopy graph to generate the monodromy group (`solve_monodromy`).


## Quick Start

### 1. Basic certified tracking

If you want to track a single solution path from $t=0$ to $t=1$:

```julia
using Nemo, AbstractAlgebra
using CertifiedMonodromyComputation

# 1. Set up the polynomial ring
@monodromy_setup begin
    vars = (x, y)
end
const CCi = _CCi # Alias for the coefficient ring (Complex Interval Field)

# 2. Define your system F(x, y) and the start system G(x, y)
f1 = x^2 + 3*y - 4
f2 = y^2 + 3

F = [f1 f2]
G = [x^2-1 y^2-1]

# 3. Define the start point at t=0 and the homotopy H
H = straight_line_homotopy(F, G, t)
point = [CCi(1), CCi(-1)]

# 4. Track!
track(H, point)
track(H, point; iterations_count=true) # print the number of iterations
track(H, point; show_display=false) # turn off the display

```


### 2. Certified monodromy group computation

To compute the monodromy group of a parameterized system:

```julia
using Nemo, AbstractAlgebra
using CertifiedMonodromyComputation

# 1. Set up the polynomial ring
@monodromy_setup begin
    vars = (x, y)
    params = (p, q)
end
const CCi = _CCi # Alias for the coefficient ring (Complex Interval Field)

# 2. Define your parameter system F(x, y; p, q)
f1 = p*x^2 + 3*y - 4
f2 = y^2 + q

F = [f1 f2]

# 3. Set up the initial seed
bp = [CCi(1), CCi(-1)] # values of p and q
x = [CCi(1) , CCi(1)] # a solution at (p, q) = bp


# 4. Set up a homotopy graph
v1 = vertex(bp,[x])
vs = parameter_points(v1, 2, 6) # make 6 parameter points (vertices) in C^2 including the vertex v1

# 5. Solve monodromy (Tracking)
edges = solve_monodromy(F, vs; max_roots=4) # it may require several tries to find all solutions for each vertex

# 6. GAP analysis
G = build_gap_group(4, edges) # Find a group of size 4 from edge correspondences

if G !== nothing
    println("Structure Description:")
    println(GAP.Globals.StructureDescription(G)) # C2 x C2
    println("Galois Width:")
    gw = galois_width(G)
    println(gw) # 2
end
```
