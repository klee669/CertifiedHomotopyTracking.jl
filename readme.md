# CertifiedHomotopyTracking.jl
![alt text](src/CHTsticker_page-0001.jpg "Logo Title Text 1")

**CertifiedHomotopyTracking.jl** is a Julia package for **certified homotopy tracking** and **monodromy group computation**. It uses interval arithmetic (via [Nemo.jl](https://github.com/Nemocas/Nemo.jl) and Arb) to provide mathematically rigorous results and integrates with [GAP](https://www.gap-system.org/) for group-theoretic analysis.


## Table of Contents
* [Features](#features)
* [Quick Start](#quick-start)
  * [1. Basic certified tracking](#1-basic-certified-tracking)
  * [2. Certified monodromy group computation](#2-certified-monodromy-group-computation)


## Features

* **Certified tracking:** Uses interval arithmetic and the Krawczyk test to certify solution paths.
* **Certified monodromy computation:** Automatically tracks the complete homotopy graph to generate the monodromy group (`solve_monodromy`).


## Quick Start

### 1. Basic certified tracking

If you want to track a single solution path from $t=0$ to $t=1$:

```julia
using CertifiedHomotopyTracking

# 1. Set up the polynomial ring
@variables x y
const PREC_BITS = 256
const CC = AcbField(PREC_BITS) # Complex Field (acb)


# 2. Define your system F(x, y) and the start system G(x, y)
f1 = x^2 + 3*y - 4
f2 = y^2 + 3

F = [f1, f2]
G = [x^2-1, y^2-1]

# 3. Define the start point at t=0 and the homotopy H
H = straight_line_homotopy(F, G, [x, y])
point = [CC(1), CC(-1)]

# 4. Track!
y, res_boolean = track_path(H, point)
max_norm(hcat(evaluate_H(H, y, CC(1)))) # check the residual of the result!

```


### 2. Certified monodromy group computation

To compute the monodromy group of a parameterized system:

```julia
using CertifiedHomotopyTracking

# 1. Set up the variables
@variables x y
@variables p q
const PREC_BITS = 256
const CC = AcbField(PREC_BITS) # Complex Field (acb)


# 2. Define your parameter system F(x, y; p, q)
f1 = p*x^2 + 3*y - 4
f2 = y^2 + q
F_exprs = [
    f1, f2
]

x_vars = [x,y]
p_vars = [p,q]

# 3. Set up the initial seed
bp = [CC(1), CC(-1)] # values of p and q
x0 = [CC(1) , CC(1)] # a solution at (p, q) = bp


# 4. Set up a homotopy graph
v1 = vertex(bp,[x0])
vertices = [v1]

for i in 1:3
    rand_u = [CC(cis(rand()*2*pi)) for _ in 1:2]
    push!(vertices, vertex(rand_u))
end

compiled_homotopy = compile_edge_homotopy(F_exprs, x_vars, p_vars; homogeneous=false);

# 5. Solve monodromy (Tracking)
edges = solve_monodromy(compiled_homotopy, vertices; max_roots=4)

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


