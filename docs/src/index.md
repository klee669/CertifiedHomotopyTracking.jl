# CertifiedHomotopyTracking.jl

CertifiedHomotopyTracking.jl provides interval-certified tools for tracking homotopy paths 
and relevant computations (e.g., computing monodromy actions and certified local approximations of algebraic varieties). It pairs Symbolics/Nemo-style model construction with ACB interval
arithmetic, with interactions to HomotopyContinuation.jl and GAP.

Main functionalities are like below:

- [A Posteriori Certification](posteriori/certification.md)
- [Tracking](tracking/path-tracking.md)
- [Monodromy](monodromy/solving.md)
- [Variety Approximations](variety-approximations.md)

There are relevant constructions and helpers for each task:

- [Homotopies and Systems](homotopies-and-systems.md)
- [Krawczyk Test](krawczyk-test.md)
- [Visualization](visualization.md)

## Quick Start

```@setup quick_start
using CertifiedHomotopyTracking
```

```@repl quick_start
@variables x y t;
gamma = 1 + im;
compiled = compile_homotopy([(1 - t) * gamma * (x^3 - 1) + t * (x^3 + 2y - 2),
                             (1 - t) * gamma * (y^2 - 1) + t * (y^2 + x - 2)],
                            [x, y], t);
CC = AcbField(256);
res = track_path(compiled, [CC(1), CC(1)])
success(res)
solution(res)
```

See the sidebar for the currently documented public interface.
