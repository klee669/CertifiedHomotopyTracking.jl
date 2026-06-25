# Homotopies and Systems

```@docs
straight_line_homotopy
compile_homotopy
compile_system
compile_edge_homotopy
make_edge_system
SpecializedHomotopy
CompiledHomotopy
HomotopySourceData
evaluate_H
evaluate_Jac
evaluate_dt
system_with_precision
```

## Evaluation Example

```@repl homotopy_evaluation
using CertifiedHomotopyTracking;

@variables x y;
CC = AcbField(128);
F = [x^2 + y - 2, y^2 + x - 2];
G = [x^2 - 1, y^2 - 1];

H = straight_line_homotopy(F, G, [x, y]; CCRing=CC, gamma=CC(1, 1));
point = [CC(1), CC(1)];
t = CC(0.25);

evaluate_H(H, point, t)
evaluate_Jac(H, point, t)
evaluate_dt(H, point, t)
```

## Direct Homotopy Example

```@repl direct_homotopy
using CertifiedHomotopyTracking;

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
