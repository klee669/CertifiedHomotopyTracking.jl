# Krawczyk Test

These lower-level helpers expose the interval Krawczyk tests used by the path
tracker, Moore-box refinement, and variety certification routines.

## `krawczyk_test`

```@docs
krawczyk_test
```

```@repl krawczyk_test_example
using CertifiedHomotopyTracking;

@variables x y;
CC = AcbField(128);
compiled = compile_system([x^2 + y^2 - 2, x^3 - y], [x, y]);
sys = SpecializedHomotopy(compiled, CC);

x0 = [CC(1.001), CC(0.999)];
t = CC(0);
A = compute_preconditioner(sys, x0, t);

krawczyk_test(sys, x0, t, 1e-2, A)
krawczyk_test(sys, x0, t, 1e-3, A)
```

## `compute_preconditioner`

```@docs
compute_preconditioner
```

## `refine_moore_box`

```@docs
refine_moore_box
```

```@repl krawczyk_test_example
refined, radius, A_refined, success = refine_moore_box(sys, x0, t, 1e-4, A);
success
refined
radius
krawczyk_test(sys, refined, t, radius, A_refined)
```

## `krawczyk_operator`

```@docs
krawczyk_operator
```
