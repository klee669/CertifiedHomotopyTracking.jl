# Tracking Results

## Result Accessors

```@setup tracking_results_example
using CertifiedHomotopyTracking
```

```@repl tracking_results_example
@variables x y;
CC = AcbField(256);
F = [x^2 + 3y - 4, y^2 + 3];
G = [x^2 - 1, y^2 - 1];
H = straight_line_homotopy(F, G, [x, y]; CCRing=CC);
res = track_path(H, [CC(1), CC(-1)]);

success(res)
solution(res)
certified_region(res)
approximate_solution(res; digits=5)
```

```@docs
TrackResult
success
solution
certified_region
projective_solution
near_infinity
input_start
refined_start
projective_input_start
projective_refined_start
approximate_solution
```
