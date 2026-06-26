# Path Tracking

The certified tracking implementation follows Krawczyk-based homotopy tracking
ideas from [duff2024certified,guillemot2024validated](@cite).

## Basic Usage

```@setup path_tracking_example
using CertifiedHomotopyTracking
```

```@repl path_tracking_example
@variables x y;
CC = AcbField(256);
F = [x^2 + 3y - 4, y^2 + 3];
G = [x^2 - 1, y^2 - 1];
H = straight_line_homotopy(F, G, [x, y]; CCRing=CC, gamma=CC(0.5, 0.5));

res = track_path(H, [CC(1), CC(-1)]; show_progress=false)
success(res)
solution(res)
evaluate_H(H, certified_region(res), CC(1))
```

```@docs
track_path
```

```@raw html
<object data="../assets/path2.pdf" type="application/pdf" width="100%" height="520">
  <a href="../assets/path2.pdf">Open example path PDF</a>
</object>
```

## References

```@bibliography
Pages = [@__FILE__]
```
