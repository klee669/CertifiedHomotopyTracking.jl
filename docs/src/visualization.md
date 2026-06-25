# [Visualization](@id visualization)

```@repl visualization_api_example
using CertifiedHomotopyTracking;

@variables x y;
CC = AcbField(128);
F = [x^2 + y - 2, y^2 + x - 2];
G = [x^2 - 1, y^2 - 1];
H = straight_line_homotopy(F, G, [x, y]; CCRing=CC, gamma=CC(1, 1));
res = track_path(H, [CC(1), CC(1)]; visualize=true);

length(path_boxes(res))
filename = tempname() * ".tex";
export_path_tikz(res, filename; axes=(:t, 1)) == filename
```

```@docs
PathBox
PathVisualization
path_boxes
export_path_tikz
export_path_obj
```
