# Visualization

Set `visualize` to keep the certified boxes, and export them to TikZ with either
two or three axes. The following example draws the certification boxes in
`(t, Re(x1), Im(x2))` coordinates and overlays the numerical HC.jl trace.

```@setup posteriori_visualization
using CertifiedHomotopyTracking
```

```@repl posteriori_visualization
@variables x y;
CC = AcbField(128);
F = [x^2 + 3*y - 4, y^2 + 3];
G = [x^2 - 1, y^2 - 1];
H = straight_line_homotopy(F, G, [x, y]; CCRing=CC, gamma=CC(0.5, 0.5));

cert = certify_posteriori(
    H,
    [CC(1), CC(-1)];
    max_step_size=Inf,
    max_depth=6,
    visualize=(axes=(:t, (1, :real), (2, :imag)), show_trace=true),
)

cert
cert.success
cert.total_boxes
length(path_boxes(cert))
export_path_tikz(
    cert,
    "posteriori_visualization.tex";
    axes=(:t, (1, :real), (2, :imag)),
    show_trace=true,
)
```

```@raw html
<object data="../assets/posteriori_visualization.pdf" type="application/pdf" width="100%" height="420">
  <a href="../assets/posteriori_visualization.pdf">Open posteriori visualization PDF</a>
</object>
```
