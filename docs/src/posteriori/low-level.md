# Numerical Trace

The following lower-level calls expose the intermediate HomotopyContinuation.jl
objects and numerical trace used by [`certify_posteriori`](@ref). This is
mostly useful for debugging or inspecting what is being certified.

```@setup posteriori_low_level
using CertifiedHomotopyTracking
```

```@repl posteriori_low_level
@variables x y;
CC = AcbField(128);
F = [x^2 + 3*y - 4, y^2 + 3];
G = [x^2 - 1, y^2 - 1];
H = straight_line_homotopy(F, G, [x, y]; CCRing=CC, gamma=CC(0.5, 0.5));

tracker = prepare_posteriori_tracker(H);
H_hc = posteriori_hc_homotopy(tracker);
trace = collect_hc_trace(H_hc, ComplexF64[1, -1]; t_start=0.0, t_target=1.0);

trace.success
trace.status
trace.accepted_steps
trace.rejected_steps
length(trace.trace)
trace.trace[1]
trace.trace[end]
```

```@docs
PosterioriTracker
prepare_posteriori_tracker
posteriori_hc_system
posteriori_hc_homotopy
collect_hc_trace
```
