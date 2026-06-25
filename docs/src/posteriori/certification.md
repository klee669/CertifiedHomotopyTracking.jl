# Certification

The high-level entry point is [`certify_posteriori`](@ref). It first
collects a numerical HomotopyContinuation.jl trace and then certifies the trace
with CHT interval checks. By default, each trace segment is certified in an
adaptively chosen local parameter: CHT may use `t`, `Re(x_i)`, or `Im(x_i)`
instead of forcing the whole argument to use `t`.

```@repl posteriori_certification
using CertifiedHomotopyTracking;

@variables x y;
CC = AcbField(128);
F = [x^2 + 3*y - 4, y^2 + 3];
G = [x^2 - 1, y^2 - 1];

H = straight_line_homotopy(F, G, [x, y]; CCRing=CC, gamma=CC(0.5, 0.5));
cert = certify_posteriori(H, [CC(1), CC(-1)]; max_step_size=Inf, max_depth=6)

cert.success
success(cert)
solution(cert)
certified_region(cert)
cert.total_boxes
cert.failed_segments
```

```@docs
certify_posteriori
PosterioriPathResult
solution(::PosterioriPathResult)
certified_region(::PosterioriPathResult)
```

## Threading

Start Julia with multiple threads before enabling threaded certification:

```bash
julia --threads=4 --project
```

or:

```bash
JULIA_NUM_THREADS=4 julia --project
```

Then pass `threading=true`. `ntasks` controls the maximum number of concurrent
segment-certification tasks used by CHT.

```julia
cert = certify_posteriori(
    H,
    [CC(1), CC(-1)];
    threading=true,
    ntasks=4,
)
```

## Diagnostics

Set `diagnostics=:basic` to record certification counters and local-parameter
choices. Use `:timing` to include timing fields as well.

```@repl posteriori_certification
cert = certify_posteriori(
    H,
    [CC(1), CC(-1)];
    max_step_size=Inf,
    max_depth=6,
    diagnostics=:basic,
);

cert.diagnostics
cert.diagnostics.local_parameter_choices
```

With `diagnostics=:timing`, the same summary includes timing fields such as
`time_in_refinement`, `time_in_validation`, and
`time_in_local_parameter_search`.
