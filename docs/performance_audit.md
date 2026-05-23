# Performance Audit

이 문서는 CertifiedHomotopyTracking.jl의 현재 추적 구현을 성능 관점에서 읽은 정적 감사 결과이다. 이번 변경은 문서 추가만 포함하며, 수학적 판정식이나 public API는 변경하지 않는다.

## Scope

- 감사 대상: `src/tracking.jl`, `src/internals/*.jl`, `src/homotopy.jl`, `src/monodromy.jl`, `examples/irreducible_surface_optimization/surface_optimization.jl`.
- 비교 대상: 첨부된 `algpath-main-2.zip`의 `src/traits/track/*`, `src/prog.rs`, `src/taylor_model/*`, `src/reckless/*`, `src/flint/*`.
- 기준: 인증 조건은 유지하고, 최적화 전 benchmark를 먼저 추가하며, 각 최적화 PR은 작고 되돌리기 쉬워야 한다.

## Current Structure

핵심 파일은 다음처럼 나뉜다.

- `src/tracking.jl`: public `track`, `track_path`, `track_path_adaptive_precision`와 메인 추적 루프.
- `src/internals/moore_box.jl`: Moore box refinement.
- `src/internals/krawczyk.jl`: Krawczyk validation, preconditioner, Taylor-model step validation.
- `src/internals/predictors.jl`: tangent velocity, linear predictor, Hermite predictor, `TaylorModel3` predictor construction.
- `src/internals/systems.jl`: `CompiledHomotopy`, `HCSystem`, compiled `evaluate_H`, `evaluate_Jac`, `evaluate_dt`.
- `src/internals/homotopy_constructor.jl`: `Symbolics.build_function` 기반 compiled homotopy construction.
- `src/internals/taylor_model.jl`: order-3 Taylor model arithmetic.
- `src/internals/linear_algebra.jl`: symbolic `jac`, symbolic evaluation, Acb matrix inverse helpers.
- `src/internals/interval_arithmetic.jl`: Acb/Arb utility layer.
- `src/monodromy.jl`: monodromy edge loop에서 `track_path` 또는 legacy `track` 호출.

현재는 두 tracking path가 공존한다.

- Legacy symbolic path: `track(H, x)`와 `tracking_without_predictor(H, x)`.
- Compiled path: `track_path(sys::HCSystem, x_start_input)`, `track_path_adaptive_precision(sys, x_start_input)`.

성능 개선의 주된 대상은 compiled path가 되어야 한다. Legacy path는 public API 보존을 위해 남기되, 우선순위는 benchmark와 regression coverage로 묶어 관리하는 편이 안전하다.

## Tracking Call Graph

### Compiled `HCSystem` Path

```text
track_path_adaptive_precision(sys, x)
  -> _normalized_precision_schedule(sys, precisions)
  -> _convert_system_precision(sys, bits)
  -> track_path(trial_sys, trial_x; kwargs...)

track_path(sys, x_start_input)
  -> optional projective setup
       has_projective_patch(sys)
       lift_to_patch(...) / repatch(...)
       homogeneous chart normalization
  -> compute_preconditioner(sys, x, t)
       -> get_mid_vec(x)
       -> evaluate_Jac(sys, x_mid, t_mid)
            -> sys.compiled.func_Jx(...)
            -> optional patch row allocation
       -> inv_acb(J_val, sys.CC)
  -> refine_moore_box(sys, x, t, r, A)
       -> krawczyk_test(sys, y, t, s, U; rho=tau)
            -> evaluate_H(sys, y, t)
            -> evaluate_Jac(sys, y .+ B*r, t)
            -> norm_inf(K)
       -> evaluate_H(sys, y, t)
       -> delta = U * fy
       -> get_mid_vec(y - delta)
       -> evaluate_Jac(sys, y, t)
       -> inv_acb(Jy, sys.CC)
  -> compute_velocity(sys, x, t, A)
       -> evaluate_dt(sys, get_mid_vec(x), t_mid)
       -> -A * Ht
  -> while t < t_target
       -> optional homogeneous chart switch
       -> refine_moore_box(sys, x, t, r, A)
       -> compute_velocity(sys, x, t, A)
       -> while h > min_h
            -> first step:
                 TaylorModel3(x[i], v[i], 0, 0, 0, h)
               later steps:
                 construct_hermite_predictor_tm(sys, x, x_prev, v, v_prev, h_prev, h)
            -> validate_step_taylor3(sys, X_tm, t, h, r, A; rho)
                 -> t_tm = TaylorModel3(...)
                 -> evaluate_H(sys, X_tm, t_tm)
                 -> evaluate_taylor!(...) for F and X
                 -> evaluate_Jac(sys, X_bound .+ B*r, T_expanded)
                 -> norm_inf(K)
            -> if accepted:
                 evaluate_taylor!(...) for X_tm
                 get_mid_vec(x_next_interval)
                 h = min(2h, 0.5)
               else:
                 h /= 2
  -> final compute_preconditioner(...)
  -> final refine_moore_box(...)
  -> _track_result(...)
```

메인 loop는 `src/tracking.jl:304`부터이고, step-size rejection loop는 `src/tracking.jl:349`부터이다. Refinement는 `src/internals/moore_box.jl:31`, Krawczyk validation은 `src/internals/krawczyk.jl:58`과 `src/internals/krawczyk.jl:84`가 중심이다.

### Legacy Symbolic Path

```text
track(H, x)
  -> evaluate_matrix(Matrix(transpose(hcat(H))), 0)
  -> jacobian_inverse(G, x)
       -> jac(system)
       -> evaluate_matrix(J, x)
       -> pseudo_inv(...) or inv(...)
  -> linear_tracking(...)
       -> refine_step(...)
            -> evaluate_matrix(...)
            -> refine_moore_box(symbolic system, ...)
            -> speed_vector(...)
       -> linear_predictor(...)
       -> krawczyk_operator_taylor_model(...)
       -> proceeding_step(...)
  -> while t < 1
       -> hermite_tracking(...)
            -> refine_step(...)
            -> hermite_predictor(...)
            -> krawczyk_operator_taylor_model_original(...)
            -> proceeding_step(...)
       -> evaluate_matrix(matrix(X), input)
       -> midpoint_complex_box(...)
  -> final refine_step(...)
```

Legacy path는 매 step마다 symbolic `jac`, `derivative`, `AbstractAlgebra.evaluate`, `hcat`, `transpose`, `Matrix` 변환이 반복된다. 인증 연구용 구현으로는 읽기 좋지만, Algpath-style 성능 개선의 출발점으로 삼기에는 allocation과 dispatch 비용이 크다.

## Main Components

### Main Tracking Loop

- `track_path` is the active compiled tracker.
- It performs correction/refinement at the current `t`, constructs a tangent or Hermite Taylor predictor, validates the step by a Krawczyk test over a Taylor-model tube, then advances to the midpoint of the validated enclosure.
- Step growth is simple: accepted step doubles `h` up to `0.5`; rejected step halves `h`.
- `track_path_adaptive_precision` is an outer retry loop over complete path attempts, not an inner adaptive precision policy.

### Refinement Step

- `refine_moore_box(sys, x, t, r_init, A_init; tau=0.125)` iterates at most 20 times.
- It first checks `krawczyk_test(sys, y, t, s, U; rho=tau)`.
- If validation passes, it tries to grow radius by powers of two.
- If validation fails, it computes `delta = U * evaluate_H(sys, y, t)`.
- If `norm_inf(delta) <= (1/64) * tau * Float64(s)`, it shrinks radius.
- Otherwise it applies a midpoint Newton correction and recomputes `U = inv_acb(evaluate_Jac(sys, y, t), sys.CC)`.

The criteria are certification-sensitive and should not be altered in performance PRs unless there is a separate design note and tests.

### Validation / Krawczyk Step

- Point validation: `krawczyk_test(sys, x, t, r, A; rho=0.7)`.
- Step validation: `validate_step_taylor3(sys, X_tm, t_start, h, r, A; rho=0.7)`.
- Both build a unit complex interval box `B`, evaluate `F`, evaluate `J` on an expanded interval box, form:

```text
K = -(A * F) / r + (I - A * J) * B
```

- Validation accepts when `norm_inf(K) < rho`.

This matches the Moore/Krawczyk flavor that must remain stable. Allocation cleanup should preserve this expression exactly.

### Polynomial / Jacobian Evaluation

- Compiled path uses `Symbolics.build_function` in `compile_homotopy` and `compile_edge_homotopy`.
- `evaluate_H`, `evaluate_Jac`, `evaluate_dt` call `sys.compiled.func_H`, `sys.compiled.func_Jx`, and `sys.compiled.func_dt`.
- For homogeneous/projective systems, `evaluate_H` appends a patch equation and `evaluate_Jac` appends a patch row.
- Legacy path uses `AbstractAlgebra.evaluate`, `jac(system)`, and symbolic `derivative`.

Compiled evaluation is conceptually the right layer to optimize. The current closure outputs are still allocated arrays, and patch augmentation allocates new arrays.

### Interval Arithmetic Backend

- The certified backend is Nemo/Arb: `AcbFieldElem` for complex balls and `ArbFieldElem` for real balls.
- `interval_arithmetic.jl` provides midpoint, norm, conversion, and double display helpers.
- There is no outward-rounded `Float64` interval backend for certified tracking.
- `convert_to_double_*` helpers are for conversion/display/SVD-style approximation and are not a certified double interval fast path.

## Likely Allocation Hotspots

These are code-inspection candidates. They should be confirmed with benchmark and allocation traces before changing source.

| Area | Candidate hotspot | Why it likely allocates |
| --- | --- | --- |
| `track_path` loop | `X_tm = [TaylorModel3(...)]` and `construct_hermite_predictor_tm(...)` | New `Vector{TaylorModel3}` every trial step, including rejected steps. |
| `track_path` accept path | `cache = TMCache(CC)` and `x_next_interval = [...]` | Cache and result vector allocated on every accepted step. |
| `krawczyk_test` | `B = [b_int for _ in 1:n]`, `x_expanded = x .+ ...`, `I_mat`, `term1`, `term2`, `K` | Fresh vectors/matrices and broadcast temporaries per validation. |
| `validate_step_taylor3` | `F_val`, `X_bound`, `B`, `X_expanded`, `I_mat`, `K` | Same pattern inside the innermost step rejection loop. |
| `evaluate_H_augmented` | `[val_sys; patch]` | Allocates a new vector for homogeneous/projective systems. |
| `evaluate_Jac` | `J_aug = Matrix{AcbFieldElem}(undef, ...)` | Allocates and copies Jacobian when homogeneous/projective. |
| `compute_velocity` | `x_mid = get_mid_vec(x)` and `-A * Ht` | Broadcasted midpoint vector and matrix-vector result. |
| `compute_preconditioner` | `get_mid_vec`, `evaluate_Jac`, `inv_acb` | Matrix conversion into Nemo matrix and back to Julia matrix. |
| `TaylorModel3` arithmetic | `c::Vector{AcbFieldElem}`, `C = [CC(0) for _ in 1:4]`, broadcasted `a.c .+ b.c` | Coefficients are heap vectors instead of fixed-size storage. |
| `TaylorModel3` multiplication | repeated `evaluate_taylor(TaylorModel3(...))` | Temporary TaylorModel objects and interval powers. |
| Projective chart logic | `[mag_complex(xi) for xi in x]`, `x ./ scale` | Allocates during chart checks/switches. |
| Legacy path | `hcat`, `transpose`, `Matrix`, `jac`, `derivative`, `evaluate_matrix` | Symbolic objects and matrices rebuilt inside loops. |
| `interval_arithmetic.jl` | `convert_to_box_int` uses `result = []` | `Vector{Any}` and push-based growth. |

The highest-leverage first target is not polynomial math; it is the repeated allocation inside Krawczyk validation and Taylor-model construction in `track_path`.

## Type-Instability Risks

These are risks, not confirmed regressions. Confirm with `@code_warntype`, JET, or SnoopCompile-style inspection after benchmarks exist.

- `evaluate_H_augmented(sys, x, t)` accepts untyped `x` and `t`, and returns either raw compiled output or augmented vectors depending on `sys.homogeneous`.
- `evaluate_Jac(sys, x, t)` returns the raw compiled Jacobian for affine systems but a manually allocated `Matrix{AcbFieldElem}` for homogeneous systems.
- `HCSystem` stores `p_start`, `p_end`, `p_const`, and `patch_vector` as `Tuple{Vararg{AcbFieldElem}}` without length in the type. Splatting these into compiled functions may be hard for inference on some call sites.
- `TaylorModel3.c::Vector{AcbFieldElem}` does not encode coefficient count in the type. Operations on `c` allocate and may hide fixed-size opportunities.
- Many signatures use `Number`, `Vector`, or unparameterized `Matrix`, especially in legacy APIs and predictor code.
- `tracking` mode is a `String`, which forces runtime string comparison in legacy `hermite_tracking`.
- `r` is mixed between `Float64`, `ArbFieldElem`, and `Number`; `validate_step_taylor3(..., Float64(r), ...)` is correctness-sensitive if precision policy changes later.
- `convert_to_box_int` returns a dynamically typed vector because it starts with `result = []`.
- `to_acb_vec` is exported but no definition was found in the current tree; this is an API hygiene issue rather than a direct tracking hotspot.

## Comparison With Algpath Concepts

### Certified Corrector-Predictor Loop

Algpath separates the loop into `Refine -> Extend -> ToChart` components. `GOneStep::one_step` first changes chart when needed, refines the Moore box, then extends with a certified Taylor-model predictor. The current `track_path` has the same conceptual pieces, but they are bundled in one function with allocation-heavy glue.

Recommended direction: do not rewrite the loop wholesale. First extract small scratch/caching structures around existing `refine_moore_box`, `compute_velocity`, and `validate_step_taylor3` while preserving current criteria.

### Tangent Predictor

Algpath computes speed as `-A * fdot` after refinement. Current `compute_velocity` does the same via `evaluate_dt` and `-A * Ht`. The first step in `track_path` is a linear Taylor model `x + v*t`.

Recommended direction: keep the formula. Optimize allocation around midpoint extraction, `evaluate_dt`, and matrix-vector multiply.

### Hermite Predictor

Algpath uses a cubic Hermite predictor when previous point and previous speed are available. Current `construct_hermite_predictor_tm` also builds a cubic predictor from `x`, `x_prev`, `v`, `v_prev`, and `h_prev`.

Recommended direction: keep the coefficients, but make construction allocation-light and cacheable. Any coefficient change should be treated as a mathematical change.

### Taylor Model Use

Algpath has generic order `N` Taylor models with explicit domains and enclosure operations. Current code has `TaylorModel3` only, with four coefficients, a remainder, and a real interval domain `h`.

Recommended direction: before increasing order or changing remainder propagation, make `TaylorModel3` type-stable and benchmarked. Higher-order Taylor models should be a separate design note.

### Double Interval Fast Path

Algpath's Reckless backend is a fast outward-rounded double interval path, selected when precision is at or below 53 bits and conversion is safe. Current code has no equivalent certified double interval backend.

Recommended direction: do not add this early. First add benchmarks and isolate arithmetic interfaces. A Julia fast path needs explicit outward rounding, inclusion conversion tests, and automatic fallback to Acb on failure.

### Arb Fallback

Algpath can try Reckless and fall back to Arb when the fast path cannot certify. Current code always uses Acb/Arb, so Arb is the primary backend rather than a fallback.

Recommended direction: after a fast path exists, fallback should be per step, not per entire path, and should return the same `TrackResult` surface.

### Adaptive / Mixed Precision

Algpath has a precision context that can increase precision during refinement or extension and mixed one-step dispatch. Current `track_path_adaptive_precision` retries the entire path at scheduled precisions and returns the first success.

Recommended direction: keep the current wrapper for API compatibility, but design an isolated `PrecisionPolicy` before inner-loop adaptive precision. Avoid hidden global precision state except a clearly isolated policy object.

## Benchmarking Requirements Before Optimization

Before any performance PR, add a benchmark harness that records:

- Julia version, package commit, system precision, number of paths.
- `time_s`, `alloc_MB`, `iterations`, `accepted_steps`, `rejected_steps`.
- `status`, `success`, `final_radius`, `final_krawczyk_norm`.
- Maximum final residual from `evaluate_H`.
- Surface optimization case from `examples/irreducible_surface_optimization/surface_optimization.jl`.

Suggested benchmark tiers:

1. `small_quadratic`: fast sanity path for CI-like runs.
2. `27_lines` or existing projective benchmark: medium path count and projective behavior.
3. `surface_optimization`: target problem; allow a short deterministic segment for default runs and a full profile for manual runs.

Every optimization PR should include a table like:

```text
case                  mode      time_s before  time_s after  alloc_MB before  alloc_MB after  status
small_quadratic       fixed     ...
surface_optimization  fixed     ...
surface_optimization  adaptive  ...
```

This audit-only document does not provide before/after numbers because it intentionally does not modify executable code.

## Proposed PR Sequence

Each PR below should stay under about 300 changed lines and include benchmark output in the PR description after PR 1 lands.

### PR 1: Baseline Benchmarks Only

- Add `benchmark/run_tracking_baseline.jl` or standardize the existing example benchmark into a committed benchmark entry point.
- Cover `small_quadratic` and one deterministic `surface_optimization` segment.
- Print TSV or JSON with time, allocations, path stats, and validation status.
- No source optimization.

### PR 2: Allocation Trace and Inference Notes

- Add developer documentation for running `@allocated`, `@timed`, and `@code_warntype` on `track_path`, `krawczyk_test`, `validate_step_taylor3`, and `construct_hermite_predictor_tm`.
- Optionally add a non-exported diagnostic helper under `benchmark/`.
- No mathematical or public API change.

### PR 3: Krawczyk Scratch Allocation Cleanup

- Replace obvious vector comprehensions and broadcast temporaries in `krawczyk_test` and `validate_step_taylor3` with local loops or small scratch buffers.
- Preserve `K = -(A * F) / r + (I - A * J) * B` and all thresholds.
- Add regression tests that compare success status and `final_krawczyk_norm` within a tight tolerance on small examples.

### PR 4: `TaylorModel3` Fixed-Size Coefficients

- Change `TaylorModel3` coefficient storage from `Vector{AcbFieldElem}` to fixed-size storage, such as explicit fields or `NTuple{4,AcbFieldElem}`.
- Keep constructors and public behavior compatible.
- Benchmark Taylor-model multiplication and `track_path`.
- Stop and write a design note if remainder propagation becomes unclear.

### PR 5: Compiled Evaluation Buffering

- Add internal `evaluate_H!`, `evaluate_dt!`, and possibly `evaluate_Jac!` wrappers if `Symbolics.build_function` output can be copied into caller-provided buffers safely.
- Avoid patch-vector concatenation allocation in homogeneous/projective systems.
- Keep existing `evaluate_H`, `evaluate_Jac`, and `evaluate_dt` APIs as allocating wrappers.

### PR 6: Preconditioner and Matrix Conversion Cleanup

- Profile `inv_acb` and `compute_preconditioner`.
- Reduce repeated Julia matrix to Nemo matrix to Julia matrix conversions where safe.
- Do not change inverse criteria or fallback behavior.

### PR 7: Legacy Path Cache or Deprecation Note

- Either cache `hcat(H)`, symbolic `jac(H)`, and `dH/dt` inside legacy `track`, or document the compiled `HCSystem` path as the performance path.
- Preserve `track(H, x)` behavior.
- Keep this isolated from compiled-path optimizations.

### PR 8: Inner Precision Policy Design

- Write a design note for `PrecisionPolicy` that can increase precision inside `refine_moore_box` or step validation.
- Keep `track_path_adaptive_precision` as a compatibility wrapper.
- Avoid package-wide mutable global state except an explicitly isolated precision policy.

### PR 9: Double Interval Fast Path Design

- Design first, code later.
- Specify outward rounding, conversion from Acb to fast intervals, certification failure semantics, and Acb fallback.
- Prototype only one validation component before integrating the full one-step tracker.

### PR 10: Surface Optimization Target Benchmark

- Promote `surface_optimization` to the main manual benchmark.
- Add deterministic parameter vertices and fixed starting solution.
- Require benchmark tables for fixed precision and adaptive precision modes before changing algorithms further.

## Highest-Confidence Next Step

Add the benchmark harness first. Then target allocation cleanup in `validate_step_taylor3` and `krawczyk_test`, because those are inside the innermost accept/reject loops and can be improved without changing mathematical criteria.
