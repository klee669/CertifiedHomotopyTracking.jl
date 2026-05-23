# Benchmarks

이 디렉터리는 패키지 동작을 바꾸지 않고 `track_path` 성능을 측정하는 benchmark suite를 담고 있습니다.

## Setup

benchmark 전용 environment를 instantiate합니다.

```sh
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
```

benchmark script는 package root를 `LOAD_PATH`에 추가하므로, 현재 checkout과 package root의 일반 dependency를 그대로 사용합니다.

## Run

기본 case 전체를 실행합니다.

```sh
julia --project=benchmark benchmark/bench_track_path.jl
```

빠른 smoke run은 다음처럼 실행합니다.

```sh
CHT_BENCH_SAMPLES=1 CHT_BENCH_SECONDS=1 julia --project=benchmark benchmark/bench_track_path.jl
```

특정 case만 실행하려면 `CHT_BENCH_CASES`를 사용합니다.

```sh
CHT_BENCH_CASES=readme_straight_line,small_dense julia --project=benchmark benchmark/bench_track_path.jl
```

사용 가능한 case:

- `readme_straight_line`: README straight-line homotopy example.
- `small_dense`: small dense quadratic polynomial system.
- `monodromy_edge`: README-style parameter system에서 나온 deterministic monodromy edge 하나.
- `surface_optimization`: `examples/irreducible_surface_optimization/surface_optimization.jl`에서 나온 deterministic edge 하나.

## Output

script는 다음 항목을 tab-separated row로 출력합니다.

- runtime: `median_time_ns`, `min_time_ns`
- allocations: `memory_bytes`, `allocations`
- tracking steps: `iterations`, `accepted_steps`, `rejected_steps`
- validation attempts: `accepted_steps + rejected_steps`
- refinement attempts: 현재 public result가 refinement iteration count를 노출하지 않으므로 `not_accessible`로 기록합니다.
