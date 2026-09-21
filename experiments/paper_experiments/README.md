# Paper experiments

This directory contains the executable Julia code and tabular data for the four
subsections of the paper ``A priori bounds for certified Krawczyk homotopy tracking''. 

## Setup

From a clone of this repository, check out the `complexity-experiment` branch
and instantiate the package environment:

```sh
git switch complexity-experiment
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Every script below imports `CertifiedHomotopyTracking` from this repository and
uses paths relative to its own `@__DIR__`. It therefore does not depend on the
caller's working directory or on files outside the repository.

## Layout and commands

1. `01_benchmark_examples`
   Constant-predictor statistics for Katsura and dense quadratic systems.

   ```sh
   julia --project=. experiments/paper_experiments/01_benchmark_examples/run.jl
   CHT_FULL=1 julia --project=. experiments/paper_experiments/01_benchmark_examples/run.jl
   ```

2. `02_predictor_order_tradeoff`
   One fixed path of dense cubic systems, predictor orders 0 through 6.

   ```sh
   julia --project=. experiments/paper_experiments/02_predictor_order_tradeoff/run.jl
   CHT_FULL=1 julia --project=. experiments/paper_experiments/02_predictor_order_tradeoff/run.jl
   ```

3. `03_univariate_validation`
   Exact scalar experiment comparing the intrinsic and fixed radii.

   ```sh
   julia --project=. experiments/paper_experiments/03_univariate_validation/run.jl
   CHT_FULL=1 julia --project=. experiments/paper_experiments/03_univariate_validation/run.jl
   ```

4. `04_computational_efficiency`
   Adaptive constant, a priori constant, Hermite, and a priori order-3
   predictors on all benchmark paths.

   ```sh
   julia --project=. experiments/paper_experiments/04_computational_efficiency/run.jl
   CHT_FULL=1 julia --project=. experiments/paper_experiments/04_computational_efficiency/run.jl
   ```

The default invocation is a smoke run. `CHT_FULL=1` selects dimensions 3-6
and all paths. Full runs, especially constant prediction in dimensions 5 and 6,
can take a long time. `CHT_MAX_PATHS=N`, `CHT_ONLY_N=N`, and
`CHT_ONLY_FAMILY=katsura|random_dense` can restrict experiment 4. Set
`CHT_GAMMA_MODE=scalar_start` or `diagonal_start` to select the gamma trick.

## Data policy

The scripts seed the random dense systems and gamma trick. Timing is
machine-dependent, and changing a seed or gamma mode changes iteration counts.
The `beltran_leykin` column in experiment 4 is archival comparison data; this
repository does not contain the external alpha-theory tracker used to produce
that column.
