# %%
# Load package and tools

using CertifiedHomotopyTracking
using Printf

include("track_examples.jl")

# %%
# Build benchmark cases

cases = make_cases()

# %%
# Warm up one case

warmup_result = track_case(cases[1])
@show succeeded(warmup_result)
@show warmup_result.final_t
@show warmup_result.iterations
@show warmup_result.accepted_steps
@show warmup_result.rejected_steps

# %%
# Run a quick timed measurement on the simple case

simple_row = timed_track_case(cases[1])
print_csv_header()
print_csv_row(simple_row)

# %%
# Inspect time and allocations with base Julia macros

case = cases[2]
track_case(case)
@time track_case(case)
@allocated track_case(case)

# %%
# Run all quick cases and print a compact table

rows = [timed_track_case(case) for case in cases]

@printf("%-20s %6s %8s %10s %14s %6s %6s %6s %8s %12s\n",
    "case", "bits", "success", "time_s", "alloc_bytes", "iter", "acc", "rej", "final_t", "radius")
for row in rows
    @printf("%-20s %6d %8s %10.4g %14d %6d %6d %6d %8.4g %12.4g\n",
        row.case,
        row.precision_bits,
        row.success,
        row.elapsed_seconds,
        row.allocated_bytes,
        row.iterations,
        row.accepted_steps,
        row.rejected_steps,
        row.final_t,
        row.final_radius,
    )
end

# %%
# Print CSV-like output for copying into before/after files

print_csv_header()
for row in rows
    print_csv_row(row)
end

# %%
# Optional: run only the projective case again

projective_case = last(cases)
projective_row = timed_track_case(projective_case)
print_csv_header()
print_csv_row(projective_row)

# %%
# Suggested before/after workflow
#
# git checkout main
# julia --project=. benchmark/track_examples.jl > before.csv
#
# git checkout optimize-taylor-cache
# julia --project=. benchmark/track_examples.jl > after.csv
#
# Compare before.csv and after.csv manually, in a spreadsheet, or with a small script.
