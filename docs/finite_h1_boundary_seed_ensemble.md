# Multi-master-seed finite H1 boundary robustness audit

The nested-grid boundary-resolution audit verifies that a finite
collapse--recovery bracket survives 25 → 49 → 97 barrier points under one
independent master seed. That is a meaningful confirmation, but it does not
show whether the bracket result is robust across independent stochastic streams.

This audit repeats the same predeclared boundary-resolution design across a
fixed set of master seeds and retains every seed-specific raw continuation.

## Predeclared ensemble

The first ensemble uses:

```text
20260630, 20260631, 20260632, 20260633, 20260634
```

The ensemble includes the first independent confirmation seed (`20260630`) and
four additional seeds. It does not select seeds after observing whether they
support a loop.

Within a master-seed run, replicate index and route/grid seeds remain paired as
in the boundary-resolution audit. Across master-seed runs, the base seed changes
all derived stochastic streams.

## Fixed design across seeds

Every seed uses exactly the same finite design:

```text
endpoint padding fraction: 0.5
stage generations:          30
nested barrier grids:       25, 49, 97
maximum bracket width:      0.03 × canonical interval width
```

Every canonical-H1-applicable `(A_ref, kappa)` pair is evaluated under every
seed. There is no pair-specific endpoint choice and no outcome-based seed,
replicate, or parameter filtering.

## Why seed-level and pooled results differ

The output gives both:

- **pooled replicate support:** the fraction of all seed × replicate records
  satisfying the resolution-stable predicate; and
- **master-seed support distribution:** minimum, median, mean, and maximum of
  each seed-run's support probability.

A high pooled fraction can hide a whole seed run that performs poorly. Therefore
the audit also records:

```text
all_master_seed_runs_fully_resolution_stable
all_master_seed_runs_fully_resolution_stable_h1_mechanism
```

These are true only when *every* declared seed has 100% support among its own
finite replicates. They are intentionally stronger than a high pooled average.

## Interpretation

- **All seeds fully support stable brackets:** strong Type S evidence that the
  finite loop and high-trait mechanism are robust to the declared seed ensemble
  and grid refinement.
- **Most seeds support but one or more do not:** report the seed-level range;
  do not replace it with a pooled statement alone.
- **Pooled support is low or seed-level outcomes vary strongly:** the finite
  boundary is stochastic under this closure and needs a probabilistic, rather
  than deterministic, interpretation.

Even complete support does not identify an exact finite bifurcation. The output
is still a numerical bracket under specified thresholds and the declared finite
closure.

## Running

```bash
python scripts/run_finite_h1_boundary_seed_ensemble.py \
  --profile standard \
  --ensemble-master-seed 20260630 \
  --ensemble-master-seed 20260631 \
  --ensemble-master-seed 20260632 \
  --ensemble-master-seed 20260633 \
  --ensemble-master-seed 20260634 \
  --endpoint-padding-fraction 0.5 \
  --stage-generations 30 \
  --barrier-points 25 \
  --barrier-points 49 \
  --barrier-points 97 \
  --maximum-normalized-bracket-width 0.03 \
  --output-dir artifacts/h1_boundary_seed_ensemble \
  --prefix standard_seed_ensemble_v1
```

The JSON contains the complete raw continuation output for each master seed.
The CSV is an aligned pair-level ledger with pooled and seed-level summaries.
The manifest records the seed ensemble, common design, source revision, and
no-selection policy.
