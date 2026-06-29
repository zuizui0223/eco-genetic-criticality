# Nested-grid finite H1 boundary refinement audit

The endpoint-expanded sweep audit established that a finite collapse/recovery
loop can be observed after the continuation range is widened. Its first reported
collapse and recovery locations are still grid points, not exact finite critical
barriers. If every replicate reports the same grid point, that can mean the
finite response is highly robust, but it can also mean the chosen barrier grid
is the limiting resolution.

This audit resolves that ambiguity by repeating the same finite continuation on
nested barrier grids at one common endpoint range.

## Confirmation range

The first endpoint study found that an endpoint padding of one-half the canonical
interval width closes a threshold-defined loop for all standard parameter pairs.
The confirmation audit therefore uses a **single common** padding fraction:

```text
f = 0.5
```

for every canonical-H1-applicable `(A_ref, kappa)` pair. It does not replace the
original endpoint results or choose a pair-specific range based on where that
pair first closed.

For an interval of width `W`, every grid spans:

```text
[theta_lower - 0.5 W, theta_upper + 0.5 W].
```

## Nested grids

The first resolution ladder is:

```text
25, 49, 97 barrier points.
```

These correspond to 24, 48, and 96 subdivisions. Every coarse grid point is
therefore also present on the next finer grid. The barrier step halves at each
refinement without changing the sweep endpoints.

A new master seed should be used for the confirmation run so that agreement is
not merely a replay of the endpoint-discovery stochastic stream. Within that
new run, seeds and replicate indices are paired across grid resolutions.

## Brackets, not point estimates

For each route transition, the audit records adjacent barriers.

For the rising route, if terminal interaction first falls to or below the low
threshold between two successive barriers, the collapse is reported as:

```text
[previous_barrier, first_crossing_barrier].
```

For the falling route, if terminal interaction first rises to or above the high
threshold while barriers decrease, recovery is reported as:

```text
[first_crossing_barrier, previous_barrier].
```

Thus every finite boundary is an interval with width equal to one declared grid
step. The finite loop gap is likewise an interval:

```text
collapse_lower - recovery_upper
  <= finite_gap
  <= collapse_upper - recovery_lower.
```

A loop is bracket-supported only when the **lower bound** of this gap is strictly
positive. This is intentionally more conservative than checking only the two
first-crossing grid points.

## Resolution-stability predicate

A replicate is labelled `resolution_stable_loop_supported=True` only when both
two finest grids show a bracket-supported loop, all of the following hold, and
the finite closure already satisfies the standard trait-switch mechanism:

1. the collapse-bracket midpoint moves by no more than one penultimate-grid
   step;
2. the recovery-bracket midpoint moves by no more than one penultimate-grid
   step; and
3. the widest finest-grid boundary bracket is at most `0.03` of canonical
   interval width.

The corresponding stronger predicate
`resolution_stable_h1_loop_mechanism_supported=True` keeps the same geometric
requirements and explicitly requires finite potential high-trait switch support
on both finest grids.

These are predeclared operational Type S checks. They demonstrate grid-resolution
robustness for a particular finite closure; they do not prove an exact finite
bifurcation or a universal hysteresis theorem.

## Running

```bash
python scripts/run_finite_h1_boundary_resolution_audit.py \
  --profile standard \
  --master-seed 20260630 \
  --endpoint-padding-fraction 0.5 \
  --stage-generations 30 \
  --barrier-points 25 \
  --barrier-points 49 \
  --barrier-points 97 \
  --maximum-normalized-bracket-width 0.03 \
  --output-dir artifacts/h1_boundary_resolution \
  --prefix standard_resolution_v1
```

The JSON preserves every route and all three resolution levels. The CSV is a
pair-level index. The manifest fixes the common range, nested grid design,
thresholds, resolution-stability criterion, seed, source revision, and
no-selection policy.
