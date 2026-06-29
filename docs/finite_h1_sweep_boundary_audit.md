# Finite H1 endpoint-expanded sweep boundary audit

The stage-duration audit established whether rising/falling route memory remains
when each continuation stage is allowed more time. It does not guarantee that
the finite continuation reaches a barrier far enough from the canonical interval
to observe both route transitions.

This audit expands the continuation endpoints and asks whether the finite
high-start route collapses on the rising sweep and the finite low-start route
recovers on the falling sweep.

## Endpoint design

For a one-large canonical H1 interval

```text
(theta_lower, theta_upper)
```

with width

```text
W = theta_upper - theta_lower,
```

the sweep at endpoint-padding fraction `f` uses

```text
[theta_lower - f W, theta_upper + f W].
```

The first declared ladder is:

```text
f = 0.1, 0.5, 1.0, 2.0.
```

Every `(A_ref, kappa)` pair is run at every declared `f`. No observed finite
outcome determines which endpoint values are retained.

## Finite boundary observations

At each endpoint expansion, the existing continuation audit records:

- `rising_collapse_barrier`: the first rising-route grid barrier whose terminal
  mean interaction is at or below `low_state_threshold`;
- `falling_recovery_barrier`: the first falling-route grid barrier whose terminal
  mean interaction is at or above `high_state_threshold`; and
- `jump_boundary_gap = rising_collapse_barrier - falling_recovery_barrier` when
  both are observed.

A `finite_loop_closed=True` observation requires both boundaries and a strictly
positive gap. The stronger `finite_h1_loop_mechanism_supported=True` also
requires the finite potential high-trait switch predicate.

These are **grid-limited** estimates. The output stores

```text
barrier_step = total_sweep_span / (barrier_points - 1)
```

for every observation. A reported collapse or recovery barrier is an estimate
to within the declared grid resolution, not an exact finite critical value.

## First boundary bracket

For each replicate, the audit saves the smallest padding fraction where a
closed finite loop is observed, together with its two boundary estimates,
finite gap, and grid step. This answers two distinct questions:

1. whether the prior continuation range was too short to observe transitions;
and
2. how far outside the canonical interval the finite closure must be swept
before both transitions appear under the declared thresholds.

It does not use the smallest successful padding as a cell-selection criterion;
all earlier unsuccessful endpoint sweeps remain in the raw output.

## Reading results

- **Loop closes at small padding:** the prior missing jump proxy was primarily
  endpoint-range limited.
- **Loop closes only at large padding:** the finite sweep requires much wider
  forcing than the canonical interval. Report that finite/canonical difference.
- **Only one route transition appears:** finite route memory exists, but the
  current endpoint ladder does not close a threshold-defined loop. Expand the
  ladder or revise thresholds only as a new, declared sensitivity design.
- **No transitions appear:** do not infer a finite jump boundary. The result is
  route memory without a detected threshold-defined loop in the sampled range.

## Running

```bash
python scripts/run_finite_h1_sweep_boundary_audit.py \
  --profile standard \
  --stage-generations 30 \
  --barrier-points 25 \
  --endpoint-padding-fraction 0.1 \
  --endpoint-padding-fraction 0.5 \
  --endpoint-padding-fraction 1.0 \
  --endpoint-padding-fraction 2.0 \
  --output-dir artifacts/h1_sweep_boundaries \
  --prefix standard_sweep_v1
```

### GitHub Actions run

The **H1 Finite Sweep Boundaries** workflow is manual (`workflow_dispatch`).
Use it for the first boundary-localization study so the full design, code
revision, seed, threshold defaults, and raw output remain together as one
artifact-backed record.

1. Open **Actions** and select **H1 Finite Sweep Boundaries**.
2. Use `standard` for the first scientific run; `quick` only verifies runner
   installation and artifact shape.
3. Keep endpoint padding as `0.1,0.5,1.0,2.0`, stage duration as `30`, and
   grid resolution as `25` for the predeclared initial endpoint design.
4. Use master seed `20260629`; leave replicates blank to retain the standard
   profile's declared count and direct comparability with the duration ladder.
5. Use `full` only after explicitly enabling the full-profile approval input.

The workflow uploads CSV, JSON, and manifest for 90 days, and also uploads
partial output on failure. Every padding fraction remains in the artifact even
if no finite threshold-defined loop is detected.

The CSV is a compact pair-level index. JSON preserves every endpoint sweep and
finite trajectory, including unavailable cells and sweeps with no detected
boundary. The manifest freezes thresholds, stage duration, range ladder, grid
resolution, seed, and source revision.
