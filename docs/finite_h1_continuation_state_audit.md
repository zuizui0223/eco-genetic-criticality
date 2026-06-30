# Finite H1 full-state continuation hold audit

## Why this audit exists

The H1 hysteresis and boundary-resolution campaigns use a continuation protocol.
At every barrier stage, the complete terminal finite state becomes the starting
state for the next stage:

```text
population
interaction q
high-allele frequency p
realised trait-bin abundance
```

A later branch-locked H1-H2-H3 campaign restarted from a canonical q value only.
That restart lost the other three state variables and all high starts collapsed.
It therefore cannot be interpreted as a failure of the finite continuation loop
itself, nor as a valid H2/H3 test conditional on that loop.

This audit tests the missing bridge before H2 or H3 is attempted again.

## Design

For every declared master seed, parameter pair, and replicate:

1. Reproduce the validated H1 finite boundary calibration using common endpoint
   padding `0.5`, 30 generations per stage, and nested grids `25,49,97`.
2. Take the finest-grid conservative interior:

   ```text
   (falling recovery bracket upper, rising collapse bracket lower)
   ```

3. Choose the finest-grid barrier nearest the interior midpoint, but only if it
   lies strictly inside that open interval.
4. Replay both continuation routes to that exact barrier:

   ```text
   rising route  -> high route state
   falling route -> low route state
   ```

5. Carry the full terminal state of each route into a fresh 30-generation hold
   at the same barrier, using fresh deterministic hold seeds.

The primary predicate requires both source and held states to have high-minus-low
interaction larger than the predeclared separation threshold, with potential
high-trait viability present only on the high route.

## Interpretation

Passing demonstrates that the finite H1 route memory remains after a fresh hold
when the whole continuation state is transferred. It does **not** prove an
exact finite bifurcation, and it does not yet establish any H2 warning order or
H3 fragmentation effect.

Failing would be informative too: it would show that the currently observed
finite loop is continuation-protocol memory that does not survive an additional
fixed-barrier interval, even with full-state transfer.

## Standard run

```bash
python scripts/run_finite_h1_continuation_state_audit.py \
  --profile standard \
  --master-seed 20260630 \
  --master-seed 20260631 \
  --master-seed 20260632 \
  --master-seed 20260633 \
  --master-seed 20260634 \
  --endpoint-padding-fraction 0.5 \
  --stage-generations 30 \
  --hold-generations 30 \
  --barrier-points 25 \
  --barrier-points 49 \
  --barrier-points 97 \
  --interaction-separation-threshold 0.05 \
  --maximum-normalized-bracket-width 0.03 \
  --output-dir artifacts/h1_continuation_state \
  --prefix standard_full_state_hold_v1
```

The CSV is a parameter-pair ledger. JSON preserves every calibration record,
selected interior barrier, source state, hold terminal, and unavailable case.
The manifest records the carried variables and no-selection policy.