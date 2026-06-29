# Finite H1 hysteresis stage-duration robustness audit

The finite H1 continuation audit reports route memory after a declared number
of generations at each barrier stage. A large rising-versus-falling difference
can be scientifically interesting, but it is not automatically quasi-static
hysteresis: it may be a transient left over because a stage was too short to
relax.

This audit asks whether the route-memory result survives when the same
one-large continuation protocol is repeated at progressively longer stage
durations.

## Why this is separate from the targeted barrier grid

A rising/falling continuation spans the full canonical bistable interval. Its
trajectory is defined by `(area_reference, interaction_feedback)`, barrier grid,
seeds, and stage duration—not by an arbitrary fixed raw barrier selected for a
phase-diagram row. Therefore this audit creates exactly one midpoint-labelled
one-large family per `(A_ref, kappa)` pair rather than repeating equivalent
continuations for every targeted fixed-barrier position.

Pairs without a strict one-large canonical bistable interval are retained as
**unavailable**. They are not converted into finite failures.

## Paired duration ladder

Use an ordered ladder, for example:

```text
5, 10, 30, 80 generations per barrier stage.
```

Each duration repeats the same rising and falling barrier sequence. Within one
`(A_ref, kappa, replicate)` family, the pair-specific master seed and replicate
index are fixed across all durations. The route's stage seeds are therefore
paired by route and stage index. The purpose is not to claim independent random
replicates across durations; it is to isolate the effect of allowing more time
per continuation stage.

## Predeclared quantities

For every duration, retain the existing finite H1 continuation outputs:

- route-memory support;
- route-memory plus potential high-trait switch support;
- maximum interaction gap at shared internal barriers;
- potential-switch barrier count; and
- thresholded jump proxy when it exists.

The new summary also records:

```text
longest_pair_gap_change
  = abs(maximum_gap_at_longest_duration
        - maximum_gap_at_penultimate_duration)
```

and compares it with `gap_stability_tolerance`.

A replicate has `convergence_robust_hysteresis_supported=True` only when:

1. finite route-memory support is present at both longest durations; and
2. `longest_pair_gap_change <= gap_stability_tolerance`.

The stronger `convergence_robust_h1_mechanism_supported=True` additionally
requires the potential high-trait switch predicate at both longest durations.

These are intentionally operational Type S predicates. They do **not** prove
that every stage reached an equilibrium or that arbitrary ecological systems
have hysteresis.

## How to read the outcome

- **Support persists and gap stabilizes:** stronger finite evidence that the
  continuation result is not simply disappearing as stages lengthen.
- **Support disappears at longer stages:** earlier route memory was duration
  dependent; report it as finite transient memory, not quasi-static hysteresis.
- **Support persists but gap changes materially:** the result is sensitive to
  relaxation duration and needs a longer or finer ladder before a stable finite
  description is claimed.
- **Jump proxy remains unavailable:** route-memory may still be present, but the
  chosen low/high interaction thresholds did not locate clean finite jump
  boundaries. Do not infer jump locations from a missing proxy.

## Running

```bash
python scripts/run_finite_h1_hysteresis_duration_audit.py \
  --profile standard \
  --stage-generations 5 \
  --stage-generations 10 \
  --stage-generations 30 \
  --stage-generations 80 \
  --output-dir artifacts/h1_hysteresis_duration \
  --prefix standard_duration_v1
```

The JSON retains all duration-specific finite continuation records. The CSV is
a compact pair-level index. The manifest freezes the duration ladder, thresholds,
seed, source revision, and no-selection rule.
