# H2-R trait-loss-only calibration

## Purpose

H2-A remains the original fixed-absolute-threshold proposition. Its stationary
mutation-primary run was right-censored, not a false result. H2-R is a separate
conditional proposition about baseline-relative genetic change under a declared
monotone deterioration of interaction support.

Before any H2-R warning statistic is evaluated, this calibration selects a
support-deterioration schedule using **post-baseline realised high-trait loss
only**. It is deliberately blind to `H_alpha`, `H_gamma`, relative-warning
times, lead counts, and lead-time magnitudes.

## State path kept fixed

Every attempted record follows the already merged primary chain path:

```text
frozen mutation-H1 primary cell
-> new-seed nested H1 calibration
-> rising high-route replay and full-state hold
-> conservation-preserving equal-isolated projection
-> externally specified monotone barrier schedule
```

H1, H3, the symmetric allele-mutation map, and the full-state projection rule
are not changed. The deterioration schedule modifies only the barrier in the
interaction update after the projected baseline.

## Deterioration schedule

Let \(\theta_0\) be the H1 interior anchor and \(w_\theta\) the canonical
bistable-barrier interval width for that cell. With horizon \(T\) and normalized
increase \(d\), generation \(g=1,\ldots,T\) uses

\[
\theta_g=\theta_0+w_\theta d\frac{g}{T}.
\]

Generation zero remains at the anchor. Thus a loss already present at baseline
is ineligible for calibration and is never counted as a deterioration-induced
trait-loss event.

## Fixed calibration grid

```text
independent calibration master seeds: 20260910–20260914
replicates: 5 per primary cell per seed
horizons: 60, 120
normalized total barrier increases: 0.15, 0.30, 0.45
primary cells: 12
```

Each H1-conditioned source is evaluated at all six schedules with common random
numbers. This makes schedule differences paired within a source; it does not
change the source state or assign a schedule a favourable seed.

## Schedule selection rule

For each cell and schedule, calculate, separately in every calibration seed
block,

\[
P(0<\tau_{T}\le T\mid\text{H1 source prepared, projection supported,
baseline realised high trait present}).
\]

A schedule is eligible only when every seed block has a probability in
\([0.30,0.70]\). Among eligible schedules, choose the pooled probability closest
to 0.50; ties go to the shorter horizon and then the smaller normalized barrier
increase.

No schedule is selected when no candidate qualifies. The artifact retains all
candidates, source failures, projection failures, and baseline-ineligible
records.

## What calibration does not infer

This run does not test H2-R. It only produces the predeclared deterioration
schedule for a later independent H2-R validation. In particular, a selected
schedule says nothing about whether relative H-alpha or H-gamma decline leads
trait loss.
