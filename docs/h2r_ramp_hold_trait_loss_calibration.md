# H2-R ramp-and-hold trait-loss calibration v2

## Why a second calibration family

The first H2-R calibration (`h2r_trait_loss_only_calibration_v1`, Actions run
`28493522149`) used a linear barrier ramp over the entire 60- or 120-generation
observation window. No frozen primary cell had a schedule meeting the
all-seed-block trait-loss availability rule. That outcome remains part of the
record; it is not overwritten by this v2 campaign.

V2 does not modify H1, H3, H2-A, the symmetric mutation closure, the H1
full-state reconstruction, the equal-isolated conservation projection, the
trait-loss endpoint, or the schedule-selection criterion. It changes only the
time profile of the same normalized external barrier increase.

## V2 schedule family

For H1 interior anchor \(\theta_0\), canonical interval width \(w_\theta\),
normalized total increase \(d\), ramp duration \(R=30\), and hold duration
\(L\), generation \(g\) uses

\[
\theta_g=
\begin{cases}
\theta_0+w_\theta d\,g/R, & 1\le g\le R,\\
\theta_0+w_\theta d, & R<g\le R+L.
\end{cases}
\]

The projected baseline is generation zero at \(\theta_0\). The candidate family
is fixed before this run:

```text
ramp duration: 30 generations
hold durations: 90, 210 generations
normalized total increases: 0.15, 0.30, 0.45
candidate total horizons: 120, 240 generations
```

Thus v2 asks whether realised trait loss requires time at the already declared
final deterioration level. It is not a stronger barrier intervention than v1.

## State path and fixed run

```text
frozen mutation-H1 primary cell
-> new-seed nested H1 calibration
-> high-route full-state replay and hold
-> conservation-preserving equal-isolated projection
-> 30-generation barrier ramp
-> 90- or 210-generation barrier hold
```

The fixed manual workflow uses independent master seeds `20261010`–`20261014`
and five replicates per cell per seed. Each reconstructed source is run across
all six schedule candidates using common random numbers.

## Unchanged selection rule

For each cell and candidate, the only measured quantity is

\[
P(0<\tau_T\le R+L \mid \text{source prepared, projection supported,
baseline realised high trait present}).
\]

A candidate is selectable only when **every** seed block has probability within
\([0.30,0.70]\). The chosen schedule is the selectable candidate whose pooled
trait-loss probability is nearest 0.50; ties choose the shorter total horizon,
then the smaller normalized increase. A cell without such a candidate remains
`no_schedule_selected`.

## What v2 does not do

V2 does not evaluate H2-R. It neither computes nor writes `H_alpha`, `H_gamma`,
relative-warning times, warning-lead indicators, or lead-time values. Its sole
purpose is to choose an observable trait-loss schedule for a later independent
H2-R validation.
