# H2 reframing after the mutation-conditioned primary chain

## What remains fixed

Nothing in this document changes H1, H3, the symmetric-mutation closure, the
full-state H1 transfer rule, or the existing H2 absolute-threshold definition.

The existing H2 statement is retained as **H2-A**:

> Under declared conditions, fixed absolute \(H_\alpha\) or \(H_\gamma\)
> warning thresholds can precede realised high-trait loss.

H2-A is not labelled false after the stationary mutation-conditioned primary
chain. It is **unresolved in that campaign**, because most trajectories were
right-censored: realised trait loss and/or a fixed-threshold warning was absent
within the 30-generation horizon. An event order cannot be inferred from a
missing event.

## Why an additional proposition is needed

The primary chain established a different finite result: after a valid H1 high
full-state source was projected at fixed total area, equal isolation consistently
reduced interaction, local effective size, and realised high-trait mass relative
to one large patch. That is H3 evidence. Yet symmetric mutation maintained enough
polymorphism that both absolute genetic warning crossings and realised trait-loss
crossings were sparse in the same horizon.

Thus the stationary run cannot distinguish these possibilities:

1. genetic change does precede loss, but neither event has been observed yet;
2. mutation maintains diversity above fixed absolute thresholds while support
   deteriorates; or
3. an absolute heterozygosity threshold is not the appropriate warning statistic
   for this closure.

The correct response is not to relabel censoring as support or failure. It is to
state a new, narrower, testable dynamic proposition.

## New proposition: H2-R

**H2-R — conditional relative genetic warning.**

Let \(H_x(t)\) denote \(H_\alpha(t)\) or \(H_\gamma(t)\), for
\(x\in\{\alpha,\gamma\}\), along a trajectory starting from a polymorphic,
H1-conditioned full state. Under a **predeclared monotone decline in interaction
support**, define

\[
\tau_{\Delta H_x(r)}
=
\inf\left\{t>0: H_x(t)\le (1-r)H_x(0)\right\},
\qquad r\in\{0.05,0.10,0.20\},
\]

and let \(\tau_T\) be first realised high-trait loss. H2-R asks whether

\[
\tau_{\Delta H_x(r)}<\tau_T.
\]

The inequality is evaluated only when the baseline is eligible and both
first-passage events are observed. Otherwise the record is censored. The
baseline time \(t=0\) is never called an early warning.

This is a **dynamic hypothesis (H)**. Any finite result remains Type S for the
specified mutation and deterioration closure.

## The deterioration schedule is not selected using warning success

For an H2-R campaign, interaction support is reduced by a monotone barrier
schedule expressed relative to that cell's canonical bistable-interval width.
Before examining H-alpha or H-gamma warning outcomes, a short calibration chooses
from this predeclared family:

```text
horizons: 60, 120 generations
total normalized barrier increases: 0.15, 0.30, 0.45
trait-loss target: 0.30–0.70 in every calibration seed block
```

The selected schedule is the eligible candidate whose pooled trait-loss frequency
is closest to 0.50; ties go to the shorter horizon and then the smaller barrier
increase. If no candidate satisfies every seed block, the cell remains
`no_schedule_selected`. The rule does **not** inspect warning lead counts,
lead-time size, or H-alpha/H-gamma values.

## Consequences for interpretation

H2-A and H2-R answer different questions.

| label | question | evidence status now |
|---|---|---|
| H2-A | Does a fixed absolute diversity threshold precede loss? | unresolved because of stationary-run censoring |
| H2-R | Under declared deteriorating support, does baseline-relative diversity erosion precede loss? | untested, now formally specified |

A future H2-R result may support the relative warning, show no lead, show late
warning, or remain censored. None of these outcomes modifies the already fixed
H1 or H3 statements.
