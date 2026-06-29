# Finite H1 branch-separation audit

The canonical H1 certificate proves bistability and a high-trait viability
switch only for the one-state interaction map. The finite-bin multipatch
simulator contains density variation, realised-trait feedback, allele feedback,
finite recruitment, and optional migration. It therefore needs its own
simulation audit.

## Question

For a parameter cell where canonical H1 has two stable interaction branches and
a high-trait margin switch, does the finite coupled simulator retain a lasting
difference when started from the canonical low versus high branch?

## Matched branch pairs

For each landscape, parameter cell, replicate index, and seed, the audit makes
two runs that differ only in initial interaction:

- `low_start` initializes every patch at the canonical low stable interaction;
- `high_start` initializes every patch at the canonical high stable interaction.

Population, trait distribution, allele frequency, landscape, parameter cell,
replicate index, and random seed are held fixed. The same seed is intentional:
it makes stochastic variation comparable rather than treating unrelated random
runs as evidence of branch separation.

## Predeclared finite predicates

For every branch pair, the audit records the mean interaction over the final
`terminal_window` snapshots. Branch separation is supported when

```text
mean_q(high_start, terminal window) - mean_q(low_start, terminal window)
    > interaction_separation_threshold.
```

The stronger `finite_h1_mechanism_supported` predicate additionally requires a
terminal potential high-trait switch:

```text
high_start: potential high trait viable
low_start:  potential high trait not viable.
```

Terminal realised high-trait mass and local effective-size differences are
reported separately. They are finite-model consequences, not additional clauses
of the canonical H1 theorem.

## Noncanonical cells

Cells that lack a canonical branch-dependent high-trait certificate are retained
in the output with finite branch pairs unavailable (`None`). They are not
reported as failed finite branch separation: the audit was not applicable to
them. This preserves the distinction between a negative result in the testable
region and a cell outside the declared H1 mechanism.

## Scope and interpretation

Per-run H1 theorem-boundary audits accompany low and high starts. These show
whether either finite trajectory happens to occupy the exact canonical limit.
The finite branch-separation predicate remains Type S evidence even when a
canonical context exists. It is not a proof of bistability, and it is not yet a
full hysteresis test.

A hysteresis test is the next stricter step: it must continue a state through
rising and falling barrier paths and show different jump boundaries under the
finite closure.
