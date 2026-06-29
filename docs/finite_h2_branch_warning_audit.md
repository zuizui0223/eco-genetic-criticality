# Branch-conditioned finite H2 warning audit

The canonical H2 certificate concerns a deterministic expectation recursion with
constant effective population size and monotone realised-trait decline. The
finite coupled simulator instead has stochastic allele sampling, interaction-
dependent demography, finite trait recruitment, and potentially different H1
branches. This audit asks a narrower, directly testable question:

> When a finite replicate has already shown the finite H1 branch mechanism,
> does a genetic indicator precede realised high-trait loss on its low-start or
> high-start trajectory?

## H1 precondition

The audit reuses the finite H1 branch-separation experiment. A replicate enters
H2 denominators only when its same-seed low-start/high-start pair satisfies the
strong finite H1 predicate:

1. terminal interaction remains separated above the declared threshold; and
2. the high-start path retains potential high-trait viability while the low-start
   path does not.

Canonical-H1-inapplicable cells and finite H1 failures remain in the artifact,
but are not converted into H2 warning failures. This keeps three outcomes
separate:

- no applicable H1 branch mechanism;
- an applicable H1 mechanism that did not persist in the finite replicate; and
- a branch-conditioned H2 comparison with observed first-passage events.

## Branch-specific H2 comparisons

For each H1-preconditioned replicate, low-start and high-start paths are
reported separately. Each path compares the first realised high-trait loss time
with:

- `tau_H_alpha`;
- `tau_H_gamma`; and
- `tau_allele_loss`.

A warning leads only if

```text
tau_warning < tau_trait_realised.
```

The artifact records the signed quantity

```text
lead_time_trait_minus_warning = tau_trait_realised - tau_warning.
```

Positive values mean the warning came first, zero means a tie, and negative
values mean the warning occurred after realised trait loss.

## Censoring and denominators

A warning/trait pair is valid only when both events occur within the declared
simulation horizon. If either event is absent, the comparison is censored.
For every branch and warning type, report together:

1. H1-preconditioned replicate count;
2. valid-pair count;
3. censored-pair count;
4. lead count;
5. lead probability conditional on valid pairs; and
6. lead probability across all H1-preconditioned replicates.

The conditional probability answers whether warning leads **when both events
were observed**. The denominator-wide probability answers how often the full
branch-conditioned warning pattern occurred. Neither should be substituted for
the other.

## Interpretation

This is Type S evidence for the finite closure. It does not upgrade the
canonical H2 expectation theorem into a universal stochastic ordering. It can
show branch-specific genetic lead, no lead, ties, late warnings, or insufficient
observed event pairs. These possibilities are all reportable outcomes rather
than discarded cases.
