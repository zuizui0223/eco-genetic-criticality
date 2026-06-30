# Independent validation of the mutation-H1 polymorphism window

## Status of the preceding screen

The first mutation-window run used five master seeds
`20260630`–`20260634` and five replicates per parameter pair. It was an
exploratory Type S screen. It identified three mutation rates to validate:

```text
mu = 0.10, 0.15, 0.20
```

This document defines the next run before looking at its results.

## Fixed independent-validation protocol

The workflow `H1 Mutation Independent Validation` has no mutable inputs. It
runs exactly:

```text
profile: standard
master seeds: 20260710, 20260711, 20260712, 20260713, 20260714
mutation rates: 0.10, 0.15, 0.20
replicates: 20 per parameter pair per master seed
endpoint padding fraction: 0.5
stage generations: 30
hold generations: 30
nested barriers: 25, 49, 97
interaction separation threshold: 0.05
maximum normalised bracket width: 0.03
polymorphism epsilon: 1e-12
```

The validation seeds do not overlap with the screen seeds. The run contains

\[
3\ \text{rates} \times 5\ \text{master seeds} \times 9\ \text{parameter pairs}
\times 20\ \text{replicates}=2700
\]

seed-replicate records.

## Unit of interpretation

Do not choose a mutation rate from a pooled mean alone. The reporting unit is

\[
(\mu,\ A_{\rm ref},\ \kappa,\ \text{master-seed block}).
\]

For every mutation-rate and parameter-pair cell, report:

```text
h1_full_state_hold_supported_probability
h2_genetic_baseline_eligible_probability
screen_supported_probability
high_branch_allele_frequency.mean
high_branch_H_alpha.mean
by_master_seed.<seed>.screen_supported_probability
```

A record contributes `screen_supported=True` only when the same replicate has
both a finite H1 full-state high/low hold and an eligible polymorphic high
baseline. A high aggregate rate with one or more failing seed blocks must remain
reported as heterogeneous validation evidence.

## Boundary of inference

This run validates the numerical domain in which the symmetric-mutation closure
simultaneously supports finite H1 memory and genetic variation. It does not test
H2 lead times or H3 fragmentation outcomes.

Only after inspecting all 27 rate-by-pair cells and their five independent seed
blocks can a later branch-conditioned H2/H3 campaign be configured. Any such
campaign must retain this mutation closure, the full-state transfer rule, and
all nonqualifying cells in its artifacts.

## Artifact identity

The generated manifest contains:

```text
campaign: finite_h1_symmetric_mutation_window_v1
campaign_role: independent_validation_v1
```

The uploaded artifact also includes `validation_protocol.txt`, which records the
fixed values above.
