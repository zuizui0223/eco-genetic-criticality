# Mutation-H1 primary-analysis domain

## Why this domain is frozen first

The first symmetric-mutation run was exploratory.  A second run used independent
master seeds `20260710`–`20260714`, 20 replicates per parameter pair and seed,
and all three screened mutation rates.  Before H2/H3 trajectories are run, this
document fixes how that validation is translated into a primary analysis domain.

No cell is removed from the validation ledger.  Primary eligibility is a label
for the next dynamic campaign; excluded cells remain reported as a sensitivity
and feasibility boundary.

## Evidence source

```text
Actions run: 28436777080
campaign role: independent_validation_v1
mutation rates: 0.10, 0.15, 0.20
A_ref values: 0.8, 1.0, 1.2
kappa values: 3.0, 4.5, 6.0
independent seeds: 20260710–20260714
replicates: 20 per pair per seed
```

The same-replicate predicate is

```text
finite H1 full-state high/low hold
AND polymorphic high branch
AND H-alpha/H-gamma above their baseline-warning thresholds.
```

## Frozen selection rule

A validation cell is primary-analysis eligible exactly when

\[
\min_{s\in\{20260710,\ldots,20260714\}}
P_s(\text{joint support}) \ge 0.75.
\]

Thus a high pooled value cannot compensate for a failing independent seed block.
The rule yields **12 primary cells** and **15 retained-but-excluded cells** out
of 27.

### Primary cells

| \(\mu\) | \(A_{\rm ref}\) | \(\kappa\) | pooled support | minimum seed-block support |
|---:|---:|---:|---:|---:|
| 0.10 | 0.8 | 6.0 | 0.85 | 0.75 |
| 0.15 | 0.8 | 4.5 | 0.89 | 0.75 |
| 0.15 | 0.8 | 6.0 | 0.93 | 0.90 |
| 0.15 | 1.0 | 6.0 | 0.91 | 0.75 |
| 0.15 | 1.2 | 4.5 | 0.86 | 0.75 |
| 0.20 | 0.8 | 3.0 | 0.87 | 0.80 |
| 0.20 | 0.8 | 4.5 | 0.91 | 0.85 |
| 0.20 | 0.8 | 6.0 | 0.95 | 0.90 |
| 0.20 | 1.0 | 4.5 | 0.95 | 0.90 |
| 0.20 | 1.0 | 6.0 | 0.87 | 0.75 |
| 0.20 | 1.2 | 4.5 | 0.83 | 0.75 |
| 0.20 | 1.2 | 6.0 | 0.94 | 0.85 |

For example, \(\mu=0.20,A_{\rm ref}=1.0,\kappa=3.0\) has pooled support
0.80 but a minimum independent seed-block support of 0.60, and is therefore
not primary eligible.

## How H2/H3 must use this policy

The next branch-conditioned H2/H3 dynamic campaign must:

1. import `primary_analysis_cells()` rather than defining a new hidden subset;
2. write `domain_manifest()` into its artifact manifest;
3. run the 12 eligible cells as primary analysis;
4. retain the 15 excluded cells in a separate feasibility/sensitivity ledger;
5. treat migration only as allele-frequency mixing, not demographic rescue.

This selection rule provides Type S scope control.  It does not establish a
universal mutation rate, and it does not make an H2 timing or H3 fragmentation
claim by itself.
