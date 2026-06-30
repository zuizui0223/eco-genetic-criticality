# Finite H1 full-state fragment projection audit

## Why this bridge is necessary

The H1 full-state hold audit showed that a high finite branch persists at a
shared interior barrier when its complete state is carried:

```text
population
interaction q
high-allele frequency
realised trait-bin abundance
```

To test fragmentation, that one-large state must be transferred into a
multi-patch landscape without silently changing the amount of population,
genetic state, or realised trait material. This audit fixes that construction
before any H2/H3 outcome is interpreted.

## Declared projection rule

Let the source state have total population \(N\), trait-bin totals \(C_k\),
population-weighted interaction \(\bar q\), population-weighted high-allele
frequency \(\bar p\), and total area \(A\). Let target fragment areas be
\(a_1,\ldots,a_m\), with \(\sum a_i=A\).

The target receives:

1. **Population:** positive integer counts \(N_i\) proportional to target area,
   allocated with deterministic largest-remainder rounding, satisfying
   \(\sum_i N_i=N\).
2. **Trait abundance:** an integer matrix \(C_{ik}\) with exact margins

   \[
   \sum_k C_{ik}=N_i,\qquad \sum_i C_{ik}=C_k.
   \]

   This is a deterministic contingency allocation, not independent bin-wise
   rounding. Therefore both patch population and each trait-bin total are exact.
3. **Interaction:** copy the intensive source mean to every patch,
   \(q_i=\bar q\).
4. **Allele frequency:** copy the intensive source mean to every patch,
   \(p_i=\bar p\). Because total population is conserved, the target's
   population-weighted frequency remains exactly \(\bar p\).

Source and target total area must match. The construction does not apply
migration; equal-migrating differs only later through scenario dynamics.

## What is tested

For each declared seed, parameter pair, and replicate, the audit:

1. reproduces a 25 → 49 → 97-point H1 calibration;
2. replays the high route to a declared interior grid anchor;
3. runs the full-state high hold for 30 generations;
4. projects that source into one-large, equal-isolated, and equal-migrating
   parameter templates; and
5. asks the simulator for its actual generation-0 snapshot, then checks all
   conservation statements there.

The output keeps every seed-replicate. A record without a validated H1 full
state is unavailable, not removed.

## Interpretation

Passing establishes only that the initial state bridge is well-defined and
conservation-preserving under this declared finite closure. It does not say
fragmentation has changed interaction, traits, effective size, genetic
variation, or warning order. Those are the next H2/H3 dynamic comparisons.

## Standard run

```bash
python scripts/run_finite_h1_fragment_projection_audit.py \
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
  --output-dir artifacts/h1_fragment_projection \
  --prefix standard_fragment_projection_v1
```

CSV is a pair-level ledger. JSON records the source high state and each target
scenario's initial population, q, high-allele frequency, and all invariants.
The manifest records the allocation algorithms and no-selection policy.
