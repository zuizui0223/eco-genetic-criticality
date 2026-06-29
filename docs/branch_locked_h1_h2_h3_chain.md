# Branch-locked same-replicate H1-H2-H3 audit

This Type S campaign asks a narrower question than the separate H1, H2, and H3
audits: do a finite H1 loop, a genetic warning before realised trait loss, and
the high-branch fragmentation contrast occur in the same seed and replicate?
It does not compose separate theorem certificates into a new theorem.

## Fixed order

For every declared master seed, parameter pair, and replicate:

1. Re-run one-large H1 boundary calibration with the validated 25, 49, and
   97-point grids, padding fraction 0.5, and 30 generations per stage.
2. Use the midpoint of the conservative finest-grid interior
   `[falling recovery bracket upper, rising collapse bracket lower]` as the
   anchor barrier.
3. Initialize canonical low and high H1 starts at that anchor using the same
   derived seed as the calibration replicate.
4. Evaluate H2 on the high-start equal-isolated landscape. H-alpha and H-gamma
   are each compared with realised trait loss.
5. Evaluate H3 on the high start by comparing equal-isolated with one-large,
   and report equal-migrating relative to isolated as allele-frequency mixing.
6. Record same-replicate chain predicates separately for H-alpha, H-gamma, and
   their predeclared union.

## H1 conditioning

A record supplies an H1 anchor only when the calibration has
`resolution_stable_h1_loop_mechanism_supported=True`. A second branch-lock check
requires that the new one-large high and low starts remain separated, that the
high start retains potential high-trait viability, and that the low start lacks
it. Calibration failures and branch-lock failures remain in the output and full
denominator; they are not recoded as H2 or H3 failures.

## H2 censoring

The primary H2 endpoint is high-start equal isolation.

```text
valid pair = warning observed AND realised trait loss observed
censored   = NOT valid pair
warning lead = warning time < realised trait-loss time
```

The output reports lead probability both conditional on valid pairs and across
all H1-branch-locked records. An absent warning or absent trait loss is never
turned into an endpoint-time observation.

## H3 and migration

The high-branch fragmentation predicate requires equal isolation, relative to
one large, to lower terminal interaction, local effective size, and realised
high-trait mass. Migration is reported through FST and mean allele-frequency
modulation relative to isolation. This model has allele-frequency mixing only;
it does not model demographic rescue or recolonisation.

## Running

```bash
python scripts/run_branch_locked_h1_h2_h3_chain.py \
  --profile standard \
  --master-seed 20260630 \
  --master-seed 20260631 \
  --master-seed 20260632 \
  --master-seed 20260633 \
  --master-seed 20260634 \
  --h1-endpoint-padding-fraction 0.5 \
  --h1-stage-generations 30 \
  --h1-barrier-points 25 \
  --h1-barrier-points 49 \
  --h1-barrier-points 97 \
  --h1-maximum-normalized-bracket-width 0.03 \
  --output-dir artifacts/branch_locked_h1_h2_h3 \
  --prefix standard_branch_locked_v1
```

## GitHub Actions

The **Branch-Locked H1 H2 H3 Chain** workflow is manual (`workflow_dispatch`).
For the initial run, use the standard defaults:

```text
profile: standard
master_seeds: 20260630,20260631,20260632,20260633,20260634
h1_endpoint_padding_fraction: 0.5
h1_stage_generations: 30
h1_nested_barrier_points: 25,49,97
h1_maximum_normalized_bracket_width: 0.03
interaction_separation_threshold: 0.05
replicates: blank
allow_full_profile: false
```

It uploads CSV, JSON, and manifest files for 90 days, including partial output
on failure. The run retains every H1 calibration record, H2 censored pair, and
H3 comparison; outcomes do not control inclusion.

CSV is a pair-level ledger. JSON retains every calibration, anchor, branch
outcome, H2 comparison, H3 contrast, migration-mixing contrast, and chain
predicate. The manifest freezes the seed ensemble and the no-selection policy.