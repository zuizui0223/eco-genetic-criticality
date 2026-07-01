# Eco-genetic criticality

A theorem-first research repository for finite-population eco-genetic criticality.

## Questions

- **H1:** When can interaction feedback alter potential high-trait viability?
- **H2:** Under what declared conditions can genetic warning precede realised high-trait loss?
- **H3:** At fixed total area, when do isolation, migration rescue, and migration erosion produce different trait and genetic outcomes?

## Research architecture

```text
mathematical theorem
-> declared model projection
-> finite-population closure
-> simulation robustness test
-> empirical measurement design
```

Potential viability, realised trait occupancy, allele persistence, and genetic diversity are distinct states.

## Final finite-model status

The current H1--H3 finite-model campaign is closed. Its canonical results and
limits are recorded in [`docs/final_evidence_ledger.md`](docs/final_evidence_ledger.md).

- **H1:** mutation-conditioned interaction-memory is supported as Type S evidence
  in the declared finite closure.
- **H3:** conditional on valid H1 full-state transfer, equal isolation lowers
  interaction, local effective size, and realised high-trait mass as Type S
  evidence in the declared closure.
- **H2-A:** fixed absolute diversity thresholds \(H_\alpha,H_\gamma\le0.20\)
  are not retained as a robust canonical warning rule after a no-resimulation
  secondary audit found mixed lead/lag ordering.
- **H2-R:** baseline-relative \(H_\alpha/H_\gamma\) erosion precedes observed
  realised trait loss in one calibration-selected deterioration configuration;
  this is conditional Type S evidence, not a universal rule.

## Current model layers

- `canonical_h1_bifurcation.py` gives the specified-system H1 certificate for the one-state logistic reduction: strict bistability, branch stability, and high-trait margin change.
- `first_passage_reporting.py` and `censoring_aware_phase_diagram.py` keep H2 warning lead probabilities, valid-pair denominators, and censored replicates distinct.
- `network_migration_matrix_theory.py` gives H3 allele-floor and focal-rescue bounds for arbitrary network mixing matrices.
- `network_h3_lifecycle.py` and `network_h3_experiments.py` add a separate finite-population H3 closure with individual dispersal, extinction, rescue, recolonisation, realised trait abundance, allele copies, and replicated event summaries.

The H3 lifecycle is a declared stochastic model, not a universal claim that connectivity is beneficial. It can represent migration rescue, migration erosion, or no material effect under different declared kernels and life-cycle parameters. See `docs/h3_extinction_recolonisation_lifecycle.md`.

## Scope

This repository closes the active H1--H3 theorem and finite-bin closure program migrated from `microdonta`. Further biological closures, mutation models, threshold choices, or deterioration schedules should be developed as separately declared extensions rather than silently revising the final evidence ledger. It excludes generic RACH/rule-transition methods, Campanula-Izu case-study work, Streamlit tooling, attraction-trait models, and unrelated ABM families. See `MIGRATION_MANIFEST.md`.
