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

## Current model layers

- `canonical_h1_bifurcation.py` gives the specified-system H1 certificate for the one-state logistic reduction: strict bistability, branch stability, and high-trait margin change.
- `first_passage_reporting.py` and `censoring_aware_phase_diagram.py` keep H2 warning lead probabilities, valid-pair denominators, and censored replicates distinct.
- `network_migration_matrix_theory.py` gives H3 allele-floor and focal-rescue bounds for arbitrary network mixing matrices.
- `network_h3_lifecycle.py` and `network_h3_experiments.py` add a separate finite-population H3 closure with individual dispersal, extinction, rescue, recolonisation, realised trait abundance, allele copies, and replicated event summaries.

The H3 lifecycle is a declared stochastic model, not a universal claim that connectivity is beneficial. It can represent migration rescue, migration erosion, or no material effect under different declared kernels and life-cycle parameters. See `docs/h3_extinction_recolonisation_lifecycle.md`.

## Scope

This is the active H1–H3 theorem and finite-bin closure program migrated from `microdonta`. It excludes generic RACH/rule-transition methods, Campanula-Izu case-study work, Streamlit tooling, attraction-trait models, and unrelated ABM families. See `MIGRATION_MANIFEST.md`.