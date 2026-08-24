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

The current H1--H3 finite-model campaign is closed. Its canonical results and limits are recorded in [`docs/final_evidence_ledger.md`](docs/final_evidence_ledger.md).

- **H1:** mutation-conditioned interaction-memory is supported as Type S evidence in the declared finite closure.
- **H3:** conditional on valid H1 full-state transfer, equal isolation lowers interaction, local effective size, and realised high-trait mass as Type S evidence in the declared closure.
- **H2-A:** fixed absolute diversity thresholds \(H_\alpha,H_\gamma\le0.20\) are not retained as a robust canonical warning rule after a no-resimulation secondary audit found mixed lead/lag ordering.
- **H2-R:** baseline-relative \(H_\alpha/H_\gamma\) erosion precedes observed realised trait loss in one calibration-selected deterioration configuration; this is conditional Type S evidence, not a universal rule.

## Submission role

This repository is the **mechanistic parent** of the integrated Ecology Letters submission. It provides the theorem-guided interaction and fragmentation framework, the closed finite H1/H3 ledger, and the conditional symmetric-warning benchmark.

The companion repository, [`eco-genetic-warning-extensions`](https://github.com/zuizui0223/eco-genetic-warning-extensions), contains the independently declared directional-transition campaign and is the submission-bundle orchestrator. The code bases and evidence ledgers remain separate; only the manuscript argument and release package are integrated.

## Reproduce and package

Start with [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md). The machine-readable scientific lock is [`reproducibility/release_manifest.json`](reproducibility/release_manifest.json).

A lightweight verification is:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev,reproducibility]'
python -m pytest
python scripts/verify_release_contract.py
python -m build
```

The `Submission reproducibility` workflow additionally smoke-tests the built wheel and uploads a checksummed `eco-genetic-criticality-release-bundle`. The canonical scientific commit remains `dd8ee379d0d3518194c767d16402042525bc00dc`; later packaging maintenance does not revise its evidence.

## Maintenance command surface

The default branch is now a maintenance/release surface rather than an active simulation campaign. Historical campaign-specific GitHub Actions workflows remain immutable in the canonical scientific commit instead of staying enabled on current `main`.

After installing the package, use the single `egc` entry point for the retained command-line interfaces:

```bash
egc --list
egc theorem-boundary --help
egc h2r-independent-validation --help
```

The former `scripts/run_*.py` convenience wrappers are intentionally removed from the maintenance head. Historical protocol documents preserve the original commands as provenance; check out the canonical scientific commit when exact historical execution is required. New maintenance work should not add one-off runner wrappers or new campaign workflows to this parent repository.

Only two Actions workflows are active on the maintenance head: ordinary `CI` and `Submission reproducibility`.

## Manuscript synthesis

The paper-facing theory draft is in [`manuscript/`](manuscript/). It separates exact theorems, closure-conditional results, dynamic hypotheses, and finite Type S results; it neither alters the final evidence ledger nor introduces new simulations. The entry point is [`manuscript/README.md`](manuscript/README.md).

## Current model layers

- `canonical_h1_bifurcation.py` gives the specified-system H1 certificate for the one-state logistic reduction: strict bistability, branch stability, and high-trait margin change.
- `first_passage_reporting.py` and `censoring_aware_phase_diagram.py` keep H2 warning lead probabilities, valid-pair denominators, and censored replicates distinct.
- `network_migration_matrix_theory.py` gives H3 allele-floor and focal-rescue bounds for arbitrary network mixing matrices.
- `network_h3_lifecycle.py` and `network_h3_experiments.py` add a separate finite-population H3 closure with individual dispersal, extinction, rescue, recolonisation, realised trait abundance, allele copies, and replicated event summaries.

The H3 lifecycle is a declared stochastic model, not a universal claim that connectivity is beneficial. It can represent migration rescue, migration erosion, or no material effect under different declared kernels and life-cycle parameters. See `docs/h3_extinction_recolonisation_lifecycle.md`.

## Repository ownership boundary

This repository is the sole owner of eco-genetic criticality code and evidence.
Conversely, RACH N1–N4 channel identifiability, theorem-to-model projection, and
next-observation design are owned by
[`zuizui0223/microdonta`](https://github.com/zuizui0223/microdonta).

The former local copy of `theorem_projection_ledger.py` described microdonta
backends and Campanula rather than the H1–H3 criticality programme. It has been
removed together with its mirrored test. Cross-program context must use a link or
an explicit adapter, not a second implementation.

## Scope

This repository closes the active H1--H3 theorem and finite-bin closure program migrated from `microdonta`. Further biological closures, mutation models, threshold choices, or deterioration schedules should be developed as separately declared extensions rather than silently revising the final evidence ledger. It excludes generic RACH/rule-transition methods, Campanula-Izu case-study work, Streamlit tooling, attraction-trait models, and unrelated ABM families. Historical migration and campaign provenance remain available in Git history and the canonical scientific commit.
