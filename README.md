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
-> state-representation boundary
-> empirical measurement design in the extension repository
```

Potential viability, realised trait occupancy, allele persistence, genetic diversity, and realised ecological function are distinct states.

## Final finite-model status

The active H1–H3 finite-model campaign is closed. Its canonical results and limits are recorded in [`docs/final_evidence_ledger.md`](docs/final_evidence_ledger.md).

- **H1:** mutation-conditioned interaction memory is supported as finite Type S evidence in the declared closure.
- **H3 / fragmentation:** conditional on valid H1 full-state transfer, equal isolation lowers interaction, local effective size, and realised high-trait mass before demographic disappearance.
- **H2-A:** fixed absolute diversity thresholds `H_alpha,H_gamma <= 0.20` are not retained as a robust canonical warning rule after the no-resimulation audit found mixed lead/lag ordering.
- **H2-R:** baseline-relative `H_alpha/H_gamma` erosion precedes observed realised trait loss in the inherited calibrated symmetric domain. The extension later prospectively reproduced all six relative warning orderings in an independent fresh-seed ensemble under the same frozen loss-generating state, without turning them into universal thresholds.

## State-sufficiency result

Under the declared simulator closure, the complete explicit present state together with future forcing and the stochastic law is Markov/future-sufficient. This is a theorem about the declared model representation, not a claim that any usual ecological summary is sufficient in nature.

A constructive alignment audit shows why the distinction matters. Two states can share census, interaction and allele-frequency marginals, realised trait-bin state, `H_alpha`, `H_gamma`, and `F_ST` while differing in patchwise cross-layer alignment. Their exact next interaction transition can then differ substantially. The fixed long-horizon campaign did not establish a directional loss-incidence effect of alignment, so the result is a **representation boundary**, not a universal alignment-risk rule.

## Submission role

This repository is the **mechanistic parent** of the integrated Letter manuscript. It owns the theorem-guided interaction/fragmentation framework, the closed finite H1/H3 evidence ledger, state-sufficiency certificates, and the inherited conditional warning benchmark.

The companion repository, [`eco-genetic-warning-extensions`](https://github.com/zuizui0223/eco-genetic-warning-extensions), is now the **condition-recovery, warning-replication/portability, state-representation, natural state-sufficiency, integrated-manuscript, and submission-bundle repository**. Its empirical programme tests whether candidate ecological state variables predict their downstream endpoints before asking whether geography or fragmentation history adds residual information.

The repositories remain separate computational provenance units. They are integrated at the argument, manuscript, reproducibility-contract, and release-package levels rather than by merging code or evidence ledgers.

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

The `Submission reproducibility` workflow additionally smoke-tests the built wheel and uploads a checksummed release bundle. The canonical scientific commit remains `dd8ee379d0d3518194c767d16402042525bc00dc`; later packaging/README maintenance does not revise its scientific evidence.

The package version is currently `0.1.0`, matching the extension repository. Final release/citation metadata are coordinated from the extension release-readiness ledger and require explicit author approval before immutable tags or archive DOI creation.

## Maintenance command surface

The default branch is a maintenance/release surface rather than an active simulation campaign. Historical campaign-specific GitHub Actions workflows remain immutable in the canonical scientific commit instead of staying enabled on current `main`.

After installing the package, use the single `egc` entry point for the retained command-line interfaces:

```bash
egc --list
egc theorem-boundary --help
egc h2r-independent-validation --help
```

The former `scripts/run_*.py` convenience wrappers are intentionally removed from the maintenance head. Historical protocol documents preserve the original commands as provenance; check out the canonical scientific commit when exact historical execution is required. New maintenance work should not add one-off runner wrappers or new campaign workflows to this parent repository.

Only ordinary `CI` and `Submission reproducibility` workflows should remain active on the maintenance head.

## Manuscript synthesis

The parent paper-facing theory material is in [`manuscript/`](manuscript/). It separates exact theorems, closure-conditional results, dynamic hypotheses, and finite Type S results. The integrated manuscript and the current natural-data synthesis live in the extension repository.

## Current model layers

- `canonical_h1_bifurcation.py` gives the specified-system H1 certificate for the one-state logistic reduction: strict bistability, branch stability, and high-trait margin change.
- `first_passage_reporting.py` and `censoring_aware_phase_diagram.py` keep H2 warning lead probabilities, valid-pair denominators, and censored replicates distinct.
- `network_migration_matrix_theory.py` gives H3 allele-floor and focal-rescue bounds for arbitrary network mixing matrices.
- `network_h3_lifecycle.py` and `network_h3_experiments.py` add a separate finite-population H3 closure with individual dispersal, extinction, rescue, recolonisation, realised trait abundance, allele copies, and replicated event summaries.

The H3 lifecycle is a declared stochastic model, not a universal claim that connectivity is beneficial. It can represent migration rescue, migration erosion, or no material effect under different declared kernels and life-cycle parameters. See `docs/h3_extinction_recolonisation_lifecycle.md`.

## Scope and stop rule

This repository closes the active H1–H3 theorem and finite-bin closure programme migrated from `microdonta`. Further biological closures, mutation models, threshold choices, or deterioration schedules must be separately declared extensions rather than silent revisions of the final evidence ledger.

Do not reopen simulator parameter/seed tuning merely to obtain a preferred empirical or warning result. New biological evidence belongs in prospectively declared extension work; parent maintenance should preserve the canonical scientific commit and release provenance.
