# Reproducibility guide

This repository is the **mechanistic parent** of the integrated Ecology Letters
submission and the source of a standalone state-separation manuscript. It
contains the theorem-guided ecological framework, the closed finite-model
evidence ledger, and the preregistered fragmentation-gradient sensitivity. The
companion repository,
[`eco-genetic-warning-extensions`](https://github.com/zuizui0223/eco-genetic-warning-extensions),
contains the independently declared directional-transition campaign, owns
full-denominator predictive-warning validity, and assembles the integrated
submission bundle.

The two repositories remain separate provenance units. Reproduction must not pool their trajectories or silently replace the parent scientific commit.

## Canonical scientific state

- repository: `zuizui0223/eco-genetic-criticality`
- canonical scientific commit: `dd8ee379d0d3518194c767d16402042525bc00dc`
- package: `eco-genetic-criticality==0.1.0`
- canonical finite evidence: `docs/final_evidence_ledger.md`
- claim boundaries: `manuscript/claim_evidence_map.md`

Maintenance commits after the canonical scientific commit may improve packaging, documentation, or submission tooling, but they do not alter the closed H1–H3 evidence ledger.

## Reproduction levels

### Level 1 — install and verify the software surface

```bash
git clone https://github.com/zuizui0223/eco-genetic-criticality.git
cd eco-genetic-criticality
git checkout dd8ee379d0d3518194c767d16402042525bc00dc
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev]'
python -m pytest
```

This verifies the theorem, state-transition, censoring, projection, and finite-model contracts covered by the test suite. It does not rerun every expensive simulation campaign.

### Level 2 — verify the frozen scientific release contract

On the submission-maintenance branch or a later release containing the reproducibility tooling:

```bash
python scripts/verify_release_contract.py
python -m pip install build
python -m build
```

The verifier checks that the canonical commit is available, that its load-bearing files have the Git blob identifiers recorded in `reproducibility/release_manifest.json`, and that package identity and evidence-ledger boundaries remain intact.

### Level 3 — rerun finite campaigns

The exact campaign-specific GitHub Actions workflows and original convenience runner paths used to generate H1, H3, H2-R, and the H2-A secondary audit are preserved in the canonical scientific commit. They are intentionally not kept active on the maintenance head because the finite-model campaign is closed.

For exact historical execution, check out the canonical scientific commit. For deliberate robustness reruns against the maintained package, the installed `egc` command provides one consolidated entry point to the retained CLI modules:

```bash
egc --list
egc theorem-boundary --help
egc h1-boundary --help
egc h2r-independent-validation --help
```

For a manuscript reproduction, use the locked numerical summaries and evidence ledgers rather than tuning or rerunning a new parameter search. A full campaign rerun is a robustness exercise and may generate a new Type S evidence set; it must not overwrite the canonical ledger.

## Reproducibility boundaries

The parent result is bounded by the declared symmetric recurrent-mutation, finite trait-recruitment, full-state-transfer, fragmentation, and deterioration closures. In particular:

- H1 and H3 are finite Type S results for the declared source and projection design;
- the fresh fragmentation gradient is finite Type S state-separation evidence
  under its fixed repeated-measures design;
- H2-R is a historical event-conditioned benchmark conditional on one
  trait-loss-only calibrated domain and observed event pairs; it does not
  establish discrimination, specificity, risk separation, or predictive warning
  validity;
- H2-A is retained as a negative robustness audit for fixed absolute thresholds;
- non-events remain censored;
- migration composition bounds are not demographic or functional rescue theorems.

Directional-transition results belong to the companion repository and are not parent evidence.

## Archival release checklist

Before depositing a release in Zenodo or another archive:

1. merge only maintenance changes that do not alter the scientific ledger;
2. run `CI` and `Submission reproducibility` workflows;
3. download the `eco-genetic-criticality-release-bundle` artifact;
4. verify `MANIFEST.sha256`;
5. create an immutable repository release and archive DOI;
6. update the manuscript data/code statement with the final DOI and author-approved citation metadata.

Author order, affiliations, licensing, and CRediT statements require explicit author approval and are therefore not inferred by this repository tooling.
