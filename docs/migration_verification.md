# Migration verification

## Status: complete

`eco-genetic-criticality` is the active working repository for the theorem-first H1–H3 eco-genetic criticality program.

The transfer was merged in PR #9 on 2026-06-26. It includes the finite-bin coupled closure, theorem and certificate layers, selected examples and documentation, standalone tests, and a Python 3.10 / 3.11 / 3.12 CI matrix. Subsequent theorem development, including the moving allele-corridor theorem, was merged in PR #11.

The CI run for PR #11 completed successfully on Python 3.10, 3.11, and 3.12. This finalization PR re-runs the same standalone matrix without changing model behavior.

## Scope boundary

This repository contains the active eco-genetic theorem program:

- H1: interaction criticality and potential high-trait viability;
- H2: finite-bin realised trait persistence and genetic-lead certificates;
- H3: fragmentation, migration rescue, and migration erosion conditions.

`microdonta` retains the generic RACH/rule-transition program and historical context. It is not deleted or archived wholesale; Campanula-Izu case-study work, Streamlit tooling, attraction-trait models, and unrelated ABM families remain outside this repository.

After this verification PR is merged, new H1–H3 theory, finite-bin model, simulation, and theorem work should be developed here.