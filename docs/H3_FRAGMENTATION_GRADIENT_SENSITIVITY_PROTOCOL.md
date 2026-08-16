# H3 fragmentation-gradient sensitivity protocol

## Status

**Post-review supplementary sensitivity analysis, declared before inspecting any new gradient outcome.**

This campaign does not alter the closed H1/H3 evidence ledger or replace the canonical one-large versus equal-isolated result. It addresses a manuscript-review question: whether the very large H3 endpoint contrast could be an artefact of comparing only one large patch with four equal isolated patches.

The canonical scientific state remains commit `dd8ee379d0d3518194c767d16402042525bc00dc`. The new campaign reuses that model closure and the already frozen 12 primary mutation-H1 cells, but uses a fresh seed family and new simulations.

### Pre-outcome maintenance note

The first execution attempt stopped before producing an admissible complete gradient because the existing private parent helper `mutation_primary_h1_h2_h3_chain._prepare_mutation_high_state` still called `canonical_h1_certificate` with the function's former positional API. The current certificate API is keyword-only. Before inspecting any fragmentation-gradient outcome, that call was changed to the exactly equivalent keyword arguments (`feedback_strength`, `area`, `area_reference`, `barrier`, `trait_parameters`). The parent-chain and new-gradient targeted tests passed after this compatibility-only change. No equation, parameter, seed, source-preparation rule, landscape design, or outcome definition was changed. The failed first execution is not evidence; the first complete run after this maintenance fix is the authoritative gradient campaign.

## Question

At fixed total area, what happens to interaction state, local effective size, and realised high-trait mass as the same H1-prepared high full state is partitioned into progressively more **isolated** equal patches?

The analysis is intended to distinguish a genuine fragmentation gradient or threshold from a degenerate endpoint-only contrast.

## Fixed design

- total area: `4.0`;
- primary cells: the 12 cells returned by the frozen `primary_analysis_cells()` definition;
- mutation closure: symmetric recurrent mutation, exactly as in the parent H1/H3 campaign;
- source preparation: finite H1 boundary-resolution audit followed by the same high-state replay and 30-generation full-state hold used by `mutation_primary_h1_h2_h3_chain`;
- outcome horizon after projection: 30 generations;
- migration: `0.0` in every gradient landscape;
- patch counts: **1, 2, 3, 4, 6, 8, 12, 16**;
- patch area: `4.0 / patch_count`;
- fresh master seeds: **20260820, 20260821, 20260822, 20260823, 20260824**;
- source replicates: 20 per master seed and primary cell;
- attempted source replicates per primary cell: 100;
- repeated-measures unit: one prepared source replicate projected across all eight patch counts.

No outcome from the historical H1/H3 artifact is used as a new gradient observation. The historical result is used only to define the frozen primary-cell set and the model/protocol identity.

## Analytical coordinate

For each primary cell and patch count, record the canonical feedback coordinate

\[
K_{\mathrm{fragment}} = \frac{\kappa A_{\mathrm{patch}}}{A_{\mathrm{ref}}}
= \frac{\kappa(4/n)}{A_{\mathrm{ref}}}.
\]

For the stated canonical interaction map, `K = 4` is the analytical boundary below which the three-fixed-point structure is unavailable. The finite gradient is **not** assumed to collapse exactly at `K=4`; the coordinate is recorded only to connect the finite sensitivity to the theorem-guided mechanism.

## Primary outcomes

For every projection-supported patch-count outcome record:

1. mean final interaction across patches;
2. mean final local effective size across patches;
3. realised high-trait mass mean;
4. potential high-trait viability;
5. realised high-trait persistence;
6. `H_alpha` and `H_gamma` as descriptive genetic states.

The primary gradient display uses the first three outcomes because those are the three quantities in the canonical H3 fragmentation claim.

## Paired estimands

Within each source replicate, patch count 1 is the paired reference. For each larger patch count and each primary metric, calculate

\[
R_x(n)=1-\frac{x_n}{x_{n=1}}.
\]

If a reference value is zero, the paired reduction for that metric is missing rather than imputed.

For each primary cell and patch count report:

- projection-supported denominator;
- median and interquartile range of each raw metric;
- median and interquartile range of each paired reduction relative to patch count 1;
- probability that **all three** primary metrics are below their paired patch-count-1 values;
- potential-high-trait viability probability;
- realised-high-trait persistence probability.

Pooled summaries across the 12 cells are descriptive only; cell-specific curves remain visible in the supplementary figure.

## Interpretation rules

- This is new finite Type S sensitivity evidence, not a theorem or empirical effect size.
- It may support or weaken the claim that the canonical H3 endpoint contrast represents a broader fragmentation response.
- It must not overwrite the historical 1,055-replicate H3 ledger.
- Failure to prepare an H1 high full state remains in the denominator and is never treated as a fragmentation failure.
- Patch-count outcomes from the same prepared source are repeated measures and are not independent replicates.
- No warning endpoint is selected, tuned, or evaluated in this campaign.
- `K=4` is an analytical reference for the specified canonical map, not a fitted finite threshold.

## Planned supplementary outputs

- `h3_fragmentation_gradient_records.csv`: one row per attempted source × patch count;
- `h3_fragmentation_gradient_cell_summary.csv`: cell-level gradient summaries;
- `h3_fragmentation_gradient_pooled_summary.csv`: pooled descriptive summaries;
- `figure_s_fragmentation_gradient.svg`: three-panel paired-gradient figure for interaction, local effective size, and realised high-trait mass;
- self-describing JSON metadata with seed family, patch counts, scientific commit, and protocol boundary.

The manuscript will cite this analysis only after the full declared grid has completed and the artifact has passed denominator and reproducibility checks.
