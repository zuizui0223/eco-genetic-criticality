# Claim-to-evidence map

This table is the manuscript's wording guardrail. Every substantive claim in
`main_text.md` must map to one row below. It does not add evidence; it records
the appropriate label and the existing source.

| ID | manuscript claim | label | existing support | allowed wording | prohibited wording |
|---|---|---|---|---|---|
| T1 | The canonical sigmoid interaction map has one fixed point when \(K\le4\); for \(K>4\), strict three-fixed-point geometry occurs exactly on an open barrier interval. | T | `docs/canonical_h1_bifurcation.md`, “Exact fixed-point geometry” | “for the declared canonical map”, “exactly” | “all positive-feedback systems” |
| T2 | In that map, the low/high fixed points are locally stable and the middle point is unstable in the strict bistable interval. | T | `docs/canonical_h1_bifurcation.md`, multiplier \(Kq(1-q)\) | “locally stable fixed points” | “globally stable ecological states” |
| C1 | A high-trait component is potentially unavailable on the low branch and viable on the high branch when the declared margin changes sign. | C | `docs/canonical_h1_bifurcation.md`, “Specified-system H1 certificate”; `docs/eco_genetic_hypothesis_program.md`, P1 | “conditional on the declared performance surface” | “the trait necessarily persists/collapses” |
| T3 | Deterministic row-stochastic mixing preserves a common allele-frequency floor before finite sampling. | T | `docs/network_migration_matrix_theory.md`, “Common-floor theorem” | “for the declared migration update” | “migration preserves diversity” |
| T4 | A focal destination has lower bound \(\sum_jM_{ij}b_j\); rescue at a target requires that bound to exceed the target. | T | `docs/network_migration_matrix_theory.md`, “Focal rescue condition” | “frequency rescue condition” | “demographic rescue condition” |
| T5 | Under unbiased post-selection transmission with positive conditional variance, expected heterozygosity declines by twice the transmission variance. | T | `docs/eco_genetic_hypothesis_program.md`, G0 | “under unbiased post-selection transmission” | “genetic diversity always declines” |
| T6 | A contraction bound \(\kappa A M<1\) certifies no bistability for the declared response map. | T | `docs/eco_genetic_hypothesis_program.md`, P0 | “sufficient no-bistability certificate” | “necessary and sufficient bistability criterion” |
| C2 | A patchwise support threshold implies non-additivity of total area under the declared threshold closure. | C | `docs/eco_genetic_hypothesis_program.md`, P2 | “if the closure requires \(A_j>A_c\)” | “fragmentation always harms persistence” |
| C3 | If \(N_e\) increases with interaction state and transmission variance decreases with \(N_e\), low interaction branches erode expected local diversity faster. | C | `docs/eco_genetic_hypothesis_program.md`, G1 | “conditional on the stated monotonicities” | “interaction universally causes genetic erosion” |
| H1-S | A mutation-conditioned high full state retained interaction memory in the declared finite continuation closure. | S | `docs/final_evidence_ledger.md`, H1 | “supported as Type S in the declared closure” | “proved H1 generally” |
| H3-S | Equal isolation lowered interaction, local effective size, and realised high-trait mass conditional on H1 source preparation and conservation-aware projection. | S | `docs/final_evidence_ledger.md`, H1 and H3 | “finite result conditional on source preparation” | “fragmentation always lowers these quantities” |
| H3-G | In the fresh fragmentation gradient, potential high-trait viability was present in 1,037/1,037 one-patch outcomes and absent in 1,037/1,037 outcomes at every tested subdivision, while realised occupancy persisted in approximately 99.6–100% of supported trajectories and realised high-trait mass followed a non-monotone response after its initial loss. | S | `docs/H3_FRAGMENTATION_GRADIENT_SENSITIVITY_RESULTS.md` | “finite state-separation result under the fixed gradient design” | “universal lag”, “all state variables decline monotonically” |
| H2-R | Relative decline in \(H_\alpha\) or \(H_\gamma\) preceded observed realised trait loss in all 35 valid pairs for each of six predeclared endpoints in the selected domain. This is a historical event-conditioned ordering benchmark, not a standalone submission claim. | S / historical non-headline | `docs/final_evidence_ledger.md`, H2-R; `docs/h2_relative_warning_reframing.md` | “all observed valid event pairs in the selected closure” | “all trajectories”, “universal early-warning law”, “discriminative”, “specific”, “risk separating”, “predictively valid warning” |
| H2-A | Fixed thresholds \(H_\alpha,H_\gamma\le0.20\) were not retained as robust rules because the secondary audit observed lags. | S / negative robustness result | `docs/final_evidence_ledger.md`, H2-A; `docs/h2_relative_warning_reframing.md` | “not retained in the selected closure” | “proved false in all systems” |

## Standalone submission claim set

The standalone parent manuscript's headline claim IDs are `T1`, `T2`, `C1`,
`T3`, `T4`, `H1-S`, `H3-S`, and `H3-G`. Its organizing conclusion is that
potential viability, realised occupancy, demographic state, diversity, and
allele persistence are distinct target-dependent states and need not share one
fragmentation response.

`H2-R` and `H2-A` remain in this map only as historical non-headline evidence.
The H2-R observation is conditioned on realised-loss event pairs and does not
establish discrimination, specificity, separation of event from non-event risk,
or predictive warning validity. Those validity questions and their
full-denominator audit are owned by
[`zuizui0223/eco-genetic-warning-extensions`](https://github.com/zuizui0223/eco-genetic-warning-extensions).

## Non-negotiable reporting language for H2

For any warning-versus-trait statement, report all of the following whenever
available:

1. number of attempted trajectories;
2. number with an available source/projection trajectory;
3. number with trait loss observed;
4. valid-pair count;
5. lead, tie, and lag counts; and
6. censoring count or mechanism.

A finite-horizon non-event is not assigned the terminal generation. This rule is
shared by the manuscript and `docs/first_passage_reporting.md`.
