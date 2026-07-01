# Figure and table plan

This plan uses only existing equations, committed evidence ledgers, and already
completed workflow artifacts. It does not authorize a new parameter scan or
simulation campaign.

## Main figures

### Figure 1 — Claim architecture from interaction feedback to genetic evidence

**Purpose.** Establish the manuscript's key distinction between the causal chain
and its evidential layers.

**Panel A: conceptual chain**

```text
patch size
  -> interaction state q
  -> potential high-trait viability
  -> realised high-trait occupancy / local N_e
  -> H_alpha, H_gamma, and first-passage events
```

**Panel B: labels on arrows**

- exact canonical-map theorem for interaction branch geometry;
- conditional trait-margin implication;
- finite trait–allele closure for realised occupancy and local \(N_e\);
- censoring-aware finite assessment for warning order.

**Source.** Deterministic redraw from `docs/eco_genetic_hypothesis_program.md`
and `manuscript/claim_evidence_map.md`.

**Do not imply.** That every arrow is proved by a single theorem.

### Figure 2 — Canonical interaction-map geometry

**Purpose.** Visualize the exact H1 theorem for the specified sigmoid map.

**Panel A.** Plot \(F(q)=\operatorname{logit}(q)-Kq+\kappa\theta\) under a
representative \(K>4\) configuration, showing three roots.

**Panel B.** Plot the two stable branches and the unstable middle branch over
\(\theta\), marking \(\theta_-\) and \(\theta_+\).

**Panel C.** Show the sign of the declared high-trait viability margin on the
low and high branches as a schematic condition, not a new empirical estimate.

**Source.** `canonical_h1_bifurcation.py` and
`docs/canonical_h1_bifurcation.md`; the existing \(\kappa=8\),
\(A/A_{\rm ref}=1\) hysteresis illustration may be reused.

**Caption boundary.** “Exact for the stated one-state logistic reduction.”

### Figure 3 — What migration theorems do and do not establish

**Purpose.** Separate allele-frequency mixing from demographic or trait rescue.

**Panel A.** Directed network with source frequencies \(p_j\), weights
\(M_{ij}\), and target patch \(i\).

**Panel B.** Common-floor inequality \(p'_i\ge p_{\min}\).

**Panel C.** Focal bound \(p'_i\ge\sum_jM_{ij}b_j\), with a target line
\(p_{\rm target}\).

**Source.** Deterministic schematic from `docs/network_migration_matrix_theory.md`.

**Caption boundary.** The figure does not depict individual movement, abundance,
extinction, or recolonisation.

### Figure 4 — Predeclared finite-warning design

**Purpose.** Make the anti-selection logic auditable.

```text
H1-prepared full state
-> equal-isolated projection
-> trait-loss-only ramp-and-hold calibration
-> one selected domain
-> fresh seeds
-> six relative-warning endpoints
-> censored first-passage report
```

**Source.** `docs/h2_relative_warning_reframing.md` and the final ledger.

**Required annotations.**

- calibration sees no \(H_\alpha\), \(H_\gamma\), warning time, or lead time;
- fresh validation seeds differ from calibration seeds;
- non-events remain censored;
- absolute threshold audit is secondary, not a selection stage.

### Figure 5 — Relative versus absolute diversity-warning results

**Purpose.** Present the finite H2 result without overstating it.

**Panel A.** Flow/count plot: 100 attempted sources → 83 available trajectories
→ 35 observed trait-loss trajectories + 48 right-censored trajectories.

**Panel B.** Six relative endpoints: each 35 valid pairs, 35 leads, 0 ties,
0 lags.

**Panel C.** Fixed absolute audit: \(H_\alpha\le0.20\) shows 14 leads / 6 lags
among 20 valid pairs; \(H_\gamma\le0.20\) shows 8 leads / 8 lags among 16 valid
pairs.

**Source.** `docs/final_evidence_ledger.md` and the committed H2-A audit CSVs.

**Caption boundary.** The relative result is conditional Type S evidence in one
selected configuration and for observed event pairs.

## Main tables

### Table 1 — Claim hierarchy and scope

Columns: identifier, statement, label (T/C/H/S), assumptions, manuscript
section, and forbidden overclaim. Source: `manuscript/claim_evidence_map.md`.

### Table 2 — State variables and non-equivalences

Rows: interaction state, potential viability, realised high-trait occupancy,
allele persistence, \(H_\alpha\), \(H_\gamma\), \(F_{ST}\). Columns: definition,
where introduced, and what it is not equivalent to.

### Table 3 — H2 calibration and independent validation ledger

Rows: calibration selection, source preparation, projection support, trait loss,
relative-warning endpoints, absolute-threshold audit. Include numerical counts
and the corresponding denominator.

## Supplementary figures

- **Figure S1:** Derivative geometry \(F'(q)\) and the \(K=4\) threshold.
- **Figure S2:** Hysteresis path for the existing \(\kappa=8\),
  \(A/A_{\rm ref}=1\) example; label it as a numerical confirmation of the
  analytic branch geometry rather than an independent finite result.
- **Figure S3:** First-passage bookkeeping with lead, tie, lag, and censored
  outcomes.
- **Figure S4:** Seed-block denominators and observed trait-loss counts for the
  H2-R independent validation.

## Implementation order

1. Draw Figures 1 and 4 as conceptual diagrams from existing definitions.
2. Generate Figure 2 directly from the canonical analytic functions.
3. Generate Figure 5 from the existing committed H2 ledger CSVs and artifacts.
4. Do not create a new “best warning metric” figure: the six endpoints are not
   independent and the completed campaign did not select a winner.
