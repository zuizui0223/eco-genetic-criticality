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
- finite fragmentation-gradient assessment of distinct state responses.

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

### Figure 4 — Preregistered fragmentation-gradient design

**Purpose.** Make source preservation and the repeated-measures fragmentation
contrast auditable.

```text
H1-prepared full state
-> one, two, three, four, six, eight, twelve, or sixteen equal isolated patches
-> same prepared source represented at every patch count
-> interaction / local N_e / realised high-trait mass
-> potential viability versus realised occupancy
```

**Source.** `docs/H3_FRAGMENTATION_GRADIENT_SENSITIVITY_PROTOCOL.md` and
`docs/H3_FRAGMENTATION_GRADIENT_SENSITIVITY_RESULTS.md`.

**Required annotations.**

- the source full state is held fixed within each repeated-measures comparison;
- master seeds and patch counts were declared before outcome inspection;
- 1,037 sources support every patch-count projection;
- the first stopped execution is not evidence.

### Figure 5 — State separation across the fragmentation gradient

**Purpose.** Show why one eco-genetic summary cannot represent all target states.

**Panel A.** Retained interaction and local effective size versus patch count,
showing continued declines after the first subdivision.

**Panel B.** Retained realised high-trait mass versus patch count, showing the
initial loss followed by partial recovery.

**Panel C.** Potential viability and realised occupancy: viability changes from
1,037/1,037 at one patch to 0/1,037 at every subdivision, while occupancy
persists in approximately 99.6–100% of supported trajectories at generation 30.

**Source.** The verified fragmentation-gradient publication artifact summarized
in `docs/H3_FRAGMENTATION_GRADIENT_SENSITIVITY_RESULTS.md`.

**Caption boundary.** These are finite Type S responses under one declared
closure, not a universal fragmentation dose–response or lag law.

## Main tables

### Table 1 — Claim hierarchy and scope

Columns: identifier, statement, label (T/C/H/S), assumptions, manuscript
section, and forbidden overclaim. Source: `manuscript/claim_evidence_map.md`.

### Table 2 — State variables and non-equivalences

Rows: interaction state, potential viability, realised high-trait occupancy,
allele persistence, \(H_\alpha\), \(H_\gamma\), \(F_{ST}\). Columns: definition,
where introduced, and what it is not equivalent to.

### Table 3 — Fragmentation-gradient state ledger

Rows: patch count, interaction retained, local effective-size retained, realised
high-trait-mass retained, potential viability, and realised occupancy. Include
the repeated-measures denominator and the single retained exception at eight
patches.

## Supplementary figures

- **Figure S1:** Derivative geometry \(F'(q)\) and the \(K=4\) threshold.
- **Figure S2:** Hysteresis path for the existing \(\kappa=8\),
  \(A/A_{\rm ref}=1\) example; label it as a numerical confirmation of the
  analytic branch geometry rather than an independent finite result.
- **Figure S3:** First-passage bookkeeping with lead, tie, lag, and censored
  outcomes for the historical H2-R benchmark. Label it non-headline and
  event-conditioned.
- **Figure S4:** Historical H2-R denominators: 100 attempted, 83 available, 35
  event trajectories, and 48 non-event trajectories. State explicitly that the
  parent ordering result does not establish discrimination, specificity, risk
  separation, or predictive warning validity.
- **Figure S5:** Historical fixed-threshold H2-A lead/tie/lag audit.

## Implementation order

1. Draw Figures 1 and 4 as conceptual diagrams from existing definitions.
2. Generate Figure 2 directly from the canonical analytic functions.
3. Generate Figure 5 from the verified fragmentation-gradient publication
   artifact; do not rerun or retune the gradient.
4. Keep H2-R/H2-A displays supplementary and historical. Do not create a new
   “best warning metric” figure or infer predictive validity from event-conditioned
   ordering.
