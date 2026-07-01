# Interaction thresholds, fragmentation, and conditional genetic warning: a theorem-guided finite-model framework

**Draft v0.1 — theory and existing-results synthesis**

## Abstract

Ecological fragmentation can alter interaction support, trait persistence, and
genetic variation, but these processes are often summarized as though they were
a single mechanism. We present a theorem-guided framework that separates the
analytical statements available for specified ecological maps from conditional
consequences of a declared life cycle and from finite stochastic simulation
results. For a canonical one-state positive-feedback map, we derive an exact
criterion for strict bistability and identify the barrier interval in which low
and high locally stable interaction states coexist. Coupling this branch geometry
to a declared high-trait viability margin gives a specified-system result: the
potential high-trait component is unavailable on the low interaction branch but
viable on the high branch. For network migration, we establish a common-floor
bound and a sharp focal-rescue bound for any destination-by-source row-stochastic
mixing matrix. We then formulate a finite trait–allele closure in which potential
viability, realised trait occupancy, local diversity, metapopulation diversity,
and allele persistence are distinct states. In that closure, equal isolation from
an H1-prepared full state lowered interaction, local effective size, and realised
high-trait mass. A predeclared, trait-loss-only calibration identified one
ramp-and-hold deterioration configuration for an independent genetic-warning
assessment. In the selected configuration, baseline-relative declines in
within-patch and metapopulation diversity preceded every observed realised
trait-loss event across six predeclared endpoints; trajectories without loss in
the finite horizon remained censored. Fixed absolute diversity thresholds showed
mixed lead–lag ordering and were therefore not retained as a robust warning rule.
The framework clarifies which conclusions are theorems, which are conditional on
a closure, and which are finite Type S evidence.

**Keywords:** positive feedback; bistability; fragmentation; genetic diversity;
early warning; first-passage time; censoring; eco-evolutionary dynamics

---

## 1. Introduction

Abrupt ecological change is frequently discussed in terms of interaction loss,
trait loss, demographic decline, and genetic erosion. These are related but not
interchangeable processes. An interaction state may make a trait potentially
viable without guaranteeing its realised occupancy; a change in local diversity
need not equal a change in metapopulation diversity; and a finite trajectory that
does not reach an event within an observation window is censored rather than
negative evidence for that event.

The central aim of this study is therefore methodological as much as ecological:
to construct an explicit logical bridge from interaction feedback to trait and
genetic outcomes without allowing results from a finite simulation to silently
become a general theorem. We organize the argument around the causal chain

\[
\text{patch size}
\rightarrow \text{interaction intensity}
\rightarrow \text{trait-space topology}
\rightarrow \text{local effective size}
\rightarrow \text{genetic diversity}.
\]

The chain is a hypothesis program, not a theorem in itself. Each arrow requires
an equation, a life-cycle map, or an explicitly declared finite closure. The
benefit of this separation is that it identifies where exact mathematics ends and
where model-dependent inference begins.

We make four contributions using material already established in the associated
repository. First, we state the exact fixed-point geometry of a canonical
positive-feedback reduction and show how it supports a branch-dependent
high-trait viability result. Second, we state migration bounds that are valid for
arbitrary row-stochastic network mixing matrices. Third, we specify a finite
trait–allele experiment in which potential viability, realised occupancy,
allele-frequency persistence, and diversity metrics are tracked separately.
Fourth, we formulate warning assessment as a censoring-aware first-passage
comparison and report a predeclared finite validation of relative genetic
warnings. The study does not claim a universal biological early-warning law or an
estimate of any real mutation parameter.

**[Literature insertion point: ecological critical transitions, positive
frequency dependence, fragmentation genetics, and early-warning theory.]**

---

## 2. Claim taxonomy and state separation

### 2.1 Four evidence labels

We distinguish four kinds of statements throughout. A **Type T** result is a
mathematical theorem under explicit assumptions. A **Type C** result is
conditional on a declared ecological or life-cycle closure. A **Type H** result
is a dynamic hypothesis. A **Type S** result is numerical evidence obtained from
a declared finite simulation and is neither a theorem nor a biological parameter
estimate.

This taxonomy is substantive rather than cosmetic. It prevents three common
collapses: treating a canonical reduction as a theorem for every ecological
system; treating an observed event order in finite time as a universal causal
rule; and treating an unobserved event in a finite horizon as evidence that it
cannot occur.

### 2.2 Distinct state variables

The model tracks five conceptually distinct objects:

1. interaction state \(q\), which controls feedback-mediated support;
2. potential high-trait viability, represented by a performance-margin
   calculation;
3. realised high-trait abundance or occupancy in the finite recruitment model;
4. allele frequencies and their persistence; and
5. genetic summaries, including local diversity \(H_\alpha\), metapopulation
   diversity \(H_\gamma\), and differentiation where relevant.

No implication among these states is assumed without a stated update rule. In
particular, branch-dependent potential viability is not identical to realised
trait occupancy, and a decrease in \(H_\alpha\) is not identical to a decrease in
\(H_\gamma\).

---

## 3. Analytical result I: interaction feedback and branch-dependent trait viability

### 3.1 Canonical interaction map

Consider the one-state interaction reduction

\[
q_{t+1}=f(q_t)=\operatorname{sigmoid}\!\left[
\kappa\left(\frac{A}{A_{\rm ref}}q_t-\theta\right)
\right],
\qquad q_t\in[0,1],
\]

where \(A\) is the focal patch area, \(A_{\rm ref}\) is a reference area,
\(\kappa>0\) controls feedback steepness, and \(\theta\) is an interaction
barrier. Let

\[
K=\kappa A/A_{\rm ref}.
\]

A fixed point solves

\[
F(q)=\operatorname{logit}(q)-Kq+\kappa\theta=0.
\]

### 3.2 Exact branch geometry

**Theorem 1 (canonical fixed-point geometry; Type T for the stated map).**
If \(K\le4\), the canonical map has a unique fixed point. If \(K>4\), define

\[
q_- = \frac{1-\sqrt{1-4/K}}{2},
\qquad
q_+ = \frac{1+\sqrt{1-4/K}}{2},
\]

and

\[
\theta_- = \frac{Kq_- - \operatorname{logit}(q_-)}{\kappa},
\qquad
\theta_+ = \frac{Kq_+ - \operatorname{logit}(q_+)}{\kappa}.
\]

There are exactly three fixed points if and only if

\[
\theta_-<\theta<\theta_+.
\]

Within this open interval, the low and high fixed points are locally stable and
the middle fixed point is unstable. At either endpoint, the system is at a
saddle-node and is not counted as strictly bistable.

The proof follows directly from

\[
F'(q)=\frac{1}{q(1-q)}-K
\]

and the fact that \(1/[q(1-q)]\ge4\) on \((0,1)\); a complete proof is provided
in Supplementary Mathematical Result S1.

This theorem is deliberately narrow. It is exact for the stated one-dimensional
sigmoid map. It is not a theorem for arbitrary positive-feedback ecological
systems, nor is it a theorem for the full finite trait–allele simulator.

### 3.3 Trait-mode lifting

Let \(m_H(q)\) denote the maximum viability margin of the designated high-trait
region under interaction state \(q\). Once the branch geometry and performance
surface are declared, the following is immediate.

**Conditional Result 1 (branch-dependent potential high-trait viability).** If
\(q_{\rm low}<q_{\rm high}\) are the stable interaction branches and

\[
m_H(q_{\rm low})<0<m_H(q_{\rm high}),
\]

then the high-trait component is potentially unavailable on the low branch and
potentially viable on the high branch.

This result concerns potential viability only. Realised high-trait abundance
requires a finite recruitment and inheritance closure and is not supplied by the
one-state theorem.

### 3.4 Patchwise non-additivity

The branch geometry also exposes why total habitat area alone is insufficient for
a fragmentation claim. When an ecological closure requires a patchwise condition
\(A_j>A_c\), redistributing a fixed total area among patches can move some patches
below the support region even when total area is unchanged. This is a Type C
statement because the threshold \(A_c\) must be derived or supplied by the
particular ecological closure.

---

## 4. Analytical result II: migration bounds on a patch network

### 4.1 Declared migration operator

For a network of patches, let allele frequencies update by

\[
p'_i=\sum_jM_{ij}p_j,
\]

where \(M_{ij}\) is the contribution of source patch \(j\) to destination patch
\(i\). We assume that \(M\) is destination-by-source row-stochastic:

\[
M_{ij}\ge0,
\qquad
\sum_jM_{ij}=1.
\]

This form permits asymmetric source–sink structure, stepping-stone corridors,
and distance-derived kernels. It describes composition after deterministic
mixing, not demographic rescue or long-run persistence.

### 4.2 Common-floor theorem

**Theorem 2 (common-floor preservation; Type T for the declared migration
operator).** If every source patch satisfies \(p_j\ge p_{\min}\), then every
destination satisfies \(p'_i\ge p_{\min}\).

Indeed,

\[
p'_i=\sum_jM_{ij}p_j
\ge\sum_jM_{ij}p_{\min}=p_{\min}.
\]

Thus deterministic migration cannot lower a common allele-frequency floor before
finite sampling. The statement does not imply that migration increases a focal
patch below the floor, that it offsets genetic drift, or that it produces
demographic or trait rescue.

### 4.3 Focal-rescue condition

If \(b_j\) is an independently established lower bound for source frequency in
patch \(j\), then the destination-specific bound is

\[
p'_i\ge\sum_jM_{ij}b_j.
\]

A target \(p_{\rm target}\) is therefore certified by deterministic mixing only
when

\[
\sum_jM_{ij}b_j\ge p_{\rm target}.
\]

This condition is sharp for the declared information and shows why nonzero
migration is not itself a rescue theorem: source composition and network weights
matter.

---

## 5. From analytical statements to a finite eco-genetic closure

The analytical results establish necessary structure, but not realised
population-level outcomes. We therefore use a finite-bin, two-kernel recruitment
closure in which trait bins are sampled stochastically and allele-linked
recruitment may feed back into interaction support. The finite experiment
separates a source-generation stage, full-state transfer, a conservation-aware
projection into an equal-isolated landscape, and subsequent stochastic
trajectories.

The relevant finite H1 and H3 results are Type S. A mutation-conditioned
high-state source could retain interaction memory under the declared continuation
closure. Conditional on such source preparation and conservation-preserving
projection, equal isolation lowered final interaction, local effective size, and
realised high-trait mass relative to a single large patch. These are results for
the stated finite trait-recruitment, symmetric allele-mutation, and
full-state-transfer closures; they are not extensions of Theorems 1 or 2.

The full-state transfer is important. Reconstructing only a low-dimensional
summary would confound fragmentation with a change of source state. By retaining
the prepared high-state composition before projection, the finite comparison
addresses the stated counterfactual: what changes when the same prepared source
is isolated into equal fragments?

---

## 6. Genetic warning as a censored first-passage problem

### 6.1 Event definitions

Let \(\tau_T\) be first post-baseline realised high-trait loss. For a diversity
summary \(H_x(t)\), where \(x\in\{\alpha,\gamma\}\), and a relative decline
fraction \(r\), define

\[
\tau_{\Delta H_x(r)}=
\inf\{t>0:H_x(t)\le(1-r)H_x(0)\}.
\]

The relative-warning question is

\[
\tau_{\Delta H_x(r)}<\tau_T,
\qquad r\in\{0.05,0.10,0.20\}.
\]

A crossing at \(t=0\) is not counted as an early warning. If either a warning or
trait-loss event is not observed within the finite horizon, the comparison is
censored rather than assigned the terminal generation.

### 6.2 Calibration before warning assessment

A deterioration schedule must not be selected because it produces favourable
warning order. The analysis therefore separated schedule calibration from
warning assessment. Calibration used post-baseline realised high-trait loss only
and selected a schedule only when every calibration seed block had a conditional
trait-loss frequency in the predeclared interval \([0.30,0.70]\). Diversity
values, warning times, leads, lags, and lead times were excluded from this stage.

The initial linear-ramp family selected no cell. A predeclared ramp-and-hold
family then selected one configuration:

```text
mutation rate = 0.10
A_ref = 0.8
kappa = 6.0
equal-isolated landscape
ramp = 30 generations
hold = 90 generations
normalized barrier increase = 0.15
```

The five calibration seed-block trait-loss frequencies were \(0.50,0.40,0.40,
0.50,0.50\). The other eleven primary cells were not carried forward to warning
validation.

### 6.3 Conditional and unconditional reporting

For every warning–trait comparison, the analysis retains the valid-pair count,
the conditional lead probability among observed pairs, and the all-replicate
observed-lead fraction. These quantities answer different questions. A
conditional lead probability of one is not interpretable without the number of
valid pairs and the amount of censoring. Conversely, a low unconditional fraction
can arise because few trajectories reached both events rather than because the
observed ordering was unfavourable.

---

## 7. Existing finite results

### 7.1 Relative warnings

Independent validation used fresh master seeds `20261110`–`20261114`, with 20
replicates per seed, in the locked configuration above. Of 100 attempted sources,
83 produced an H1-prepared and projection-supported trajectory. Thirty-five of
those 83 trajectories reached realised trait loss within the scheduled horizon;
the remaining 48 available trajectories are right-censored for trait loss.

For every one of the six predeclared endpoints—\(H_\alpha\) and \(H_\gamma\) at
relative declines \(r=0.05,0.10,0.20\)—all 35 valid same-replicate comparisons
had a diversity warning before realised trait loss. There were zero ties and zero
lags among valid pairs. This is Type S evidence that, in the selected finite
closure and among observed event pairs,

\[
\tau_{\Delta H_x(r)}<\tau_T,
\qquad x\in\{\alpha,\gamma\},\quad r\in\{0.05,0.10,0.20\}.
\]

The result does not establish an ordering in the censored trajectories, does not
generalize to the eleven calibration-unselected cells, and does not establish a
universal genetic early-warning law.

### 7.2 Absolute thresholds

The earlier absolute-threshold formulation used

\[
H_\alpha\le0.20,
\qquad H_\gamma\le0.20.
\]

These thresholds were retained unchanged for a secondary audit of the raw
independent-validation series. The audit found mixed ordering: among valid pairs,
\(H_\alpha\le0.20\) yielded 14 leads and 6 lags, whereas
\(H_\gamma\le0.20\) yielded 8 leads and 8 lags. Because observed lags prevent a
uniform ordering claim, neither fixed threshold is retained as a robust canonical
warning rule in the selected closure. This is a negative robustness result, not a
proof that any absolute threshold is impossible in every model or biological
system.

---

## 8. Discussion

The framework makes three distinctions that are often obscured in discussions of
eco-genetic criticality.

First, the existence of an interaction threshold is not equivalent to a result
about realised trait disappearance. The canonical sigmoid map supplies exact
branch geometry, and the trait-margin condition supplies potential viability,
but finite recruitment and inheritance are needed before realised trait occupancy
can be assessed.

Second, migration statements depend on what is being moved. The row-stochastic
frequency update provides exact composition bounds, but it does not describe
individual abundance, trait propagules, extinction–recolonisation, or demographic
rescue. Those require their own life-cycle map.

Third, an apparent genetic warning depends on its event definition and on its
observation process. The finite result retained here is not that a universal
fixed diversity threshold forecasts trait collapse. Rather, it is a conditional
relative-warning ordering observed after trait-loss-only schedule calibration in
one declared finite closure. The contrast between the relative and absolute
formulations is informative: the relevant signal in this closure was erosion
relative to a trajectory's own baseline, not a stable absolute diversity cut-off.

The main limitation is therefore also a guide for future work. The symmetric
mutation operator, finite trait-recruitment closure, full-state transfer rule,
and selected deterioration path define the boundary of the Type S result.
Changing any of them is an extension project, not a parameter tweak inside the
present claim. The present manuscript does not attempt such an extension.

---

## 9. Conclusion

A transparent eco-genetic criticality program can combine exact mathematics with
finite stochastic experiments without conflating their evidential status. For a
canonical positive-feedback map, branch geometry and stability are analytically
specified. For network mixing, common-floor and focal-rescue bounds are exact for
the declared matrix update. A finite trait–allele closure then permits explicit,
censoring-aware tests of realised trait and genetic event ordering. Within the
selected finite deterioration configuration, baseline-relative diversity erosion
preceded every observed realised trait-loss event, whereas fixed absolute
thresholds did not provide a robust rule. This distinction between theorem,
conditional result, and finite evidence is the central contribution of the
framework.

---

## Data and code availability

All equations, finite closures, tests, workflows, and evidence ledgers described
here are contained in the `eco-genetic-criticality` repository. The manuscript
draft introduces no new simulation output.

## Citation note

This first drafting pass intentionally contains no external literature citations.
A dedicated literature pass should fill the marked contextual citation slots
without changing the mathematical statements or finite-result ledger.
