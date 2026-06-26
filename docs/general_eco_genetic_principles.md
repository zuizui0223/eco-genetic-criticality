# General mathematical principles for eco-genetic criticality

## Why this comes before an ABM

This document fixes the theorem layer before choosing a particular ecological
response curve, a Wright--Fisher life cycle, a selfing model, or an empirical
species. The goal is to distinguish:

```text
mathematical implication
from
special-model corollary
from
biological closure requiring validation.
```

The earlier RACH channel theorems remain the observation layer: they determine
what can be identified from a performance map. The present principles are a
dynamics layer: they determine what follows from finite transmission and
collective interaction feedback.

The relationship to the three central dynamic hypotheses
`H_critical`, `H_genetic_lag`, and `H_fragmentation` is specified in
`docs/eco_genetic_hypothesis_program.md`.

---

## 1. States and ordering

For a patch, write:

```text
A       patch size
q       interaction availability
z       individual trait
W(z;q) trait performance conditional on q
p       allele frequency after the adult selection step
P'      random next-generation allele frequency.
```

The temporal order is fixed:

```text
ecological state -> individual selection -> finite transmission -> genetic state.
```

The theory never treats genetic drift as a causal substitute for selection. It
acts after the individual fitness map has generated the post-selection frequency.

---

## Theorem G0 — finite transmission variance identity

Let

```text
H(p) = 2p(1-p)
```

be biallelic gene diversity. Let `p*` be post-selection frequency and `P'` the
random next-generation frequency. Assume only that the first two conditional
moments exist.

Then exactly:

```text
E[H(P') | p*]
= H(E[P' | p*]) - 2 Var(P' | p*).
```

### Proof

```text
H(P') = 2P' - 2(P')^2.
```

Taking expectations and using

```text
E[(P')^2] = Var(P') + E[P']^2
```

gives the identity. ∎

### Corollary G0.1 — unbiased finite transmission erodes expected diversity

If transmission is unbiased,

```text
E[P' | p*] = p*,
```

and has nonzero variance,

```text
Var(P' | p*) > 0,
```

then

```text
E[H(P') | p*]
= H(p*) - 2 Var(P' | p*)
< H(p*).
```

This is the general drift principle. It does **not** assume binomial
Wright--Fisher sampling. It applies to any finite transmission kernel satisfying
those moment conditions.

### What it does not say

It does not say `H_{t+1}<H_t` across the whole generation. Selection can change
`p` before transmission, including toward 1/2. The strict loss is relative to
heterozygosity **after selection and before finite transmission**.

Wright--Fisher is a special case:

```text
Var(P'|p*) = p*(1-p*)/(2N_e),
E[H(P')|p*] = [1-1/(2N_e)]H(p*).
```

---

## Theorem P0 — global feedback contraction bound

Let a patch interaction update be

```text
T_theta(q) = g{kappa(Aq-theta)},
```

where `g` is differentiable and has a global derivative bound

```text
sup_x |g'(x)| <= M < infinity.
```

Then the global Lipschitz constant of `T_theta` is at most

```text
L = kappa A M.
```

If

```text
kappa A M < 1,
```

then `T_theta` is a contraction. Therefore it has exactly one fixed point for
every barrier `theta`.

### Proof

By the chain rule,

```text
|T'_theta(q)| <= kappa A M.
```

The contraction mapping theorem gives unique fixed point and global convergence
under iteration. ∎

### Critical interpretation

```text
kappa A M < 1
```

is a **universal no-bistability certificate** for this response class.

```text
kappa A M >= 1
```

does not prove bistability. It only means this universal uniqueness proof no
longer applies. Additional response-shape conditions are required to establish
saddle nodes.

For logistic `g`, `M=1/4`, so the familiar value

```text
A_c = 4/kappa
```

is a special-model corollary, not the general theorem.

---

## Theorem P1 — trait-mode lifting from ecological branches

Let a high-trait region be `Z_H` and define its viability margin at interaction
state `q` by

```text
m_H(q) = max_{z in Z_H} [W(z;q)-tau].
```

Suppose an ecological system has two stable interaction branches

```text
q_L < q_H.
```

If

```text
m_H(q_L) < 0 < m_H(q_H),
```

then

```text
Omega_tau(q_L) ∩ Z_H = empty set,
Omega_tau(q_H) ∩ Z_H is not empty.
```

Thus the high-trait mode is branch- and therefore history-dependent.

### Proof

The sign of `m_H` is exactly the statement of whether any trait in `Z_H` clears
threshold `tau`. Apply it separately at the two stable ecological states. ∎

This theorem does not assume an externally invented `q_required`; the threshold
is derived from the declared trait-performance function.

---

## Theorem P2 — patchwise threshold implies habitat non-additivity

Assume a collective interaction mechanism is possible only within a patch when

```text
A_j > A_c.
```

For a landscape with patch sizes `A_1,...,A_m`, total area alone is insufficient:

```text
sum_j A_j > A_c
```

does not imply the mechanism is possible. If

```text
max_j A_j <= A_c,
```

no patch can express the mechanism, even if total area exceeds `A_c`.

This is a conditional theorem: the biological model must first justify a
patchwise threshold. Once that condition is established, area is non-additive for
that mechanism.

---

## Theorem G1 — conditional eco-genetic branch ordering

Let effective reproductive size be a life-cycle-derived function

```text
N_e = Psi(A,q,xi),
```

where `xi` includes selfing, reproductive skew, sex ratio, seed bank, clonality,
migration, and other demography.

Suppose an ecological system has stable states `q_L<q_H`, and in the relevant
parameter region:

```text
Psi(A,q_L,xi) < Psi(A,q_H,xi),
Var(P'|p*,N_e) decreases as N_e increases.
```

Then

```text
Var(P'|p*,q_L) > Var(P'|p*,q_H),
```

and by G0, for unbiased transmission,

```text
E[H(P')|p*,q_L] < E[H(P')|p*,q_H].
```

Thus ecological branch switching changes expected genetic-erosion rate.

### Scope

The monotonic relation between interaction availability and effective size is a
biological closure, not a universal theorem. In some systems interaction can
increase reproductive skew and reduce `N_e`; such a system does not satisfy this
hypothesis and cannot invoke G1 in this direction.

---

## Species-level diversity is deferred carefully

For patches with frequencies `p_j` and weights `w_j`:

```text
H_alpha = sum_j w_j 2p_j(1-p_j)
p_bar   = sum_j w_j p_j
H_gamma = 2p_bar(1-p_bar)
F_ST    = 1-H_alpha/H_gamma, when H_gamma>0.
```

G0 and G1 directly govern finite within-patch transmission and therefore the
local `H_alpha` layer. General theorems for `H_gamma` and `F_ST` require explicit
migration, mutation, extinction-recolonisation, and between-patch selection.
They are not claimed yet.

---

## Relation to special models and repository branches

- `patch_interaction_bifurcation_theory` is a logistic special case and may be
  retained later as a worked corollary after P0 is adopted.
- `patch_genetic_drift_theory` is a Wright--Fisher and chosen-N_e-closure special
  case of G0 and G1.
- Neither special model should be presented as the general theorem.
- ABMs come after this layer to test which assumptions survive finite populations,
  spatial structure, selfing, reproductive skew, and migration.