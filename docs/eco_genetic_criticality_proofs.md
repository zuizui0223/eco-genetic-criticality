# Formal proof program for eco-genetic criticality

## Scope

This document upgrades the three central hypotheses into a sequence of precise
mathematical statements. It distinguishes:

- results proved for a **canonical logistic interaction model**;
- results valid only under explicitly named conditions;
- statements that remain dynamic hypotheses to be tested in a multi-patch model.

The target causal sequence is

```text
patch size A
-> interaction state q
-> trait-space topology Omega_tau
-> effective reproductive size N_e
-> genetic diversity.
```

No empirical species is assumed here.

---

# I. Canonical model for H_critical

Let the patch interaction state satisfy

```text
q_{t+1} = sigma[kappa(Aq_t-theta)],
```

where

```text
sigma(x)=1/(1+exp(-x)),
A>0,
kappa>0,
theta is an external interaction barrier.
```

Let trait performance be a continuous function

```text
W: Z x [0,1] -> R,
```

and define viable trait space at threshold `tau` by

```text
Omega_tau(q) = {z in Z : W(z;q) >= tau}.
```

Let `Z_H` be a declared high-investment trait region and define its viability
margin

```text
m_H(q) = max_{z in Z_H}[W(z;q)-tau].
```

The trait mode exists exactly when `m_H(q)>=0`.

## Theorem C1 — exact interaction saddle-node threshold

For the canonical logistic map:

```text
A <= 4/kappa
```

implies exactly one fixed point for every `theta`, whereas

```text
A > 4/kappa
```

implies a nonempty interval of `theta` with three fixed points.

### Proof

At fixed points the slope is

```text
kappa A q(1-q).
```

The equation for tangency is

```text
kappa A q(1-q)=1.
```

There are two distinct interior solutions iff `kappa A>4`:

```text
q_-=[1-sqrt(1-4/(kappa A))]/2,
q_+=[1+sqrt(1-4/(kappa A))]/2.
```

Substitution into

```text
theta(q)=Aq-logit(q)/kappa
```

produces two ordered saddle-node barriers. The outer fixed points are stable and
the middle fixed point is unstable because the slope is below one outside
`(q_-,q_+)` and above one inside. ∎

## Theorem C2 — canonical proof of H_critical

Assume `A>4/kappa`, choose a barrier inside the three-fixed-point interval, and
let the stable states be `q_L<q_H`. If

```text
m_H(q_L)<0<m_H(q_H),
```

then

```text
Omega_tau(q_L) intersect Z_H = empty set,
Omega_tau(q_H) intersect Z_H is not empty.
```

Therefore the high-investment trait mode is history dependent and switches
discontinuously when the selected stable interaction branch is destroyed.

### Proof

The sign of `m_H(q)` is equivalent to the existence of a viable trait in `Z_H`.
C1 gives the two stable interaction branches and their discontinuous collapse /
recovery transitions. Apply the margin condition to the two branches. ∎

### Meaning

This is a theorem for the canonical model, not a theorem that every positive
frequency-dependent interaction has a critical patch size. The general theorem
P0 in `general_eco_genetic_principles.md` provides only a no-bistability
certificate; C1 uses the logistic shape to go further.

---

# II. Canonical fragmentation corollary

## Theorem C3 — fixed-total-area fragmentation can remove the high trait mode

Let total area `A_total` be split into `m` isolated equal patches, so each has

```text
A_patch=A_total/m.
```

Assume the same canonical interaction map and the same margin sign condition as
C2 whenever a high branch exists. If

```text
A_total > 4/kappa
and
A_total/m <= 4/kappa,
```

then the one-patch landscape is capable of the C1 bistable interaction mechanism,
whereas no member of the equal-patch landscape is capable of it.

Consequently, the equal-patch landscape cannot maintain the high-trait mode by
that mechanism.

### Proof

C1 requires a patch area strictly greater than `4/kappa` for the three-fixed-point
interval. The single patch clears this threshold; every equal subpatch does not.
Without a high interaction branch, C2's high-trait branch state cannot exist in
any subpatch. ∎

### Scope

C3 concerns isolated equal patches and this mechanism. It does not say all
fragmentation causes trait loss, nor does it yet claim a species-wide genetic
diversity outcome.

---

# III. H_genetic_lag is not universal

Let `h_t` be expected local gene diversity in a declared patch or weighted local
alpha-diversity quantity. Suppose its life-cycle recursion can be written

```text
h_{t+1}=lambda_t h_t,
lambda_t>0.
```

The multiplier includes selection and finite transmission. It is not assumed
constant, neutral, or Wright--Fisher.

Let

```text
T = tau_trait
```

be the first time that the high-trait mode is absent and let

```text
tau_H = inf{t>=0 : h_t <= h_warn}
```

for a warning threshold chosen before inspecting the trajectory.

## Proposition L1 — exact first-passage criterion

A genetic lead exists iff

```text
exists t<T such that
product_{s=0}^{t-1} lambda_s <= h_warn/h_0.
```

### Proof

Repeated substitution gives

```text
h_t=h_0 product_{s=0}^{t-1}lambda_s.
```

Insert this into the definition of `tau_H<T`. ∎

## Corollary L1.1 — genetic lag is not a universal theorem

For the same initial diversity, warning threshold, and trait-collapse time,
there exist positive multiplier sequences with a genetic lead and sequences with
no genetic lead.

### Construction

Choose all `lambda_t` sufficiently close to one before `T`; the product stays
above `h_warn/h_0`, so no lead occurs. Choose early multipliers sufficiently
small; their product crosses the threshold before `T`, so a lead occurs. ∎

Thus H_genetic_lag must be a **conditional theorem** or model-specific dynamic
hypothesis, never an unconditional consequence of finite drift.

## Theorem L2 — uniform sufficient condition for a lead

If, before trait collapse,

```text
lambda_t <= lambda_bar < 1
```

for every `t<T`, then a lead is guaranteed when

```text
ceil[log(h_warn/h_0)/log(lambda_bar)] < T.
```

This is sufficient, not necessary.

## Theorem L3 - decay plus persistence sufficient condition

L2 assumes that the realised trait-collapse time `T` is already known. To turn
H_genetic_lag into a model theorem, `T` must itself be bounded from the trait
dynamics.

Assume two independently proved bounds over the same time axis:

```text
h_t <= h_0 lambda_bar^t,        0 < lambda_bar < 1,
N_H,t >= L_t > 0                for every t<T.
```

The first is a genetic-erosion upper bound. The second is a realised high-trait
persistence lower bound. The second condition implies

```text
tau_trait_realised >= T.
```

Let

```text
t_warn^* = ceil[log(h_warn/h_0)/log(lambda_bar)].
```

If

```text
t_warn^* < T,
```

then

```text
tau_H <= t_warn^* < T <= tau_trait_realised,
```

and therefore

```text
tau_H < tau_trait_realised.
```

### Meaning

This is the missing bridge between the phase-boundary witness and a theorem. A
simulation cell with many leads suggests where such a bridge might hold, but it
does not supply the bridge. A proof must separately establish:

```text
genetic decay bound
trait persistence bound
```

for the same closure and parameter region.

The implementation in `causal_model.eco_genetic_lag_theory` exposes this as a
certificate, not as an inference engine. It checks the logical implication once
the two bounds have been supplied.

---

# IV. What remains for the full three-hypothesis dynamic model

The next model must explicitly specify

```text
q_{j,t+1}=G(A_j,q_{j,t},mu_{j,t},N_{j,t};theta),
mu_{j,t+1}=S(mu_{j,t};W(z;q_{j,t})),
N_{e,j,t}=Psi(A_j,q_{j,t},mu_{j,t},xi),
P_{j,t+1}~K(p^*_{j,t},N_{e,j,t},xi),
```

plus migration and extinction/recolonisation.

Only then can it ask whether C2/C3 survive finite population stochasticity and
whether the multiplier sequence satisfies L1 or L2. Simulation is used to map
parameter regions and identify counterexamples; it never upgrades a simulated
pattern into a general proof.
