# Eco-genetic criticality program: theorems, hypotheses, and dynamic simulations

## Central causal chain

```text
patch size
-> interaction intensity
-> trait-space topology
-> population size / effective reproductive size
-> genetic diversity.
```

The chain is a scientific hypothesis program, not itself a theorem. Each arrow
requires a declared state equation or life-cycle map.

```text
Type T  mathematical theorem under explicit assumptions
Type C  conditional theorem once an ecological closure is supplied
Type H  substantive dynamic hypothesis about an ecosystem model
Type S  simulation result for a declared model, not proof of T/C/H.
```

The current roadmap is:

```text
general theorem layer
-> canonical logistic corollary
-> potential viability
-> finite realised trait abundance and occupancy
-> stochastic genetic first-passage experiments
-> future mutation/recolonisation model
```

The standard/full phase-diagram profiles use a declared finite-bin,
two-kernel, coupled-feedback closure. They are Type S experiments and never
retroactively strengthen earlier theorem claims. A model-specific phase diagram
is not a theorem.

---

## Theorem layer

### G0 — finite transmission variance

For unbiased post-selection transmission with positive conditional variance:

```text
E[H(P') | p*] = H(p*) - 2 Var(P'|p*) < H(p*).
```

This is Type T. It makes no ecological claim about patch size or interaction.

### P0 — no-bistability certificate

For

```text
q_next = g{kappa(Aq-theta)}
```

with `M=sup|g'|`,

```text
kappa*A*M < 1
```

certifies one fixed point. Its converse does not prove bistability.

### P1 — trait-mode lifting

Given established stable states `q_L<q_H`, high-trait region `Z_H`, and

```text
m_H(q) = max_{z in Z_H}[W(z;q)-tau],
```

then

```text
m_H(q_L)<0<m_H(q_H)
```

implies branch-dependent potential high-trait viability. This is Type T once the
branches and performance map are declared.

### P2 — patchwise non-additivity

If an ecological mechanism requires `A_j>A_c` in each patch, total area alone
cannot guarantee it. This is Type C because the patchwise threshold must first
be derived or supplied by the ecological closure.

### G1 — conditional eco-genetic ordering

If `N_e=Psi(A,q,xi)` increases with q and transmission variance decreases with
N_e, low-q branches erode expected local diversity faster. This is Type C, not a
universal interaction-system statement.

---

## H_critical

> Positive frequency-dependent interaction can generate a patch-size transition
> across which a high-investment trait mode changes discontinuously.

A theorem for a specified system must establish:

1. a rigorous branch transition;
2. high-trait margin sign change across stable branches;
3. absence/presence of the potential trait component across them.

For the explicit one-state logistic reduction, the repository now supplies a
Type T canonical certificate: `K=kappa*A/A_ref>4` and an exact open barrier
interval give three fixed points (stable-low, unstable-middle, stable-high).
Combining that branch result with `m_H(q_L)<0<m_H(q_H)` certifies H1 for that
specified map. See `docs/canonical_h1_bifurcation.md`. This does not prove H1
for arbitrary positive-feedback systems or for the full coupled simulator.

The dynamic experiment sweeps patch size, barrier, feedback strength, and trait
cost/benefit parameters. It distinguishes no transition, smooth transition,
bistability under the declared model, and noise-sensitive transitions.

---

## H_genetic_lag

> Under identifiable conditions, genetic indicators can change before the
> realised high-trait mode disappears.

This is Type H. It is not implied by G0/G1 and is not universal.

The simulation predeclares:

```text
tau_trait_potential
tau_trait_realised
tau_allele_loss
tau_H_alpha
tau_H_gamma
tau_FST
```

and reports the valid-pair probability of inequalities such as

```text
tau_H_alpha < tau_trait_realised.
```

Events that do not occur are censored, not converted into terminal-generation
values. The warning threshold must be set before inspecting outcomes.

---

## H_fragmentation

> At fixed total area, subdivision can prevent maintenance of a high-investment
> trait mode and change associated genetic structure when more patches fall below
> the interaction-support region.

The ecological component is Type C: an interaction threshold alone only removes
that mechanism. High-trait loss additionally requires the declared trait margin
condition.

The genetic component is Type H. `H_alpha`, `H_gamma`, and `F_ST` are distinct:
fragmentation can lower local `H_alpha`, increase `F_ST`, and preserve or change
`H_gamma` differently.

Dynamic comparisons keep total area fixed across:

```text
one large patch
m equal isolated patches
m equal patches with controlled migration.
```

Each report retains q, N, N_e, potential viability, realised high-trait abundance,
allele frequency, H_alpha, H_gamma, F_ST, and event times separately.

---

## Simulation closure and stop rules

The finite-bin experiment layer declares:

```text
n_{j,t+1}(z_k) ~ Multinomial(N_{j,t+1}, pi_{j,t}(z_k))
```

with optional allele-linked two-kernel recruitment and optional trait/allele
feedback into q. This supports model-specific results only.

Report counterexamples rather than discarding them when any occur:

```text
positive feedback but no bistability in the tested range
bistability without high-trait margin sign change
potential viability loss without realised occupancy loss
realised trait collapse without a genetic lead
fragmentation lowers H_alpha but preserves H_gamma
interaction raises census N but lowers N_e through skew.
```

These are scientific results, not failed simulations.

---

## Current relation in one sentence

The theorem layer specifies what must be proved under explicit assumptions; the
finite occupancy phase-diagram layer tests which of those assumptions and event
orderings arise in the declared stochastic model.
