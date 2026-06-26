# Three-hypothesis proof-status ledger

This ledger distinguishes theorem statements, conditional theorem bridges,
model-specific simulation support, and claims that are false if read
universally.

## H1 — interaction criticality and high-trait potential viability

### Strong canonical statement

For the one-dimensional canonical map

```text
q_(t+1) = sigmoid[kappa(A q_t - theta)],
```

criticality, bistability, and the saddle-node threshold are analytically
characterised. If the high-trait viability margin has opposite signs on the
low- and high-interaction branches, high-trait potential viability is
branch-dependent.

**Status: proved for the declared canonical map.**

### Finite-bin simulator statement

The finite-bin multi-patch closure has a coupled state

```text
(q, N, p, trait abundance).
```

The canonical threshold does not transfer automatically. Phase diagrams show
transition-like regions, but those are simulation results.

**Status: model-specific simulation support; no global full-closure theorem.**

### Universal reading

```text
patch-size reduction always eliminates high traits
```

is false. Whether it occurs depends on feedback, trait margin, demographic
bounds, and connectivity.

## H2 — genetic warning before realised high-trait loss

### Universal reading

```text
genetic warning always occurs before realised trait loss
```

is false. The same recursion family admits lead, tie, trait-first, and censored
paths.

### Conditional results

L3 proves a deterministic sufficient condition:

```text
genetic decay upper bound + realised trait persistence lower bound
=> tau_H < tau_trait_realised.
```

L4 proves a finite-population high-probability version:

```text
E[H_t] <= H_0 lambda_bar^t
+ finite recruitment persistence bounds
=> lower bound on P(tau_H <= t < tau_trait_realised).
```

The finite-bin closure now derives `pi_min`, `n_min`, and the Wright--Fisher
sampling part of `lambda_bar` from a declared region. The remaining difficult
premises are a pre-sampling H-alpha expansion bound and a region-invariance
argument.

The stochastic refuge certificate supplies the latter for the no-migration,
single-refuge specialisation over a finite horizon.

**Status: conditional theorem for declared submodels; simulation witnesses for
the broader finite-bin landscape; not a global simulator theorem.**

## H3 — fragmentation non-additivity at fixed total area

### Mechanism statement

In the canonical interaction model, equal partition can put every subpatch
below the interaction capacity for a high-interaction branch even when a single
patch of the same total area can support it.

**Status: proved as a conditional capacity theorem in the canonical model.**

### Finite-bin simulator statement

One-large, equal-isolated, and equal-migrating scenarios produce distinct
realised trait and genetic outcomes at equal total area in the phase-diagram
pilot.

**Status: model-specific simulation support.**

### Universal reading

```text
fragmentation is always worse than one large patch
```

is false. Migration can rescue some outcomes; its effect is parameter- and
closure-dependent.

## Current scientific claim

The defensible central claim is not a universal ecological law. It is:

```text
Positive interaction feedback plus finite recruitment can generate parameter
regions in which fixed-total-area fragmentation changes trait persistence and
genetic first-passage ordering. Canonical and conditional theorems identify
sufficient mechanisms; finite-bin simulations map candidate regions.
```

## Related-work discipline

A supplied related paper should be compared at the level of its exact model,
state variables, and claim type before asserting novelty or contradiction. No
article-specific claim is recorded in this ledger unless the article text has
been checked.
