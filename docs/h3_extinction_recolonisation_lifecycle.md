# H3 finite extinction–recolonisation lifecycle

## Why this closure exists

The earlier H3 migration theorem concerns allele frequencies only.  A convex
combination can preserve a common allele floor, but that fact alone cannot show
that a patch is demographically rescued, that a realised high-trait population
returns, or that an empty patch is recolonised.

`causal_model.network_h3_lifecycle` supplies a separate, declared stochastic
closure for those questions.  It is a **Type S** finite-population experiment,
not a general theorem about fragmented ecosystems.

## Patch state

At each census, each patch stores

\[
(N_i, T_i, C_i),
\]

where \(N_i\) is census population size, \(T_i\) is realised high-trait
abundance, and \(C_i\) is the number of high-allele copies among \(2N_i\)
diploid copies.  Empty patches have state \((0,0,0)\).

The model reports population, realised high-trait abundance, allele frequency,
\(H_\alpha\), \(H_\gamma\), and \(F_{ST}\) separately.  Thus a high trait,
an allele, and gene diversity cannot be silently treated as the same state.

## Life cycle

For each generation the model performs:

1. Adult survival, separately for high- and low-trait individuals.
2. Individual emigration with probability \(d\).
3. Directed dispersal through the source-to-destination kernel
   \(K_{ij}=P(\text{destination }j\mid\text{source }i)\).
4. Persistence, extinction, rescue, or recolonisation using declared abundance
   thresholds.
5. Density-limited recruitment up to the patch capacity, with an optional
   high-trait recruitment multiplier.

Allele copies among resident adults, migrants, and recruits are binomially
sampled conditional on their parental patch frequency.  This is an explicit
finite sampling closure, not a claim of exact Mendelian linkage between the
allele and trait bins.

## Status events

A patch can be labelled:

- `persisted`: it clears the persistence threshold without a threshold-crossing
  rescue;
- `rescued`: its resident adults lie below the persistence threshold but arrival
  lifts candidate abundance above it;
- `extinct`: an occupied patch fails the persistence threshold;
- `recolonised`: a previously empty patch receives at least the colonisation
  threshold and clears local persistence;
- `empty`: an empty patch receives too few migrants to found a colony.

The distinction between `rescued` and `recolonised` is important: a migration
path can maintain an existing small population without being capable of
founding a new colony, and vice versa.

## H3 questions now testable

At matched total capacity, the closure can compare one large patch, isolated
partitions, and networks with controlled kernels.  It can test whether
fragmentation changes:

- the probability of metapopulation extinction;
- realised high-trait loss and reappearance;
- local versus global diversity; and
- the frequency of demographic rescue versus colonisation.

Those comparisons should report the landscape kernel, capacities, thresholds,
replicate count, and censoring of non-occurring events.  A result under one
parameter set is model-specific evidence, not proof that connectivity is
always beneficial: the same framework can show migration rescue, migration
erosion, or no material effect.

## Relation to the H3 network theorem

The migration-matrix theorem uses a **destination-by-source row-stochastic**
frequency update.  This lifecycle uses a **source-to-destination row-stochastic**
transport kernel for individual emigrants.  The different conventions are
intentional: one tracks frequency mixing at a destination; the other allocates
finite individuals leaving a source.  Do not interchange the matrices without
transposing and rechecking their stochasticity conditions.
