# Finite coupled-chain audit

The canonical H1, H2, and H3 certificates are intentionally narrow. Their
composition requires explicit links from local support to interaction branch and
from interaction to effective population size. The finite coupled-chain audit
does not assume those links are true. It tests their joint occurrence in the
existing finite-bin multipatch simulator under matched landscapes.

## Question

For a fixed total area and one parameter cell, do equal isolated fragments,
relative to one large patch, jointly show:

1. lower mean final interaction;
2. lower mean **local** final effective size;
3. lower realised high-trait mass; and
4. an H-alpha warning before realised high-trait loss within the isolated
   replicate?

A replicate satisfies `finite_chain_supported` only when all four conditions
hold with the predeclared numerical tolerance. This is a stringent joint
simulation predicate, not a theorem certificate.

## Matched comparison

Every cell runs the existing three landscapes:

- `one_large`;
- `equal_isolated`; and
- `equal_migrating`.

All three variants receive the same parameter cell, replicate index, derived
seed, and initialisation rule. Differences therefore compare declared spatial
conditions rather than unrelated random-number streams.

The `equal_migrating - equal_isolated` contrasts are reported separately. A
positive migration contrast means the migrating model has a larger final value
for that metric; it is not automatically labelled rescue. Explicit rescue and
recolonisation remain the domain of the separate H3 lifecycle model.

## First-passage censoring

The joint predicate requires a valid isolated pair

```text
tau_H_alpha < tau_trait_realised.
```

If either event is absent within the declared horizon, the pair is censored and
the replicate cannot support the joint predicate. The audit reports:

- the valid-pair denominator;
- the censored-pair count;
- the conditional H-alpha lead probability among valid pairs; and
- the unconditional joint-chain support probability across all replicates.

These quantities must be read together. A small joint-support probability can
mean a failed mechanistic link, no warning lead, or simply insufficient observed
event pairs.

## Canonical H1 context and scope

For each scenario the artifact also records the canonical H1 certificate at
that scenario's local patch area. This provides analytic context for the local
area change, but it does not assert that the finite coupled simulator lies in
the canonical theorem limit. Per-replicate H1 theorem-boundary audits are
therefore stored and summarised separately.

## Interpretation

The audit advances validation from a composed set of assumptions to a single
finite stochastic model. It can support, qualify, or contradict the proposed
fragmentation-to-warning chain under the declared dynamics. It cannot establish
a universal ecological law or replace the individual canonical theorems.
