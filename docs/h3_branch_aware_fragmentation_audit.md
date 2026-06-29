# H3 branch-aware fragmentation and migration audit

H3 asks what equal subdivision and migration do at matched total area. Earlier
H3 comparisons did not distinguish whether a change occurred in a region where
finite H1 branch structure was present. This audit makes that condition
explicit.

## Shared branch starts

For each parameter cell, the **one-large** local canonical H1 certificate
provides a low and high stable interaction value. The finite simulator then
starts every patch in each landscape from those same two interaction values:

- `one_large`;
- `equal_isolated`; and
- `equal_migrating`.

Every low/high start and landscape variant shares the same parameter cell,
replicate index, seed, and initial trait/allele state. The only intended changes
are initial branch label and landscape structure.

Using one-large branches for all landscapes is deliberate. Equal fragments may
fall outside their own local canonical bistable region; treating their missing
local branch certificate as a different initial-condition definition would make
the spatial comparison circular.

## One-large finite H1 precondition

A replicate enters H3 branch-aware contrasts only if the one-large finite pair
shows both:

1. high-start minus low-start terminal mean interaction exceeds the declared
   separation threshold; and
2. high-start retains a potential high-trait component while low-start lacks it.

Thus H3 results are conditional on a finite H1 mechanism existing in the
matched unfragmented reference. Canonical-H1-inapplicable cells and finite
one-large precondition failures remain in the output but do not become evidence
for fragmentation failure, migration failure, or missing genetic warning.

## Outcomes reported within each branch

For low-start and high-start separately, every landscape records:

- terminal mean interaction;
- terminal mean **local** effective size;
- realised high-trait mass;
- potential high-trait viability;
- H-alpha, H-gamma, and FST;
- H-alpha, H-gamma, and allele-loss first-passage order relative to realised
  high-trait loss, including valid-pair counts and censoring; and
- H1 theorem-boundary scope.

The audit first asks whether the low/high branch distinction remains in each
landscape. It then reports, within a branch:

```text
equal_isolated - one_large
equal_migrating - equal_isolated
```

for interaction, local effective size, realised high-trait mass, H-alpha,
H-gamma, and FST where FST is defined.

## Interpretation of migration

The multipatch simulator represents migration as allele-frequency mixing. It
does not include the thresholded population establishment and source-to-
destination individual movement required to call an event rescue or
recolonisation. Consequently this audit uses only neutral language:

- migration increases or decreases a reported finite-model outcome;
- migration preserves or erodes a branch-conditioned quantity; or
- migration has no material effect under the declared closure.

Demographic rescue and recolonisation must still be tested with the separate H3
extinction--recolonisation lifecycle model.

## Reading results

The artifact has three nested denominators:

1. all stochastic replicates;
2. replicates where the one-large finite H1 mechanism persisted; and
3. warning/trait pairs where both first-passage events occurred.

Keep them separate. A low branch-aware warning rate can reflect failure of the
one-large H1 precondition, branch collapse after subdivision, no observed event
pair, a late warning, or a genuinely absent warning. Those are different
scientific outcomes.

This is Type S evidence for the declared finite closure. It does not prove that
fragmentation or connectivity has the same effect in every ecological system.
