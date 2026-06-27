# H3 fixed-total-capacity phase-diagram protocol

`h3_phase_diagram.py` converts the finite extinction–recolonisation lifecycle
into a reproducible landscape comparison.  Every cell uses a stated total
capacity, initial global population, initial high-trait abundance, and initial
allele-copy total.

## Matched landscapes

The standard comparison always contains:

```text
one_large
 equal_isolated
 equal_complete_network
```

`one_large` holds the entire capacity in one patch.  The other two divide that
same capacity evenly among the declared number of patches.  The connected case
uses symmetric source-to-destination transport after emigration; the
emigration probability remains a separate life-cycle parameter.

## Exact initial composition

The runner first calculates global totals, then apportions those totals across
patches proportional to capacity.  It does this successively for census
population, high-trait individuals, and diploid high-allele copies.  Therefore
partitioning cannot change the global starting composition through independent
rounding in each patch.

## Required outputs

Every row includes the scenario, total capacity, patch count, survival,
emigration, thresholds, generations, and replicate denominator, plus:

- metapopulation extinction probability;
- realised high-trait-loss probability;
- recolonisation probability;
- rescue probability;
- medians of final occupied patches, high-trait patches, H_alpha, H_gamma, and
  F_ST.

Write both CSV and JSON with `write_h3_phase_diagram_artifacts`.  The JSON is
for a self-describing analysis record; the flat CSV is for plotting and model
comparison.

## Interpretation rules

Do not interpret a lower extinction probability in the connected landscape as
a general migration benefit.  The result is conditional on the declared
survival, recruitment, colonisation threshold, initial state, and transport
kernel.  Report all parameter cells, including cells where connectivity causes
no detectable change or lowers diversity through homogenisation.

A recommended first grid is survival × emigration with at least 100 replicates,
then a sensitivity grid over colonisation threshold and high-trait recruitment
multiplier.  The resulting support and counterexample regions are the H3 phase
diagram; they are model-specific evidence, not a theorem or empirical estimate.
