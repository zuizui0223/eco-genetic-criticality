# Empirical measurement contract for H1–H3

The repository can establish specified-system theorems and model-specific
simulation results without biological field data.  It cannot estimate an
empirical interaction threshold, warning lead, or connectivity effect unless
those measurements are actually available.

`audit_empirical_columns(columns)` makes this distinction explicit before any
parameter fitting.

## H1: interaction criticality

Required columns:

```text
patch_id
patch_area
interaction_state
trait_value
performance
```

`interaction_state` is the measured feedback variable or a predeclared proxy.
`performance` must be a fitness-relevant response, such as survival, seed set,
recruitment, or a declared viability component.  A trait observation alone does
not identify a branch transition.

## H2: genetic warning before realised trait loss

Required columns:

```text
patch_id
time
realised_high_trait_abundance
sample_size
high_allele_copies
```

Repeated time points are essential.  `sample_size` and `high_allele_copies` are
required instead of a frequency alone, so finite-sampling uncertainty is
recoverable.  Warning thresholds and event definitions must be fixed before
examining event order.

## H3: fragmentation, rescue, and erosion

Required columns:

```text
patch_id
time
census_population
realised_high_trait_abundance
high_allele_copies
sample_size
source_patch_id
destination_patch_id
dispersal_count
```

The last three columns describe observed, estimated, or independently derived
movement.  A geographic distance matrix can be used only after a declared model
relates distance to dispersal counts or kernel weights.  Connectivity must not
be inferred from allele frequencies alone when claiming demographic rescue.

## What readiness means

A `ready=True` audit means the table has the minimum named columns required to
attempt the corresponding empirical analysis.  It does not certify adequate
sample size, measurement quality, causal identification, or model fit.  Those
remain study-specific design and diagnostic questions.
