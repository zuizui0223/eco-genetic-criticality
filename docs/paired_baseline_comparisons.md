# Paired trait, genetic, and full eco-genetic comparisons

This runner evaluates three **within-simulator causal ablations** under identical
landscapes, initial state, parameter cell, replicate index, and derived random
seed.  It is not a claim that the ablations reproduce every external
trait-only, population-genetic, or metapopulation model.

## The three variants

The full model has declared interaction support

\[
S_t = \alpha q_t + \beta x_{H,t} + \gamma p_t,
\]

where \(q_t\) is interaction, \(x_{H,t}\) is realised high-trait mass, and
\(p_t\) is high-allele frequency.  The simulator also has the declared
allele-to-demography contribution `high_allele_growth * p_t`.

| Variant | Interaction support | Trait recruitment | Allele demographic contribution |
|---|---|---|---|
| `trait_only` | \((\alpha+\gamma)q_t + \beta x_{H,t}\) | resident trait only | set to 0 |
| `genetic_only` | \((\alpha+\beta)q_t + \gamma p_t\) | resident trait only | declared full setting |
| `full_eco_genetic` | \(\alpha q_t + \beta x_{H,t} + \gamma p_t\) | declared full setting | declared full setting |

The removed feedback weight is reassigned to current interaction \(q_t\), so
all three variants retain the same total interaction-support weight.  This is
important: a contrast should not merely say that one model received less total
support.

In both ablations the allele state is still simulated.  For `trait_only` it is
an observed drifting/selected state with no causal path into interaction, trait
recruitment, or population growth; this also blocks the indirect
\(p \rightarrow N \rightarrow\) density \(\rightarrow q\) path.  For
`genetic_only`, trait occupancy remains an ecological response to interaction,
but it cannot support interaction or receive allele-linked recruitment.

## Run a comparison

```bash
python scripts/run_paired_baseline_comparisons.py \
  --profile standard \
  --output-dir artifacts/baseline_comparison/standard \
  --master-seed 20260627
```

For a targeted pilot:

```bash
python scripts/run_paired_baseline_comparisons.py \
  --profile quick \
  --scenario equal_migrating \
  --replicates 10 \
  --generations 20 \
  --master-seed 20260627 \
  --output-dir artifacts/baseline_comparison/pilot
```

The runner writes:

- a flat CSV with per-cell model summaries and paired contrasts;
- a full JSON with every matched replicate outcome; and
- a manifest with the baseline equations, profile, scenario set, seed rule, and
  output names.

## Reading the paired contrasts

For each focal ablation, the CSV contains `full_minus_trait_only` or
`full_minus_genetic_only` contrasts.  They include mean paired differences in:

- realised high-trait persistence;
- realised high-trait mass;
- final Hα; and
- final high-allele frequency.

The `full_trait_mass_greater_probability`, tie probability, and less probability
are descriptive fractions across matched replicates.  They are not p-values.

The Hα genetic-lead field remains conditional on uncensored event pairs.  A
missing conditional probability means no valid pair, not evidence that genetic
lead has probability zero.

## Scope

This comparison establishes what the declared coupled simulator does under
controlled channel removals.  It does not prove that coupling has the same
magnitude in any biological system.  Empirical support still requires the H1–H3
measurement contract and an appropriate field or genomic dataset.
