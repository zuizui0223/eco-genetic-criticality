# Finite H1--H3 validation campaign

The project now has several complementary finite-model audits. Running them
with independently chosen grids would make it unclear whether their apparent
agreement or disagreement reflects biology, stochasticity, or merely different
parameter choices. The finite validation campaign runs the existing audits under
one immutable `ExperimentSpec`, one set of three matched landscapes, and one
master-seed rule.

It is a **reproducibility and falsification layer**, not a new theorem and not a
new ecological model.

## Included audits

The campaign writes separate raw CSV/JSON artifacts for all of these checks:

1. **Finite H1 branch separation**: same-seed low-start/high-start runs at a
   fixed barrier.
2. **Finite H1 hysteresis continuation**: stateful rising/falling barrier
   paths; terminal population, interaction, allele state, and realised trait
   bins carry into the next stage.
3. **Finite H2 branch-conditioned warning**: H-alpha, H-gamma, and allele-loss
   first passage relative to realised high-trait loss, only after the finite H1
   branch-mechanism precondition is met.
4. **Finite coupled chain**: matched one-large/equal-isolated/equal-migrating
   tests of the proposed fragmentation -> interaction -> local effective size
   -> realised trait -> H-alpha warning chain.
5. **H3 branch-aware fragmentation and migration**: how spatial structure
   changes branch retention, trait, diversity, and warning outcomes conditional
   on a finite H1 mechanism in the one-large reference.

The `ledger.csv` and `ledger.json` align these audit summaries by parameter
cell. They are for navigating agreement, scope limits, censoring, and
counterexamples. Each raw artifact remains the authoritative place for its
replicate-level interpretation.

## Running

A small smoke run is appropriate for checking installation and artifact shape:

```bash
python scripts/run_finite_validation_campaign.py \
  --profile quick \
  --output-dir artifacts/finite_validation_campaign \
  --prefix smoke \
  --replicates 2 \
  --generations 8 \
  --hysteresis-stage-generations 3
```

A moderate exploratory run can use `--profile standard`. The `full` profile is
opt-in and should be run outside normal CI because it evaluates a larger grid
and more replicates. Override `--replicates`, `--generations`, and
`--master-seed` only as explicit declared design choices; every override is
stored in the manifest.

## Manifest and no-selection rule

The command writes a manifest containing:

- profile and complete `ExperimentSpec`;
- all audit thresholds and continuation settings;
- the three landscape IDs;
- optional code revision;
- the output paths; and
- the policy that every declared cell, inapplicable cell, finite-precondition
  failure, censored event pair, and counterexample is retained.

Do not replace a censored first-passage event with the terminal generation. Do
not report a conditional warning lead probability without its H1-preconditioned
replicate count and valid-pair denominator. Do not treat allele-frequency
mixing in the multipatch model as demographic rescue or recolonisation.

## Reading the campaign

A parameter cell can legitimately show any combination of outcomes. For
example, canonical H1 may be applicable but finite branch separation may fail;
finite branch separation may persist but no H2 event pair may be observed; or
migration may raise H-alpha while eroding FST. The campaign is designed to keep
these distinctions visible rather than reduce them to one global success score.

The publishable claim remains calibrated to the specified system: analytic
canonical statements are Type T; explicit closure-dependent statements are Type
C; and all finite campaign patterns are Type S.
