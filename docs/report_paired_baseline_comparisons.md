# Figure-ready paired-baseline comparison report

`paired_baseline_report.py` reads a flat paired-baseline CSV and creates annotated
heat maps. It does not rerun simulations; the source CSV remains the numerical
record and the report is a reproducible visualisation.

## Run

Create a paired comparison artifact first:

```bash
python scripts/run_paired_baseline_comparisons.py \
  --profile standard \
  --output-dir artifacts/baseline_comparison/standard \
  --master-seed 20260629
```

Then make figures from its CSV:

```bash
python scripts/report_paired_baseline_comparisons.py \
  artifacts/baseline_comparison/standard/baseline_comparison_standard.csv \
  --output-dir artifacts/baseline_comparison/standard/report
```

To make a smaller report, choose metrics explicitly with repeated `--metric`
arguments.

## Default figures

When these columns are present, the default report creates heat maps for:

- full model minus trait-only realised high-trait mass;
- full model minus genetic-only realised high-trait mass;
- full model minus trait-only final Hα;
- full model minus genetic-only final Hα; and
- canonical-H1 theorem-limit probability for each baseline.

Each figure keeps scenario and area-reference values separate. It never averages
across them.

## Reading contrasts beside theorem scope

A `full_minus_trait_only` or `full_minus_genetic_only` column is a matched,
within-simulator comparison: seed, initial state, landscape, and parameter cell
are held fixed. Positive values mean the full eco-genetic model was larger in
that comparison.

This does not automatically make the contrast a canonical-H1 theorem result.
Read the matching `scope.<baseline_id>.*` metrics beside it. A theorem-limit
probability of one means the relevant baseline remained in the single-patch,
no-migration canonical-H1 limit for every replicate in that cell. Lower values
mark controlled departures under the finite simulator.

`NA` is retained rather than converted to zero. For a conditional genetic-lead
metric, it usually means no uncensored warning-versus-trait-loss pair was
observed.

## Outputs

The output directory contains one PNG per selected metric × scenario ×
area-reference combination and a `REPORT.md` index with source metric plus
observed and missing parameter-cell counts.
