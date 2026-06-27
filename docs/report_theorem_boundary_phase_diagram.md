# Reporting theorem-boundary phase diagrams

The report generator reads an existing flat CSV artifact.  It never reruns the
simulation or changes its aggregation.  This distinction matters: the CSV and
manifest are the numerical record; generated PNG files are visual summaries of
that record.

## Command

```bash
python scripts/report_theorem_boundary_phase_diagram.py \
  artifacts/theorem_boundary/standard/theorem_boundary_standard.csv \
  --output-dir artifacts/theorem_boundary/standard/report
```

By default, the report produces individual heat maps for:

1. realised high-trait persistence probability;
2. conditional Hα genetic-lead probability;
3. mean maximum canonical-update residual; and
4. canonical H1 theorem-limit probability.

Each figure is indexed by landscape scenario and area-reference value.  The
horizontal axis is interaction barrier; the vertical axis is interaction
feedback.

## Interpreting `NA`

The conditional genetic-lead metric is defined only among replicates in which
both the Hα warning and realised trait-loss events are observed before the
censoring horizon.  A cell labelled `NA` therefore means that no valid event
pair was available.  It is not a probability of zero and must not be plotted or
reported as one.

## Figure logic

Read outcome maps beside scope maps:

- **High persistence with zero residual** is an outcome in the declared
  canonical H1 limit.
- **High persistence with nonzero residual** is a robustness result after a
  named departure (density, trait feedback, allele feedback, migration, or
  multiple patches).
- **Low persistence with nonzero residual** does not refute H1; it identifies
  a region that requires stochastic-model interpretation and eventually
  empirical evaluation.

`REPORT.md` lists every generated figure and the number of observed versus
missing parameter cells.  Keep it beside the source manifest when preparing
figures for the manuscript.
