# Censoring-aware first-passage reporting for H2

A claim that a genetic warning precedes realised trait loss compares two
first-passage times,

\[
\tau_{\rm warning} < \tau_{\rm trait}.
\]

In a finite simulation horizon either event may fail to occur.  Such runs are
**censored**: they are not assigned the terminal generation as an artificial
event time.

The repository now reports three quantities for every warning-versus-trait
comparison:

1. **Valid-pair count**: the number of replicates where both event times are
   observed.
2. **Conditional lead probability**:
   \[
   P(\tau_{\rm warning}<\tau_{\rm trait}\mid\text{both events observed}).
   \]
3. **Unconditional observed-lead fraction**:
   \[
   P(\text{an observed warning lead}),
   \]
   whose denominator is every replicate and therefore retains censoring in the
   summary.

These quantities answer different questions.  A conditional lead probability
of 1 based on one valid pair is weak evidence if most replicates are censored;
an unconditional observed-lead fraction of zero can mean either no observed
leads or no valid event pairs.  The valid-pair count and observability must
therefore always accompany either probability.

`causal_model.first_passage_reporting.compare_first_passage_times` is the
shared implementation.  `EnsembleSummary` exposes the alpha-warning versus
trait-absence comparison and retains the older unqualified
`genetic_lead_probability` only as an alias for the all-replicate observed-lead
fraction.  New analysis should use the explicitly named conditional and
unconditional fields.

## Phase-diagram artifacts

`censoring_aware_phase_diagram_rows(results)` adds a parallel artifact format
for `CellResult` objects from the finite-bin parameter grid.  Each parameter /
scenario cell contains a separate comparison of realised trait loss against
`tau_H_alpha`, `tau_H_gamma`, `tau_FST`, and `tau_allele_loss`.  The output is
flat enough for CSV while retaining the following fields for every comparison:

```text
replicate_count
valid_pair_count
censored_pair_count
valid_pair_probability
lead_count
conditional_lead_probability
unconditional_observed_lead_fraction
time_differences
median_time_difference
```

The original `CellResult.to_csv_row()` remains unchanged for compatibility.
Use the censoring-aware rows for any H2 phase-diagram interpretation or figure
that compares warning order across parameter cells.
