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
