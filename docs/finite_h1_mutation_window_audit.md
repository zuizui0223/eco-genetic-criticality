# Symmetric-mutation H1 polymorphism-window audit

## Why a new closure is needed

The zero-mutation finite closure supports H1 route memory, but the held high
branch reaches high-allele fixation:

\[
p=1,\qquad H_\alpha=H_\gamma=0.
\]

That state cannot test whether diversity decline warns before trait loss, because
diversity is already absent. It also cannot test a fragmentation-driven genetic
contrast. This audit introduces a separate, explicit genetic closure rather than
relabeling fixation as an early warning.

## Mutation closure

Each generation uses the existing life-cycle order through selection and
migration. Immediately before the existing finite-drift draw, allele frequency
is transformed by symmetric allele-state mutation:

\[
p_{\rm mut}=\mu+(1-2\mu)p.
\]

Here \(\mu\) is a per-generation probability for both high-to-low and
low-to-high allele-state transitions. This is not an empirically calibrated
biological mutation rate. It is a numerical persistence mechanism whose
consequences must be reported as Type S.

The legacy simulator is untouched. At \(\mu=0\), the mutation runner delegates
directly to the legacy implementation, preserving the same stochastic trajectory.

## What is screened

For every declared mutation rate, H1 parameter pair, master seed, and replicate,
the audit repeats the existing full-state H1 route protocol:

```text
nested 25 -> 49 -> 97 barrier calibration
-> common interior grid anchor
-> rising high and falling low route replay
-> fresh 30-generation full-state hold
```

It then computes the high-hold allele frequency and, for the one-large high
state, \(H_\alpha=H_\gamma=2p(1-p)\).

A replicate is `screen_supported` only when both conditions hold:

```text
finite H1 full-state hold supported
and
polymorphism_epsilon < p < 1 - polymorphism_epsilon
and
H-alpha > h_alpha_warning_threshold
and
H-gamma > h_gamma_warning_threshold
```

No H2 time ordering and no H3 fragmentation effect is evaluated in this screen.

## Screen versus validation

The workflow defaults to five master seeds and five replicates per parameter pair
for a coarse rate screen. A rate with nonzero screen support is only a candidate.
It must be revalidated with an independent master-seed ensemble and the standard
20 replicates before it can seed a later H2/H3 chain.

Rates without H1 support, rates with H1 support but fixed high branches, and
rates with polymorphism but no H1 hold all remain in the artifact. The analysis
does not delete unfavourable rates.

## Initial manual run

```bash
python scripts/run_finite_h1_mutation_window_audit.py \
  --profile standard \
  --replicates 5 \
  --master-seed 20260630 \
  --master-seed 20260631 \
  --master-seed 20260632 \
  --master-seed 20260633 \
  --master-seed 20260634 \
  --mutation-rate 0.0 \
  --mutation-rate 0.05 \
  --mutation-rate 0.10 \
  --mutation-rate 0.15 \
  --mutation-rate 0.20 \
  --endpoint-padding-fraction 0.5 \
  --stage-generations 30 \
  --hold-generations 30 \
  --barrier-points 25 \
  --barrier-points 49 \
  --barrier-points 97 \
  --interaction-separation-threshold 0.05 \
  --maximum-normalized-bracket-width 0.03 \
  --polymorphism-epsilon 1e-12 \
  --output-dir artifacts/h1_mutation_window \
  --prefix standard_mutation_screen_v1
```

The principal pair-level columns are:

```text
summary.denominators.h1_full_state_hold_supported_probability
summary.denominators.h2_genetic_baseline_eligible_probability
summary.denominators.screen_supported_probability
summary.high_branch_allele_frequency.mean
summary.high_branch_H_alpha.mean
```
