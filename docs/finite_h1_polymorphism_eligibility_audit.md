# H1 high-branch genetic eligibility audit

## Purpose

A full-state H1 high branch can be dynamically real and still be unusable for a
genetic-warning test. H2 asks whether genetic diversity declines before realised
trait loss. That comparison needs diversity at the instant the fragmentation
experiment begins. A state with \(p=1\), \(p=0\), \(H_\alpha=0\), or
\(H_\gamma=0\) cannot later show a new loss of heterozygosity.

This audit therefore checks eligibility before H2/H3 dynamics are run.

## Baseline rule

The H1 full-state high source is prepared through the already declared path:

```text
97-point finite H1 calibration
-> rising-route high state at interior anchor
-> full-state high hold
-> conservation-preserving projection to H3 landscapes
```

For each projected generation-0 landscape, the audit calculates
population-weighted high-allele frequency \(p\), \(H_\alpha\), \(H_\gamma\),
and \(F_{ST}\).

A record is **H2 dynamic-warning eligible** only if all hold:

```text
polymorphism_epsilon < p < 1 - polymorphism_epsilon
H-alpha > h_alpha_warning_threshold
H-gamma > h_gamma_warning_threshold
```

A record is **H3 genetic-contrast eligible** only if it is polymorphic and both
heterozygosities are positive. This permits later fragmentation/drift to alter
diversity or FST; it does not require initial FST to be positive.

If a warning threshold is already crossed at generation zero, it is labelled
`*_warning_preexisting`. It is **not** counted as a warning lead in H2.

## Interpretation

An eligibility probability of zero is informative: under the present finite
closure, the full-state H1 high branch reaches fixation before fragmentation is
introduced. H1 branch memory may then be well supported, while H2 warning order
and H3 genetic-diversity effects remain untestable in that exact closure.

The correct next step after zero eligibility is not to reinterpret fixation as
early warning. It is to declare a new ecological-genetic closure that preserves
polymorphism, for example by adding a documented mutation or balancing-selection
mechanism, and validate H1 again under that changed closure.

## Standard run

```bash
python scripts/run_finite_h1_polymorphism_eligibility_audit.py \
  --profile standard \
  --master-seed 20260630 \
  --master-seed 20260631 \
  --master-seed 20260632 \
  --master-seed 20260633 \
  --master-seed 20260634 \
  --endpoint-padding-fraction 0.5 \
  --stage-generations 30 \
  --hold-generations 30 \
  --barrier-points 25 \
  --barrier-points 49 \
  --barrier-points 97 \
  --interaction-separation-threshold 0.05 \
  --maximum-normalized-bracket-width 0.03 \
  --polymorphism-epsilon 1e-12 \
  --output-dir artifacts/h1_polymorphism_eligibility \
  --prefix standard_h1_genetic_eligibility_v1
```

CSV is a parameter-pair ledger. JSON preserves every source availability flag,
all three projected baselines, and every eligibility classification.