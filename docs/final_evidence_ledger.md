# Final finite-model evidence ledger

## Completion scope

This ledger closes the current finite-model campaign. It does **not** declare a
universal ecological or evolutionary law, and it does not estimate a biological
mutation rate. Every numerical result below is Type S evidence for the declared
finite trait-recruitment, symmetric allele-mutation, full-state-transfer, and
fragmentation closures.

The campaign's causal pathway is

\[
\text{patch size}
\rightarrow \text{interaction intensity}
\rightarrow \text{trait-space topology}
\rightarrow \text{local effective size}
\rightarrow \text{genetic diversity}.
\]

## Canonical evidence status

| item | finite-model question | final status | boundary of the claim |
|---|---|---|---|
| H1 | Can a mutation-conditioned high full state retain interaction memory under the declared continuation closure? | **supported, Type S** | selected mutation-H1 primary domain; not a theorem or biological estimate |
| H3 | Does equal isolation from that full state lower interaction, local effective size, and realised high-trait mass? | **supported, Type S** | conditionally on an H1-prepared source and conservation-preserving projection |
| H2-A | Does fixed absolute \(H_\alpha\le0.20\) or \(H_\gamma\le0.20\) reliably precede realised trait loss? | **not retained as a robust absolute-warning rule** | mixed ordering in one selected finite domain; no global truth value assigned |
| H2-R | Under the locked deterioration schedule, does baseline-relative diversity erosion precede realised trait loss? | **supported, Type S** | one calibration-selected configuration only; all non-events stay censored |

## H1 and H3

The mutation-primary H1/H2/H3 chain (Actions run `28456092898`) kept the
high-state source, the equal-isolated projection, and all conservation
invariants explicit. Conditional on H1 source preparation, equal isolation
lowered final interaction, local effective size, and realised high-trait mass
relative to one large patch. This remains the canonical H3 result.

The source domain used by the later H2 work was selected from the independent
mutation-H1 validation (Actions run `28436777080`), not by inspecting H2
warning outcomes.

## H2-R: retained conditional result

Ramp-and-hold trait-loss-only calibration (Actions run `28496735824`) examined
12 frozen mutation-H1 primary cells without calculating genetic-warning values.
It selected exactly one configuration because every calibration seed block had
post-baseline realised trait-loss frequency in the predeclared 0.30--0.70
interval:

```text
mutation rate = 0.10
A_ref = 0.8
kappa = 6.0
equal-isolated landscape
barrier schedule = ramp 30 generations + hold 90 generations
normalized total barrier increase = 0.15
calibration seed-block trait-loss frequencies = 0.50, 0.40, 0.40, 0.50, 0.50
```

The independent validation (Actions run `28500796310`) used fresh seeds
`20261110`--`20261114` and 20 replicates per seed. Of 100 attempted sources,
83 produced an H1-prepared, projection-supported trajectory; 35 of those 83
had post-baseline realised trait loss. For each of the six predeclared H2-R
endpoints—\(H_\alpha\) and \(H_\gamma\), each at relative declines
\(r=0.05,0.10,0.20\)—all 35 valid same-replicate pairs had warning before trait
loss, with zero ties and zero lags. The remaining 48 available trajectories did
not reach realised trait loss in the scheduled horizon and remain right-censored.

This supports only the following conditional statement:

\[
\tau_{\Delta H_x(r)} < \tau_T,
\qquad x\in\{\alpha,\gamma\},\quad r\in\{0.05,0.10,0.20\},
\]

for the selected configuration and observed event pairs. It does not generalize
to the other eleven calibration-unselected cells.

## H2-A: fixed-threshold secondary audit and retirement

H2-A existed before H2-R and used the fixed thresholds

\[
H_\alpha\le0.20,\qquad H_\gamma\le0.20.
\]

No new simulation, parameter search, threshold change, cell selection, or
schedule selection was performed for H2-A. The raw diversity series from the
already completed independent H2-R validation were reanalysed at exactly these
pre-existing thresholds. The result is stored in
`docs/evidence/h2a_fixed_threshold_secondary_audit_selected_domain_20261110_14.csv`
and its seed-block counterpart.

| fixed warning | trajectories available | warning observed | trait loss observed | valid pairs | lead | tie | lag | censored |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| \(H_\alpha\le0.20\) | 83 | 42 | 35 | 20 | 14 | 0 | 6 | 63 |
| \(H_\gamma\le0.20\) | 83 | 28 | 35 | 16 | 8 | 0 | 8 | 67 |

Observed lags preclude treating either fixed threshold as a reliable canonical
early-warning rule in this selected finite closure. This is **not** a proof that
H2-A is false in all models or biological systems; its global truth value remains
unassigned. The current repository therefore retires H2-A from the canonical
claim set rather than tuning mutation, horizon, barrier, or threshold until it
appears to work.

## Repository closure rule

The canonical claims of this repository are now H1, H3, and the conditional
H2-R result above. The H2-A audit is archived as a negative robustness result.
Any future attempt to alter the biological closure, use a different mutation
model, change the absolute threshold, or search new deterioration schedules
belongs in a new extension project with a separately declared protocol. It must
not overwrite or silently revise this ledger.
