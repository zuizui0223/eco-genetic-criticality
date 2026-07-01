# H2 reframing and final evidence status

## What remains fixed

Nothing here changes H1, H3, the symmetric-mutation closure, or the full-state
H1 transfer rule. H2-A and H2-R remain distinct propositions. No finite result
is promoted to a theorem or a biological estimate.

## H2-A: original fixed-threshold proposition

The original H2 statement is retained as **H2-A**:

> Under declared conditions, fixed absolute \(H_\alpha\) or \(H_\gamma\)
> warning thresholds can precede realised high-trait loss.

The fixed thresholds are

\[
H_\alpha\le0.20,\qquad H_\gamma\le0.20.
\]

The stationary mutation-primary chain left H2-A unresolved because most records
were right-censored. Later, the completed independent H2-R validation supplied
raw diversity series in a configuration selected *without* genetic-warning
outcomes. A no-resimulation secondary audit applied exactly these pre-existing
H2-A thresholds to those raw series:

| fixed warning | valid pairs | lead | tie | lag | censored |
|---|---:|---:|---:|---:|---:|
| \(H_\alpha\le0.20\) | 20 | 14 | 0 | 6 | 63 |
| \(H_\gamma\le0.20\) | 16 | 8 | 0 | 8 | 67 |

Thus H2-A is **not retained as a robust absolute-warning rule** in the selected
finite closure. Observed lags prevent a uniform ordering claim. This is not a
proof that H2-A is false in every model or biological system, so its global
truth value is not assigned. The repository deliberately does not tune the
threshold, mutation, barrier, or horizon after this result.

## H2-R: conditional relative genetic warning

**H2-R** is a separate dynamic hypothesis. Let \(H_x(t)\) denote
\(H_\alpha(t)\) or \(H_\gamma(t)\), for \(x\in\{\alpha,\gamma\}\), along a
trajectory starting from a polymorphic H1-conditioned full state. Under a
predeclared monotone decline in interaction support, define

\[
\tau_{\Delta H_x(r)}
=
\inf\left\{t>0:H_x(t)\le(1-r)H_x(0)\right\},
\qquad r\in\{0.05,0.10,0.20\},
\]

and let \(\tau_T\) be first realised high-trait loss. H2-R asks whether

\[
\tau_{\Delta H_x(r)}<\tau_T.
\]

A baseline time is never called a warning. Missing warning or trait-loss events
remain censored, and a tie is not a lead.

## Selection before warning outcomes

H2-R schedule calibration measured realised trait loss only. The initial
linear-ramp family selected no cell. A second predeclared ramp-and-hold family
kept the same normalized total barrier increases while allowing time at the
final barrier. In Actions run `28496735824`, exactly one of the 12 frozen
mutation-H1 primary cells met the all-seed-block trait-loss availability rule:

```text
mutation rate = 0.10
A_ref = 0.8
kappa = 6.0
equal-isolated landscape
ramp 30 + hold 90 generations
normalized barrier increase = 0.15
calibration seed-block trait-loss frequencies = 0.50, 0.40, 0.40, 0.50, 0.50
```

The other eleven cells were not carried into warning validation. The selection
neither calculated nor inspected H-alpha, H-gamma, warning lead counts, or lead
times.

## Final H2-R result

Independent validation used fresh seeds `20261110`--`20261114`, 20 replicates
per seed, and the locked configuration above (Actions run `28500796310`). Of
100 attempted H1 sources, 83 yielded an available projected trajectory; 35 of
those had realised trait loss in the scheduled horizon. For each of the six
predeclared endpoints—\(H_\alpha\) and \(H_\gamma\) at each
\(r=0.05,0.10,0.20\)—all 35 valid pairs had warning before realised trait loss,
with zero ties and zero lags. The other 48 available trajectories are
right-censored because trait loss was not observed in horizon.

H2-R is therefore **supported as Type S evidence only in this selected
configuration and for observed event pairs**. It does not establish H2-A, a
universal early-warning rule, or a conclusion about the eleven
calibration-unselected cells.

The complete final status is recorded in `docs/final_evidence_ledger.md`.
