# Independent H2-R relative-warning validation

## Why this validation has one configuration

The ramp-and-hold trait-loss-only calibration v2 (Actions run `28496735824`) was
run across all 12 frozen mutation-H1 primary cells without calculating or saving
H-alpha, H-gamma, relative-warning times, warning leads, or lead times. Exactly
one configuration satisfied the predeclared requirement that every calibration
seed block have conditional post-baseline realised trait-loss probability in
\([0.30,0.70]\):

```text
mutation rate = 0.10
A_ref = 0.8
kappa = 6.0
schedule = ramp 30 generations + hold 90 generations
total normalized barrier increase = 0.15
calibration trait-loss rates by seed block = 0.50, 0.40, 0.40, 0.50, 0.50
```

The other eleven primary cells remain calibration-unselected. This validation
makes **no** claim about them. It does not rerun the schedule selection or scan
another cell after inspecting diversity outcomes.

## Fresh validation state path

```text
frozen selected calibration domain
-> new-seed nested mutation-H1 calibration
-> high-route full-state replay and hold
-> conservation-preserving equal-isolated projection
-> locked 30-generation ramp + 90-generation hold schedule
-> relative H-alpha/H-gamma first-passage comparisons
```

The validation workflow uses fresh master seeds `20261110`–`20261114` and 20
replicates per seed. Calibration seeds (`20261010`–`20261014`) are never reused
for this outcome.

## Endpoints all reported

For each trajectory, each diversity measure
\(x\in\{\alpha,\gamma\}\), and each predeclared decline fraction
\(r\in\{0.05,0.10,0.20\}\), define

\[
\tau_{\Delta H_x(r)}
=\inf\{t>0:H_x(t)\le (1-r)H_x(0)\}.
\]

The endpoint compares this time with post-baseline realised high-trait loss
\(\tau_T\):

\[
\tau_{\Delta H_x(r)}<\tau_T.
\]

No endpoint is selected after execution. The output gives every one of the six
`(diversity metric, relative decline)` endpoints, recording:

- trajectory availability after H1 source reconstruction and projection;
- baseline eligibility;
- warning observed;
- trait loss observed;
- valid same-replicate pair;
- warning lead, tie, or lag; and
- seed-block-specific pair and lead counts.

A missing warning or trait-loss event is retained as censored. A tie is not a
lead. A baseline crossing is not an early warning because relative warnings are
evaluated only for \(t>0\).

## Scope of conclusions

A successful lead fraction for an endpoint is finite **Type S** evidence for
this selected mutation, landscape, full-state transfer, and ramp-and-hold
closure. It does not prove a universal genetic early-warning law, validate H2-A
(the fixed absolute-threshold proposition), or generalize to the 11 cells that
were not selected by trait-loss calibration.
