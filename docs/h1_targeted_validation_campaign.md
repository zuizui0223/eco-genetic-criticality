# One-large canonical-H1-targeted validation campaign

The ordinary `quick`, `standard`, and `full` profiles use Cartesian grids of
raw interaction barriers. That is appropriate for broad exploratory phase
diagrams, but it is not a controlled test of the H1 mechanism across changing
area reference and feedback strength.

The canonical interaction map is

```text
q_next = sigmoid(kappa * ((A / A_ref) * q - theta)).
```

Therefore the same raw barrier `theta` has a different location relative to the
bistable interval when `A / A_ref` changes. In particular, `one_large` has
local area equal to total area, while equal fragments have local area equal to
total area divided by patch count.

This runner targets **one-large** canonical H1 before testing what happens after
fragmentation.

## Design

For each declared `(area_reference, interaction_feedback)` pair, derive the
strict one-large canonical bistable interval

```text
(theta_lower, theta_upper).
```

For each declared inside position `r`, where `0 < r < 1`, set

```text
theta = theta_lower + r * (theta_upper - theta_lower).
```

The default positions are `0.25`, `0.50`, and `0.75`. Thus every designed cell
is strictly canonical-bistable for one-large interaction dynamics. The same raw
barrier is then applied unchanged to `equal_isolated` and `equal_migrating`;
that is the intended spatial comparison.

The design criterion is **only** canonical one-large bistability. It does not
filter cells according to finite branch retention, trait viability, genetic
warning order, fragmentation response, migration response, or whether a
first-passage pair is censored.

## Execution

A small installation and artifact-shape check:

```bash
python scripts/run_h1_targeted_validation_campaign.py \
  --profile quick \
  --replicates 2 \
  --generations 8 \
  --inside-position 0.5 \
  --output-dir artifacts/h1_targeted_validation \
  --prefix smoke
```

A moderate targeted study:

```bash
python scripts/run_h1_targeted_validation_campaign.py \
  --profile standard \
  --output-dir artifacts/h1_targeted_validation \
  --prefix standard_v1
```

The full profile is opt-in. It should not be used in normal CI because it
launches one finite H1--H3 validation subcampaign for every
`(area_reference, interaction_feedback)` pair, with every declared inside
position.

### GitHub Actions run

The **H1-Targeted Finite Validation** workflow is manual (`workflow_dispatch`).
Use it for a run whose artifacts need to be retained with its exact code
revision:

1. Open **Actions** and choose **H1-Targeted Finite Validation**.
2. Use `quick` only as an end-to-end smoke check. Use `standard` for the first
   primary finite validation run.
3. Leave inside positions as `0.25,0.5,0.75` for the predeclared three-point
   within-interval design, unless a narrower sensitivity run is explicitly
   declared.
4. Record the master seed in the analysis log; the default is `20260629`.
5. Leave replicate and generation overrides blank to use the named profile, or
   set both deliberately and treat that run as a separately declared design.
6. Select `full` only after explicitly enabling the full-profile approval
   switch. It has a longer timeout because it evaluates many more finite runs.

The workflow uploads every manifest, design map, ledger, and subcampaign raw
artifact for 90 days. It also attempts to upload partial outputs on failure, so
an interrupted run is inspectable rather than silently discarded.

## Outputs

At the top level, the runner writes:

- `*.design.json`: every normalized target and its raw barrier;
- `*.ledger.csv` and `*.ledger.json`: all cells aligned across H1, H2, and H3;
- `*.manifest.json`: base spec, interval-derived design, seed rule, arguments,
  and no-selection policy; and
- `*.subcampaigns/`: complete raw artifacts from the unified finite validation
  campaign for each fixed `(A_ref, kappa)` pair.

Each subcampaign receives a deterministic pair-specific master seed. This avoids
silently reusing identical stochastic streams across different feedback/area
reference pairs while preserving reproducibility.

## Interpretation

One-large canonical bistability is a Type T statement about the declared
one-state reduction. Finite branch separation, hysteresis continuation, genetic
warning order, fragmentation response, and migration modulation are Type S
results for the finite closure.

The targeted design makes the *entry condition* for H1 controlled. It does not
make any finite result automatic. Equal fragments may lie outside their own
local canonical bistable regime, finite feedback may erase or alter branch
separation, and first-passage comparisons can remain censored. Those are the
results to be measured, not reasons to remove cells.
