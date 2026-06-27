# Theorem-boundary phase-diagram artifacts

`theorem_boundary_phase_diagram.py` runs the same landscape scenarios,
parameter cells, random-seed schedule, and first-passage summaries as the
existing multipatch experiment framework.  It adds a theorem-boundary audit to
every replicate.

## Why a separate artifact

The canonical H1 theorem applies only to a one-state interaction map.  The full
simulator may add density dependence, realised trait feedback, allele feedback,
finite inheritance, multiple patches, and migration.  A phase diagram that
reports only outcomes risks obscuring which rows remain in the theorem limit and
which are controlled departures from it.

Each JSON replicate contains:

- the ordinary outcome summary, including censored first-passage times;
- the maximum and mean canonical-update residual;
- density and support deviations from the canonical H1 limit; and
- named departure labels.

Each flat CSV cell contains summary fields including:

- `scope.patchwise_canonical_update_probability`;
- `scope.single_patch_canonical_theorem_limit_probability`;
- residual and deviation summaries; and
- one probability per departure label.

## Reading the result

A cell with theorem-limit probability one is a verified deterministic special
case of canonical H1.  A cell with nonzero residual is not a failed theorem.
It is a deliberately richer model.  The departure labels identify which
mechanisms are responsible, so the plotted outcome is interpreted as a
robustness result rather than as an analytic proof.

## Minimal workflow

1. Run a quick profile to confirm artifacts and labels.
2. Run a standard profile to explore interaction feedback, barriers, landscapes,
   and migration.
3. Use a full profile only for final figures after thresholds and comparison
   targets are predeclared.
4. Plot outcome metrics next to scope metrics; do not collapse theorem-limit and
   departure cells into a single unlabeled claim.
