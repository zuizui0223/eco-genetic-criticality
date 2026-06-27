# Running theorem-boundary phase diagrams

Run commands from the repository root.  The runner writes three coordinated
artifacts: a flat CSV for plotting, full replicate-level JSON, and a manifest
that records the exact profile, parameter specification, selected scenarios,
seed, tolerance, and output names.

## Smoke test

```bash
python scripts/run_theorem_boundary_phase_diagram.py \
  --profile quick \
  --output-dir artifacts/theorem_boundary/quick
```

Use this only to verify the pipeline, file layout, and departure labels.

## Exploratory analysis

```bash
python scripts/run_theorem_boundary_phase_diagram.py \
  --profile standard \
  --output-dir artifacts/theorem_boundary/standard
```

The standard profile is the default exploratory grid.  For a cheaper targeted
run, select a landscape and make any reduced-compute choice explicit in the
manifest:

```bash
python scripts/run_theorem_boundary_phase_diagram.py \
  --profile standard \
  --scenario equal_migrating \
  --replicates 20 \
  --generations 30 \
  --master-seed 20260627 \
  --output-dir artifacts/theorem_boundary/pilot
```

## Manual GitHub Actions run

The **Theorem-Boundary Phase Diagram** workflow can be started from the
repository's Actions tab.  Select a profile, one landscape or `all`, and an
explicit master seed.  The workflow uploads the CSV, JSON, and manifest as one
named artifact, so no local Python installation is required for a standard
run.

The `full` profile is intentionally guarded: its `allow_full_profile` input
must be selected explicitly.  This avoids accidentally consuming a long GitHub
Actions job for an exploratory run.

## Final figures

```bash
python scripts/run_theorem_boundary_phase_diagram.py \
  --profile full \
  --output-dir artifacts/theorem_boundary/full \
  --master-seed 20260627
```

Use the full profile only after thresholds, focal outcomes, scenario set, and
comparison analyses have been predeclared.  Existing artifact names are never
replaced unless `--overwrite` is passed explicitly.

## Required interpretation

Read each outcome beside scope fields in the CSV or JSON:

- `scope.single_patch_canonical_theorem_limit_probability = 1` identifies the
  strict canonical H1 special case.
- Nonzero canonical residual identifies a controlled departure, not a failed
  theorem.
- `scope.departure_probabilities.*` reports which mechanisms generated that
  departure.
- Genetic lead probabilities remain conditional on non-censored event pairs;
  the manifest and JSON preserve the required event-count context.
