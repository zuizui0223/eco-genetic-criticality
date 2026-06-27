# H1 theorem-boundary audit

The canonical H1 proof applies to

\[
q_{t+1}=\operatorname{sigmoid}\left[\kappa\left(\frac{A}{A_{\rm ref}}q_t-\theta\right)\right].
\]

The full simulator instead updates each patch by

\[
q_{j,t+1}=\operatorname{sigmoid}\left[\kappa\left(\frac{A_j}{A_{\rm ref}}D_{j,t}S_{j,t}-\theta\right)\right],
\]

where \(D\) is density and \(S\) can include interaction memory, realised
high-trait mass, and allele frequency.  The theorem-boundary audit reports the
absolute residual between these two updates at every patch and generation.

## What the audit certifies

`patchwise_canonical_update_certified=True` means the observed full-simulator
updates equal the canonical patchwise update within the declared numerical
tolerance.  It does not by itself imply a one-state theorem, because multiple
patches or migration can still be present.

`single_patch_canonical_theorem_limit_certified=True` additionally requires one
patch and no migration.  This is the appropriate flag for carrying the canonical
H1 theorem into the full code base.

## Departure labels

The audit records named reasons why a run is outside the strict canonical limit:

- `density_not_one`
- `trait_feedback_enabled`
- `allele_feedback_enabled`
- `support_not_equal_interaction`
- `migration_enabled`
- `multiple_patches`

These labels are intended to accompany phase-diagram rows.  A nonzero residual
is not a model failure.  It identifies the biological mechanism that has been
added beyond the theorem, and therefore the mechanism whose robustness should
be assessed empirically or with replicated simulation.
