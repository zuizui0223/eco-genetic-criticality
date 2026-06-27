# H3: network migration-matrix bounds

## Declared migration update

For a network of patches, write the allele-frequency update as

\[
p'_{i}=\sum_j M_{ij}p_j,
\]

where \(M_{ij}\) is the contribution of source patch \(j\) to destination patch
\(i\).  The matrix is **destination-by-source row-stochastic**:

\[
M_{ij}\ge0,
\qquad \sum_j M_{ij}=1.
\]

This form permits asymmetric source-sink networks, stepping-stone corridors,
and distance-derived dispersal kernels.  The previous complete-graph update

\[
p'_i=(1-m)p_i+m\sum_jw_jp_j
\]

is the special case

\[
M_{ij}=(1-m)\mathbf{1}_{i=j}+mw_j.
\]

## Common-floor theorem

If every source patch satisfies \(p_j\ge p_{\min}\), then for every
destination,

\[
p'_i=\sum_jM_{ij}p_j\ge\sum_jM_{ij}p_{\min}=p_{\min}.
\]

Thus migration cannot lower a common floor before finite sampling.  This is an
exact H3 theorem for the declared matrix update.  It does **not** show that
migration raises a patch that lies below the floor, and it does not prove
long-run persistence under drift, demographic extinction, or trait loss.

## Focal rescue condition

Let \(b_j\) be an independently established lower bound for the source
frequency in patch \(j\).  A focal destination \(i\) has the sharp migration
lower bound

\[
p'_i\ge \sum_jM_{ij}b_j.
\]

For a desired target \(p_{\rm target}\), migration rescue is certified only if

\[
\sum_jM_{ij}b_j\ge p_{\rm target}.
\]

This separates a real rescue condition from the much weaker statement that
migration is nonzero.  A poorly connected patch, or one connected mainly to
low-frequency sources, may fail the condition even when the network contains
other high-frequency patches.

## Scope boundary

`network_migration_matrix_theory.py` concerns allele-frequency composition
only.  The following remain separate H3 modelling tasks:

- migration of individuals and demographic abundance;
- dispersal of trait-bin propagules and their recruitment;
- local extinction and recolonisation;
- coupling network topology to interaction support, local \(N_e\), and realised
  high-trait occupancy.

Those components must be given their own life-cycle update before a claim about
demographic rescue or trait rescue is made.
