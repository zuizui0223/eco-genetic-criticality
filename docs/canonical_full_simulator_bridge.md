# Canonical H1 as a strict full-simulator limit

The canonical H1 theorem uses the one-state map

\[
q_{t+1}=\operatorname{sigmoid}\left[\kappa\left(\frac{A}{A_{\rm ref}}q_t-\theta\right)\right].
\]

The full multipatch simulator has additional density, trait, allele, and
population states.  It must therefore not be described as the canonical map in
general.  This bridge identifies a parameter limit in which its interaction
coordinate is exactly the canonical map.

## Declared embedding

For one patch, set:

\[
D_t=1,
\qquad S_t=q_t,
\]

where \(D_t\) is the full simulator density factor and \(S_t\) its interaction
support signal.  In code, this is obtained by:

- starting census size at carrying population and setting baseline growth so it
  remains there;
- removing trait and allele terms from `interaction_support_signal`;
- setting its interaction coefficient to one; and
- using no migration.

The full update then becomes

\[
q_{t+1}=\operatorname{sigmoid}\left[\kappa\left(\frac{A}{A_{\rm ref}}D_tS_t-\theta\right)\right]
=\operatorname{sigmoid}\left[\kappa\left(\frac{A}{A_{\rm ref}}q_t-\theta\right)\right].
\]

`canonical_full_simulator_bridge_certificate` runs both recursions and checks
that every generation agrees within a declared numerical tolerance.

## Scope

This is a consistency bridge, not a theorem for the unrestricted full model.
Once density varies, trait or allele feedback contributes to support, or patches
are coupled by migration, the canonical proof no longer automatically applies.
Those departures are precisely the axes for the subsequent finite stochastic
phase diagrams.
