# Coupled eco-genetic sufficient-condition theorem

The three canonical certificates can be composed only after declaring the
missing link between interaction state and effective population size.  The
minimal coupled skeleton is:

\[
A \longrightarrow q \longrightarrow N_e(q) \longrightarrow \mathbb{E}[H_t],
\]

with the realised high-trait component and fragmentation support tracked
separately.

## Declared assumptions

1. **H1 branch switch.** The canonical interaction map has two stable branches
   and its high-trait viability margin changes sign from the low branch to the
   high branch.
2. **Interaction–effective-size closure.**

   \[
   N_e(q)=N_0+\beta q,\qquad N_0>0,\quad\beta>0.
   \]

   This is an ecological modelling assumption, not a consequence of H1 alone.
3. **H2 low-branch ordering.** The specified expected-drift and trait-retention
   map has \(\tau_H<\tau_T\).
4. **H3 local support loss.** Fixed total area is divided so that the one-large
   landscape clears the local support threshold but each equal isolated fragment
   does not.

## Conclusion

Under all four assumptions, the `CoupledEcoGeneticCertificate` certifies the
following *sufficient-condition chain*:

\[
\text{fragmentation removes declared local support}
\Rightarrow \text{low interaction branch},
\]

\[
\text{low branch} \Rightarrow N_{e,L}<N_{e,H},
\]

and, for the supplied low-branch H2 map,

\[
\tau_H<\tau_T.
\]

The certificate calls this a fragmentation-to-genetic-warning chain.  It does
not establish that all fragmented landscapes enter the low branch, or that all
low-\(N_e\) populations show genetic lead.  Those conclusions would require
additional dynamical assumptions.

## Rescue as an explicit interruption

A separate H3 rescue certificate can satisfy a post-arrival demographic or
high-trait establishment threshold.  When it does, the coupled certificate
records `rescue_interrupts_fragmentation_chain=True`.  This does not negate the
fragmentation mechanism; it reports the declared external-arrival route that
can break it locally.

## Why this form matters

The theorem prevents a common inferential shortcut: observing a fragmented
landscape, low diversity, and trait loss does not itself prove a causal chain.
Every arrow above is an assumption or a proved property of a named canonical
map.  The finite stochastic simulator is then used to test how robust this
chain remains when demographic noise, finite allele copies, density dependence,
and network dispersal are restored.
