# Canonical H1 bifurcation certificate

## Declared system

This document concerns only the canonical one-state interaction reduction

\[
q_{t+1}=f(q_t)=\operatorname{sigmoid}\left[\kappa\left(\frac{A}{A_{\rm ref}}q_t-\theta\right)\right],
\qquad q_t\in[0,1].
\]

It holds density at one and removes the additional realised-trait and allele
feedback terms used by the finite-bin multipatch closure.  Thus it is a
specified-system theorem for H1, not a theorem about every positive-feedback
ecology and not a direct theorem for the full simulator.

Write

\[
K=\kappa A/A_{\rm ref}.
\]

Fixed points solve

\[
F(q)=\operatorname{logit}(q)-Kq+\kappa\theta=0.
\]

## Exact fixed-point geometry

\[
F'(q)=\frac{1}{q(1-q)}-K.
\]

Because \(1/[q(1-q)]\ge4\), if \(K\le4\), then \(F\) is nondecreasing on
\((0,1)\), so the canonical map has one fixed point.  In particular,
\(K>4\) is necessary before bistability is possible.

For \(K>4\), there are two turning points,

\[
q_- = \frac{1-\sqrt{1-4/K}}{2},
\qquad
q_+ = \frac{1+\sqrt{1-4/K}}{2}.
\]

There are exactly three fixed points if and only if the barrier lies strictly
inside

\[
\theta_- = \frac{Kq_- - \operatorname{logit}(q_-)}{\kappa}
< \theta <
\theta_+ = \frac{Kq_+ - \operatorname{logit}(q_+)}{\kappa}.
\]

At either endpoint there is a saddle-node; the endpoint is deliberately not
reported as strict bistability.

At a fixed point, the map multiplier is

\[
f'(q)=Kq(1-q).
\]

Therefore the low and high fixed points are locally stable and the middle fixed
point is unstable.  The module `canonical_h1_bifurcation.py` implements these
relations and returns the fixed points with their local stability labels.

## Specified-system H1 certificate

Let \(m_H(q)\) be the high-trait viability margin from the declared trait
performance surface.  The canonical H1 mechanism is certified only when all
three requirements hold:

1. \(K>4\) and \(\theta_-<\theta<\theta_+\), giving low and high stable
   interaction branches;
2. \(m_H(q_{\rm low})<0\);
3. \(m_H(q_{\rm high})>0\).

Then the potential high-trait component is branch dependent: a high-investment
mode is unavailable on the low branch but viable on the high branch.  This is
the canonical realization of P1 after the existence of both branches has been
established analytically.

## Hysteresis check

`follow_barrier_path` carries the terminal interaction state from one barrier
to the next.  With \(\kappa=8\), \(A/A_{\rm ref}=1\), the strict bistable
interval is approximately

\[
0.3667900062 < \theta < 0.6332099938.
\]

A rising barrier path started on the high branch remains high at
\(\theta=0.5\), whereas a descending path started on the low branch remains
low at the same barrier.  The two paths jump at opposite saddle-node
boundaries.  The associated test is a reproducible numerical confirmation of
the analytic branch geometry; it is not independent evidence for the full
finite-bin system.

## What remains outside this certificate

The following still require separate analysis or simulation:

- whether density feedback preserves the canonical reduction;
- whether realised high-trait abundance, rather than potential viability,
  changes discontinuously;
- how finite recruitment, genetic drift, and allele-linked feedback change the
  bistable region;
- whether a biological interaction system supplies the required state equation
  and parameter values.
