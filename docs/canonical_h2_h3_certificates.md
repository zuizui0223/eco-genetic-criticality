# Canonical H2 and H3 certificates

This document adds two specified-system theorems to the existing canonical H1
bifurcation result.  They are deliberately narrow.  Their value is that every
assumption is visible and every conclusion can be falsified by violating an
assumption; they are not claims about all ecological systems.

## H2: expected genetic warning can strictly precede realised trait loss

Consider the deterministic expectation map

\[
\mathbb{E}[H_{t+1}] = \left(1-\frac{1}{2N_e}\right)\mathbb{E}[H_t],
\qquad
T_{t+1}=r_TT_t,
\]

with \(N_e\ge1\), \(0<r_T\le1\), \(H_0>h_*\ge0\), and
\(T_0>t_*>0\).  Define first-passage times

\[
\tau_H=\inf\{t:\mathbb{E}[H_t]\le h_*\},
\qquad
\tau_T=\inf\{t:T_t\le t_*\}.
\]

Both trajectories are non-increasing.  Consequently their computed first
passages are exact for this declared map, and the condition

\[
\tau_H<\tau_T
\]

is a strict expected-genetic-lead theorem for the map.  If either event is
outside the predeclared finite horizon, it is censored (`None`) rather than
replaced by the terminal generation.

This does **not** assert that every stochastic replicate has the same event
order.  It says exactly when a finite-transmission expectation trajectory leads
a monotone realised-trait trajectory under constant \(N_e\).

## H3a: equal isolated subdivision can remove local support at fixed total area

Use the binary local-support closure

\[
\text{high trait supported in patch }j \iff A_j\ge A_c.
\]

For total area \(A\) divided into \(m\) equal isolated patches, if

\[
A\ge A_c \quad\text{and}\quad A/m<A_c,
\]

then the one-large landscape supports the high-trait mechanism while every
fragment lacks it.  This is a strict fixed-total-area fragmentation theorem for
this threshold closure.

It does not say that every fragmented population loses the trait: a different
local closure, unequal patch sizes, adaptive dispersal, or external subsidies
can change the premise.

## H3b: external immigration can rescue a thresholded recipient

Let \(R\) be residents remaining after survival and \(I\) be **external**
immigrants.  Under the post-arrival establishment rule

\[
R+I\ge C,
\]

where \(C\) is the establishment threshold:

- if \(R=0\), the recipient recolonises exactly when \(I\ge C\);
- if \(0<R<C\), it is demographically rescued exactly when \(I\ge C-R\).

For a high-trait component with local interaction support and threshold \(C_H\),
its post-arrival abundance must also satisfy

\[
R_H+I_H\ge C_H.
\]

The explicit separation of external arrivals from residents avoids falsely
calling self-loop transport a rescue.

## Relation to the finite simulator

The finite H3 lifecycle has demographic stochasticity, density dependence,
directed dispersal, and realised allele-copy states.  The canonical certificates
are not substitutes for that simulator.  They provide exact boundary cases that
should be recovered in appropriate deterministic limits; the phase diagram then
maps which H2/H3 event orderings occur when those ideal conditions are relaxed.
