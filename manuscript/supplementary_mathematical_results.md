# Supplementary mathematical results

This supplement restates existing analytical material in manuscript order. It
adds no new theorem. Every result below carries the assumptions under which it
is valid.

## S1. Canonical sigmoid fixed-point geometry

### S1.1 Statement

Consider

\[
f(q)=\operatorname{sigmoid}\!\left[
\kappa\left(\frac{A}{A_{\rm ref}}q-\theta\right)
\right],
\qquad q\in[0,1],
\]

and define \(K=\kappa A/A_{\rm ref}\). Fixed points satisfy

\[
F(q)=\operatorname{logit}(q)-Kq+\kappa\theta=0.
\]

**Result.** If \(K\le4\), there is a unique fixed point. If \(K>4\), set

\[
q_-=\frac{1-\sqrt{1-4/K}}{2},\qquad
q_+=\frac{1+\sqrt{1-4/K}}{2}.
\]

Then strict three-fixed-point geometry occurs exactly when

\[
\theta_-<\theta<\theta_+,
\]

where

\[
\theta_-=\frac{Kq_- - \operatorname{logit}(q_-)}{\kappa},
\qquad
\theta_+=\frac{Kq_+ - \operatorname{logit}(q_+)}{\kappa}.
\]

The outer fixed points are locally stable and the middle fixed point is unstable.

### S1.2 Proof

Differentiating gives

\[
F'(q)=\frac{1}{q(1-q)}-K.
\]

Since \(q(1-q)\le1/4\) for \(q\in(0,1)\),

\[
\frac{1}{q(1-q)}\ge4.
\]

Thus, when \(K\le4\), \(F'(q)\ge0\) everywhere in \((0,1)\). Together with
\(F(q)\to-\infty\) as \(q\downarrow0\) and \(F(q)\to+\infty\) as
\(q\uparrow1\), monotonicity yields exactly one zero.

When \(K>4\), the equation \(F'(q)=0\) has exactly the two roots \(q_-<q_+\)
shown above. Therefore \(F\) has one local maximum and one local minimum. The
values of \(\theta\) for which the local extrema straddle zero are exactly the
open interval \((\theta_-,\theta_+)\); in that interval, continuity gives three
zeros. At either endpoint, one double root occurs, yielding a saddle-node rather
than strict bistability.

At a fixed point,

\[
f'(q)=Kq(1-q).
\]

The outer roots lie in regions where \(Kq(1-q)<1\), whereas the middle root lies
where \(Kq(1-q)>1\). Hence the outer roots are locally stable and the middle root
is unstable. \(\square\)

### S1.3 Scope

This is a theorem for the specified one-dimensional map. It does not establish
bistability for a finite-bin trait model, for a model with density feedback, or
for an empirical interaction system.

---

## S2. Branch-dependent potential high-trait viability

Let \(Z_H\) denote the declared high-trait region and let

\[
m_H(q)=\max_{z\in Z_H}[W(z;q)-\tau]
\]

be the corresponding viability margin under interaction state \(q\).

**Conditional result.** Suppose the canonical map has stable states
\(q_{\rm low}<q_{\rm high}\). If

\[
m_H(q_{\rm low})<0<m_H(q_{\rm high}),
\]

then the high-trait region is potentially unavailable at \(q_{\rm low}\) and
potentially viable at \(q_{\rm high}\).

**Justification.** The sign condition is exactly the definition of the margin's
availability state at the two declared interaction branches. No statement about
finite realised occupancy follows without a recruitment, inheritance, and
demographic closure.

---

## S3. A sufficient no-bistability certificate

Let \(g\) be a differentiable response function with

\[
M=\sup_x|g'(x)|<\infty,
\]

and consider

\[
q_{t+1}=g\{\kappa(Aq_t-\theta)\}.
\]

**Theorem (contraction certificate).** If

\[
\kappa A M<1,
\]

then the map is a contraction on its state interval and therefore has a unique
fixed point. In particular, strict bistability is impossible.

**Proof.** The derivative magnitude is bounded by \(\kappa A M<1\), so the map
is a contraction. Banach's fixed-point theorem gives a unique fixed point.
\(\square\)

The converse is not asserted: failure of this inequality does not by itself
prove bistability.

---

## S4. Finite transmission variance identity

Let \(P'\) denote an unbiased post-selection transmission frequency conditional
on parental frequency \(p^*\), so that

\[
E(P'\mid p^*)=p^*.
\]

Let heterozygosity be \(H(p)=2p(1-p)\).

**Theorem.** If \(\operatorname{Var}(P'\mid p^*)>0\), then

\[
E\{H(P')\mid p^*\}=H(p^*)-2\operatorname{Var}(P'\mid p^*)<H(p^*).
\]

**Proof.** By the variance identity,

\[
E\{P'^2\mid p^*\}=\operatorname{Var}(P'\mid p^*)+(p^*)^2.
\]

Hence

\[
\begin{aligned}
E\{H(P')\mid p^*\}
&=2E\{P'-P'^2\mid p^*\}\\
&=2\left[p^*-(p^*)^2-\operatorname{Var}(P'\mid p^*)\right]\\
&=H(p^*)-2\operatorname{Var}(P'\mid p^*).
\end{aligned}
\]

Positive conditional variance gives the strict inequality. \(\square\)

The theorem makes no claim about how interaction state or patch area affects the
variance; those links are closure-dependent.

---

## S5. Row-stochastic migration bounds

Let

\[
p'_i=\sum_jM_{ij}p_j,
\]

with \(M_{ij}\ge0\) and \(\sum_jM_{ij}=1\).

### S5.1 Common floor

**Theorem.** If \(p_j\ge p_{\min}\) for every source patch \(j\), then
\(p'_i\ge p_{\min}\) for every destination patch \(i\).

**Proof.**

\[
p'_i=\sum_jM_{ij}p_j
\ge\sum_jM_{ij}p_{\min}
=p_{\min}.
\]

\(\square\)

### S5.2 Focal rescue lower bound

**Theorem.** If \(p_j\ge b_j\) for source-specific lower bounds \(b_j\), then

\[
p'_i\ge\sum_jM_{ij}b_j.
\]

Consequently, deterministic frequency rescue at a target \(p_{\rm target}\) is
certified only if

\[
\sum_jM_{ij}b_j\ge p_{\rm target}.
\]

**Proof.** Replace each \(p_j\) in the update with its lower bound and use
non-negativity of \(M_{ij}\). \(\square\)

Neither theorem includes finite sampling, demographic abundance, extinction,
recolonisation, or trait recruitment.

---

## S6. Conditional patchwise non-additivity

Suppose a declared ecological closure requires \(A_j>A_c\) in every patch to
maintain an interaction-supported trait mode. A fixed total area
\(A_{\rm total}=\sum_jA_j\) does not guarantee this condition after subdivision:
for any \(m\) such that \(A_{\rm total}/m\le A_c\), equal subdivision produces
patches that fail the per-patch requirement. This elementary implication is
conditional on the ecological threshold assumption and does not itself identify
\(A_c\).

---

## S7. First-passage convention used in finite results

For a trajectory-specific diversity series \(H_x(t)\), relative warning time is

\[
\tau_{\Delta H_x(r)}=
\inf\{t>0:H_x(t)\le(1-r)H_x(0)\}.
\]

Let \(\tau_T\) denote first post-baseline realised high-trait loss. A warning
lead is defined only when both event times are observed, in which case it is the
strict inequality

\[
\tau_{\Delta H_x(r)}<\tau_T.
\]

If either event is unobserved in the finite horizon, the pair is censored. A tie
is not a lead. This is a reporting convention and finite-study design principle,
not a theorem that the ordering holds.
