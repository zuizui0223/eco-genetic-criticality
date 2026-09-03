# Common scalar ecological-state theorem

Status: exact finite representation theorem plus a direct witness from the locked H3 fragmentation-gradient evidence.

## Question

The manuscript already shows that potential viability, realised occupancy, interaction, effective size, diversity and allele persistence need not move together. That observation alone can sound obvious: different variables can have different trajectories.

The stronger question is:

> **When can all declared ecological responsibilities nevertheless be represented exactly by one scalar state whose increase has one coherent meaning?**

This document gives a necessary-and-sufficient answer for finite target sets.

## Definition — directionally coherent sufficient scalar state

Let `Omega` be a finite set of ecological states and let

\[
T=(T_1,\ldots,T_m):\Omega\to\mathbb R^m
\]

be the declared target vector, with every coordinate oriented so that larger means no worse for that responsibility.

A scalar state

\[
h:\Omega\to\mathbb R
\]

is **directionally coherent and sufficient** for `T` when for each target there is a nondecreasing function `g_j` such that

\[
T_j=g_j\circ h.
\]

Thus equal scalar states must have identical target vectors, and moving upward in the scalar state may not improve one declared responsibility while worsening another.

This is deliberately stronger than arbitrary lossless coding into a real number. An injective numeric label can encode any finite state set but has no common ecological order. The theorem concerns a scalar **state/health axis**, not arbitrary enumeration.

## Product order

For target vectors `u,v in R^m`, write

\[
u\preceq v
\]

when `u_j<=v_j` for every target `j`.

A collection of target vectors is a **chain** when every pair is comparable under this coordinatewise product order.

A **crossing pair** is a pair of states `a,b` for which one target prefers `a` while another prefers `b`; equivalently, neither `T(a)<=T(b)` nor `T(b)<=T(a)` componentwise.

## Theorem S1 — exact existence criterion for one monotone scalar state

A directionally coherent sufficient scalar state exists **if and only if** the distinct target vectors

\[
\{T(\omega):\omega\in\Omega\}
\]

form a chain under product order.

### Proof — necessity

Assume such a scalar `h` exists. Take any two states `a,b`.

If `h(a)=h(b)`, factorization gives

\[
T_j(a)=g_j(h(a))=g_j(h(b))=T_j(b)
\]

for every `j`, so the target vectors are equal.

If `h(a)<h(b)`, monotonicity of every `g_j` gives

\[
T_j(a)\le T_j(b)
\]

for every target. Hence `T(a)<=T(b)` componentwise. The case `h(b)<h(a)` is symmetric.

Therefore every pair of target vectors is comparable: the image of `T` is a chain. ∎

### Proof — sufficiency

Assume the distinct target vectors form a finite chain. Order those distinct vectors as

\[
v^{(0)}\prec v^{(1)}\prec\cdots\prec v^{(r)}
\]

under product order. Define

\[
h(\omega)=k\quad\text{whenever }T(\omega)=v^{(k)}.
\]

For each target coordinate define

\[
g_j(k)=v^{(k)}_j.
\]

Because the vectors are ordered componentwise, every `g_j` is nondecreasing. By construction,

\[
T_j(\omega)=g_j(h(\omega))
\]

for every state and target. Thus `h` is directionally coherent and sufficient. ∎

## Corollary S1a — one crossing pair is a certificate of impossibility

If there exist states `a,b` and targets `i,j` such that

\[
T_i(a)>T_i(b),
\qquad
T_j(a)<T_j(b),
\]

then no directionally coherent sufficient scalar state exists for the target family.

This is stronger than saying the two variables have different slopes. It proves that **no single monotone ecological health/state axis can exactly preserve both responsibilities on those states**.

## Corollary S1b — a scalar index can exist on a restricted domain and fail after domain expansion

If a subset of states has target vectors forming a chain, it admits a common monotone scalar. Adding one state that creates a crossing pair destroys that possibility. Therefore scalar adequacy is domain-relative: success on one fragmentation range does not license reuse after the range is extended.

## Locked H3 gradient supplies a direct crossing certificate

The preregistered fresh-seed H3 gradient reports pooled retained fractions for the same 1,037 prepared sources. Compare the two-patch and sixteen-patch states:

| state | retained interaction | retained local effective size | retained realised high-trait mass |
|---|---:|---:|---:|
| 2 patches | `0.001744` | `0.221311` | `0.282918` |
| 16 patches | `0.001244` | `0.033058` | `0.393880` |

All three targets are oriented as “higher retained value = no worse.”

Between two and sixteen patches:

- interaction is **higher at 2 patches**;
- local effective size is **higher at 2 patches**;
- realised high-trait mass is **higher at 16 patches**.

Hence the two target vectors are product-order incomparable. By Corollary S1a, no directionally coherent sufficient scalar state can exactly preserve even the pair `(interaction, realised high-trait mass)` on this observed finite-model gradient, and therefore no such scalar can preserve the larger family containing them.

This conclusion is not obtained from a fitted composite index. It follows from a locked crossing already present in the finite evidence.

## A second separation: potential viability versus realised occupancy

The same gradient provides a different kind of evidence:

- one patch: potential high-trait viability present in `1037/1037` supported outcomes;
- every tested subdivision from 2 to 16 patches: potential viability present in `0/1037`;
- realised high-trait occupancy nevertheless persists at the 30-generation endpoint in approximately `99.6–100%` of supported trajectories.

This does not by itself create an order crossing, but it shows why one state label cannot be interpreted as both potential support and realised occupancy without an explicit target map. The scalar impossibility theorem and the viability/occupancy separation address complementary failures.

## What the theorem does not claim

The theorem does **not** say:

- every useful ecological index must preserve all target values exactly;
- approximate scalar summaries are impossible;
- no target-specific scalar index can be useful;
- the five named state variables are universally sufficient for nature;
- the H3 finite-model crossing is a universal fragmentation law.

It says exactly this:

> for a declared finite set of targets all oriented in one “higher is no worse” direction, one exact monotone sufficient scalar exists iff the realised target vectors form a product-order chain.

The locked H3 gradient violates that condition for at least interaction and realised high-trait mass.

## Executable obligations

`tests/test_common_scalar_state_theorem.py` must verify:

1. the constructive scalar is valid whenever target vectors form a chain;
2. every crossing pair is rejected;
3. the theorem agrees with an independent exhaustive search over all small scalar labelings and small target tables;
4. duplicate target vectors may share one scalar label;
5. the locked H3 pooled-summary CSV contains the 2-vs-16 crossing certificate;
6. adding the 16-patch state to a domain that had a monotone interaction-only ordering demonstrates the domain-expansion failure without changing the frozen evidence.
