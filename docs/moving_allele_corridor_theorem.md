# Moving allele-corridor theorem

## The obstruction behind a fixed interval

The restricted H-alpha multiplier needs frequency information above one half.
A natural but incorrect next step would be to try to retain one fixed interval

```text
[p_min, p_max],  1/2 <= p_min < p_max < 1.
```

Under a uniform high-allele fitness advantage `f > 1`, the simulator selection
map

```text
s(p, f) = p f / [p f + 1-p]
```

satisfies

```text
s(p_max, f) > p_max.
```

Therefore no nontrivial fixed upper bound below one is deterministically
invariant under that directional selection. This is not a numerical failure and
not a weakness of the proof; it is an algebraic obstruction.

## Moving corridor

Let every patch at time `t` lie in

```text
L_t <= p_j,t <= U_t,
```

with interaction in `[q_min, q_max]` and next census at least `N_min`.

For monotone high-allele selection, define the deterministic post-selection
bounds

```text
a_t = s(L_t, f(q_min))
b_t = s(U_t, f(q_max)).
```

Census-weighted migration is a convex combination, so every post-migration
frequency remains in `[a_t, b_t]`.

Choose a next corridor satisfying

```text
L_(t+1) < a_t
b_t < U_(t+1).
```

The strict gaps are deliberately left for finite sampling.

## Sampling escape bounds

Let `M_min` be the lower bound on simulator gene copies implied by `N_min` and
the effective-size closure.  Binary Chernoff/KL bounds give

```text
P(p_j,t+1 < L_(t+1))
<= exp[-M_min D(L_(t+1) || a_t)],

P(p_j,t+1 > U_(t+1))
<= exp[-M_min D(U_(t+1) || b_t)].
```

For `J` patches, add the two bounds and apply a patch union bound.  For a
finite path of corridors, add each one-step any-patch escape bound.  No
independence assumption is needed.

## What this resolves

This theorem supplies the correct replacement for a fixed-interval retention
premise in the H2 chain:

```text
frequency moving corridor
-> time-indexed H-alpha multiplier envelope
-> finite-horizon probabilistic genetic-lead bound.
```

## What it does not yet resolve

The moving corridor itself still requires declared interaction and census bounds.
The next theorem must combine the corridor endpoints with a **time-indexed**
H-alpha multiplier rather than the earlier fixed-interval multiplier.
