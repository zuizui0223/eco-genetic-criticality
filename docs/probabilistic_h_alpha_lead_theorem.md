# Probabilistic H-alpha lead theorem for finite trait recruitment

## Why L3 is not directly pathwise for finite-bin recruitment

L3 combines a deterministic upper bound on diversity decay with a deterministic
lower bound on realised high-trait abundance. That is appropriate for a
deterministic trait recursion.

The finite-bin model instead samples realised trait abundance. Even when the
high-trait recruitment probability is positive, a finite cohort can contain no
high-trait recruits with nonzero probability. Therefore a nontrivial pathwise
lower bound on realised high-trait abundance cannot generally hold for every
random realization.

The correct finite-population replacement is a high-probability theorem.

## Assumptions

Fix a time `t`. Assume:

```text
E[H_t] <= h_0 lambda_bar^t,      0 < lambda_bar < 1.
```

This is an expected H-alpha decay bound.

For every generation `s=1,...,t`, assume the realised high-trait abundance
conditionally stochastically dominates a Binomial variable with:

```text
cohort size >= n_min
high-trait recruitment probability >= pi_min.
```

Let `a` be the declared realised-occupancy threshold and define

```text
mu_min = n_min pi_min.
```

When `a < mu_min`, Chernoff's lower-tail inequality gives the one-generation
risk bound

```text
P(N_H,s <= a) <= exp[-mu_min (1-a/mu_min)^2 / 2] = epsilon_bin.
```

By a union bound, realised trait loss at or before `t` has probability at most

```text
epsilon_trait(t) <= min(1, t epsilon_bin).
```

## Theorem L4

Let `h_warn` be the predeclared H-alpha warning threshold. Markov's inequality
applied to the expected diversity bound gives

```text
P(H_t > h_warn) <= h_0 lambda_bar^t / h_warn.
```

Since `H_t <= h_warn` implies `tau_H <= t`, and trait persistence through time
`t` implies `t < tau_trait_realised`, the union bound gives:

```text
P(tau_H <= t < tau_trait_realised)
>= max(0,
       1 - h_0 lambda_bar^t / h_warn - epsilon_trait(t)).
```

This is a sufficient lower probability bound for a genetic lead at a declared
time. It is deliberately conservative.

## What it proves

Under the declared expected-diversity and finite-recruitment assumptions, it
proves a lower bound on the probability of a lead event.

It does not prove:

```text
- that every realization has a genetic lead;
- that the simulator automatically satisfies n_min, pi_min, or lambda_bar;
- that the bound is sharp;
- that a lead occurs in an empirical ecosystem.
```

## Relation to simulation

The phase-boundary pilot identifies parameter cells where leads are observed.
It can suggest candidate regions in which a researcher should attempt to prove
or externally justify the required `lambda_bar`, `n_min`, and `pi_min` bounds.
It cannot replace those bounds.
