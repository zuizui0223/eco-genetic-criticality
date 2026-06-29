# Finite H1 hysteresis continuation audit

The canonical H1 reduction has an analytic bistable barrier interval and a
continuation check. The finite-bin multipatch simulator contains density,
realised-trait and allele feedback, finite recruitment, and optional migration;
therefore it needs a separate simulation experiment before any hysteresis claim
is made for that closure.

## Continuation rather than independent runs

For each applicable scenario and parameter cell, the audit derives the strict
canonical bistable interval and constructs one barrier grid extending beyond
both endpoints. It then performs two routes:

- **rising route:** starts from the canonical high stable interaction below the
  interval and increases the barrier;
- **falling route:** starts from the canonical low stable interaction above the
  interval and decreases the barrier.

At every barrier stage, the terminal population, interaction state, allele
frequency, and realised trait-bin abundances become the next stage's initial
state. The route is therefore a stateful continuation experiment, not a set of
independent fixed-barrier simulations.

## Predeclared finite hysteresis predicate

At barriers strictly inside the canonical bistable interval, compare terminal
mean interaction at the same barrier:

```text
mean_q(rising route) - mean_q(falling route).
```

A replicate supports finite hysteresis when the maximum internal difference is
larger than the declared `interaction_separation_threshold`. The stronger
`finite_h1_hysteresis_mechanism_supported` predicate additionally requires at
least one internal barrier where the rising route retains a potential high-trait
component but the falling route does not.

The audit also reports thresholded jump proxies:

- the first rising barrier at which mean interaction is at or below the declared
  low-state threshold;
- the first falling barrier at which mean interaction is at or above the
  declared high-state threshold.

Their difference is descriptive. It is not required to call route-memory
hysteresis, because finite trajectories can be noisy or can cross the arbitrary
state thresholds without a clean jump.

## Applicability and counterexamples

The experiment is run only where the local canonical map has a strict bistable
interval and the midpoint has a canonical branch-dependent high-trait mode.
Cells outside that region are stored as unavailable, not as finite failures.

Within applicable cells, all outcomes are retained:

- sustained route separation;
- route separation without a high-trait switch;
- high-trait switching without large interaction separation;
- convergence of both histories; and
- any apparent reversal of the canonical expectation.

These are model-specific Type S results. A finite route-memory result does not
prove bistability in an arbitrary ecosystem, and a failure does not invalidate
the canonical theorem. It identifies where additional finite mechanisms change
the canonical picture.

## Scope audit

Every stage includes the existing H1 theorem-boundary audit. This makes it
visible whether a trajectory lies in the exact canonical embedding or is a
controlled departure through density, trait feedback, allele feedback, migration,
or multiple patches.
