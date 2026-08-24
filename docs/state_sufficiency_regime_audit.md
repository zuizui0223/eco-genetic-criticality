# State sufficiency and regime definition audit

## Question

Can two systems produced by different fragmentation histories be treated as the same functional-fragmentation regime once their present state is matched?

## Result 1 — exact full-state sufficiency under the declared simulator

The declared simulator is Markov in the explicit finite state

`(population, interaction, high-allele frequency, realised trait-bin distribution/abundance)`

together with the fixed future parameter schedule and random-number seed. The update at generation `t+1` uses only these current quantities and the declared parameters; no history label or earlier trajectory value enters separately.

Therefore, if two histories arrive at the same complete simulator state and subsequently receive the same forcing and stochastic law, their future trajectory distributions are identical under this model closure. With the same RNG seed, the realised trajectories are identical. The regression test `test_exact_present_state_is_future_sufficient_under_declared_closure` locks this property.

This is a model theorem, not a universal ecological theorem. If a future model adds latent memory variables (soil legacies, age structure, epigenetic state, unmeasured partner identity, propagule banks, learned behaviour, etc.), those variables become part of the full state required for sufficiency.

## Result 2 — common aggregate summaries are not sufficient

A two-patch counterexample holds fixed all of the following at generation 0:

- total and patch census sizes;
- the marginal distribution and population-weighted mean of interaction state `q`;
- the marginal distribution and population-weighted mean of high-associated allele frequency `p`;
- `H_alpha`, `H_gamma`, and `F_ST`;
- realised trait-bin state.

The only change is the joint spatial alignment of `q` and `p`:

- aligned: `q=(0.8,0.2)`, `p=(0.8,0.2)`;
- anti-aligned: `q=(0.8,0.2)`, `p=(0.2,0.8)`.

Because the interaction update is patchwise and its support signal combines local interaction and local allele state, the next interaction field differs between the two systems even though the usual aggregate summaries and both marginals are identical. The regression test `test_coarse_aggregate_state_is_not_future_sufficient_when_patch_alignment_differs` locks this counterexample.

## Consequence for functional-fragmentation regimes

A defensible operational regime cannot be defined only by habitat amount, occupancy, mean interaction support, mean genetic diversity, or separate marginal distributions. For this model, a sufficient regime representation must retain the joint spatial state—or an empirically justified lower-dimensional statistic proven to be sufficient for the target future functional-loss process.

This sharpens the urban–island convergence question:

> Different fragmentation mechanisms may converge if they generate the same future-relevant joint functional state, not merely the same coarse averages or category labels.

The empirical burden is therefore to test conditional predictive equivalence: after matching the candidate regime variables, does urban-versus-island origin add predictive information about subsequent realised functional loss? If origin still predicts the future, the proposed regime omitted a state variable or memory mechanism.

## Scope boundary

This audit does not claim that real ecological systems are Markov at the measured scale, nor that urban and island systems currently converge. It establishes (i) exact state-history equivalence inside the declared simulator, and (ii) a constructive failure of coarse aggregate equivalence even inside that simple closure.
