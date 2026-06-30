# Mutation-conditioned full-state H1–H2–H3 primary chain

## Preconditions fixed before this campaign

The chain is not an unconstrained scan. It imports the 12-cell primary domain
frozen from the independent mutation-H1 validation (Actions run `28436777080`).
A cell entered that domain only when every one of five independent seed blocks
had at least 0.75 same-replicate support for:

```text
finite H1 full-state hold
AND polymorphic high branch
AND H-alpha/H-gamma above the baseline warning thresholds
```

The remaining 15 validation cells are not discarded from the scientific record;
they are outside this campaign's primary dynamic domain.

## State path

For every new master-seed and replicate, the chain uses:

```text
mutation-conditioned nested H1 calibration
-> a finest-grid interior anchor
-> rising high-route replay
-> mutation-conditioned full-state high hold
-> conservation-preserving projection
-> one-large / equal-isolated / equal-migrating trajectories
```

The carried state contains population, interaction, high-allele frequency, and
realised trait-bin abundance. The projection preserves total population, every
trait-bin total, area, population-weighted allele frequency, and intensive
interaction. Therefore the fragmented landscapes do not receive invented
individuals or trait material.

## H2 interpretation

H2 is evaluated on the `equal_isolated` trajectory. H-alpha and H-gamma warning
first-passage times are compared with realised trait-loss first passage only when
both events occur and the projected generation-0 state has both heterozygosities
above their warning thresholds. Any missing event is retained as censored.

A generation-zero threshold crossing is baseline-ineligible, not an early
warning lead.

## H3 interpretation

The primary H3 contrast is:

```text
one_large versus equal_isolated
```

The finite fragmentation-pattern predicate requires lower terminal interaction,
lower mean local effective size, and lower realised high-trait mass in isolated
fragments. H-alpha/H-gamma differences are reported but are not silently folded
into that ecological predicate.

`equal_migrating` is a separate allele-frequency-mixing comparison. Its F_ST
modulation is reported without calling migration ecological or demographic
rescue.

## Fixed first run

The `Mutation Primary H1 H2 H3 Chain` workflow has no runtime inputs. It uses:

```text
chain master seeds: 20260810–20260814
replicates: 20 per primary cell per seed
primary cells: 12
H1 grids: 25, 49, 97
H1 stage / hold: 30 / 30
mutation closure: inherited from each frozen domain cell
```

This gives 1,200 H1-conditioned primary records before each record is projected
into three landscapes. The artifact manifest embeds the full frozen domain
ledger and the workflow writes a separate protocol identity file.

## Scope

This is Type S evidence for the declared finite symmetric-mutation closure. A
successful chain does not prove universal warning ordering, general bifurcation,
or a universal fragmentation law.
