# H3 fragmentation-gradient sensitivity results

## Status

This document reports the first complete outcome of the preregistered post-review fragmentation-gradient sensitivity. It is **new supplementary finite Type S evidence** and does not modify the historical H1/H3 evidence ledger.

Authoritative workflow:

- run: `31937210601`;
- publication artifact: `9261157020` (`h3-fragmentation-gradient-publication`);
- artifact digest: `sha256:424031d0f6bcdf75c13e03deb35324f0d3f6fd46f58ff7b34961bbd00556537c`;
- fresh master seeds: `20260820`–`20260824`;
- patch counts: `1, 2, 3, 4, 6, 8, 12, 16`;
- attempted source replicates: 1,200;
- complete H1-prepared sources: 1,037;
- repeated-measures projection rows: 9,600;
- projection-supported outcomes at every patch count: 1,037;
- checksum manifest: five publication files, all verified.

The first execution attempt is not evidence: it stopped before an admissible complete gradient because of a positional-versus-keyword API drift in an existing private source-preparation helper. That compatibility call was fixed before any gradient outcome was inspected, targeted parent-chain/gradient tests passed, and the full fixed-design campaign was rerun from the same preregistered seeds and patch counts.

## Main result

The historical one-large versus four-isolated contrast was **not unique to the four-patch endpoint**. The dominant ecological transition occurred as soon as the same prepared high state was split from one patch into two equal isolated patches.

Across the 1,037 fresh H1-prepared source replicates, the pooled paired medians at two patches retained only:

- `0.001744` of the one-patch final interaction, a **99.83% reduction**;
- `0.221311` of the one-patch local effective size, a **77.87% reduction**;
- `0.282918` of the one-patch realised high-trait mass, a **71.71% reduction**.

All three quantities were below their paired one-patch value in **1,037/1,037** supported sources at two patches.

Thus the original four-patch result is not an artefact of selecting an unusually extreme endpoint. Under this declared finite closure, the response is instead strongly threshold-like at the first subdivision.

## Replication of the historical four-patch contrast

The fresh-seed gradient reproduces the magnitude of the historical four-patch result closely.

| metric | historical locked H3 median reduction | fresh gradient, n=4 median reduction |
|---|---:|---:|
| final interaction | 99.86% | 99.86% (`0.998553`) |
| local effective size | 88.73% | 88.73% (`0.887295`) |
| realised high-trait mass | 68.87% | 69.82% (`0.698200`) |

The historical values came from 1,055 H1-qualified paired replicates; the new values come from 1,037 independently prepared sources under the fresh gradient seed family. The gradient is therefore a new sensitivity analysis, not a re-expression of the historical trajectories.

## Different state variables have different gradient shapes

The gradient does **not** support a claim that every ecological quantity declines smoothly and monotonically with fragment number.

### Interaction

The cell-specific median retained interaction decreased with patch count in all 12 frozen primary cells. The pooled paired median was already `0.001744` at two patches and declined further to `0.001244` at 16 patches.

### Local effective size

The cell-specific median retained local effective size also decreased monotonically with patch count in all 12 cells. Pooled median retained fractions were:

| patches | retained local effective size | paired median reduction |
|---:|---:|---:|
| 1 | 1.0000 | 0.00% |
| 2 | 0.2213 | 77.87% |
| 3 | 0.1488 | 85.12% |
| 4 | 0.1127 | 88.73% |
| 6 | 0.0744 | 92.56% |
| 8 | 0.0579 | 94.21% |
| 12 | 0.0413 | 95.87% |
| 16 | 0.0331 | 96.69% |

### Realised high-trait mass

Realised high-trait mass behaved differently. It fell sharply at the first split but did not continue monotonically downward. The pooled retained fraction was `0.2829` at two patches, `0.3018` at four patches, and `0.3939` at 16 patches. In other words, it remained far below the one-patch reference but partially recovered as the landscape was subdivided further.

The minimum cell-specific median occurred at two, three, four, or eight patches depending on the primary cell. This non-monotonic occupancy response is important: it prevents the supplementary sensitivity from being described as a single smooth fragmentation dose-response.

## Potential viability and realised occupancy separate immediately after fragmentation

A particularly strong state-separation result emerged from the fresh gradient:

- at one patch, potential high-trait viability was present in **1,037/1,037** projection-supported outcomes;
- at every tested patch count from 2 through 16, potential high-trait viability was present in **0/1,037** outcomes;
- nevertheless, realised high-trait occupancy persisted at the 30-generation endpoint in approximately 99.6–100% of supported trajectories, depending on patch count.

This is not a contradiction. Potential viability and realised occupancy are separate states in the model. The carried full state can retain realised high-trait occupancy for a finite interval after the interaction environment no longer supports a potential high-trait component. The gradient therefore reinforces the manuscript's distinction among ecological state, potential viability, realised occupancy, and demographic/genetic state.

## Relation to the analytical K=4 boundary

For the canonical interaction map, `K = kappa*A/A_ref > 4` permits the three-fixed-point geometry, while `K <= 4` rules it out. The finite gradient did **not** collapse exactly when this analytical coordinate crossed four.

At two patches, all 12 primary cells still had `K_fragment > 4`; the cell range was `7.5–15.0`. Nevertheless, final potential high-trait viability was absent in every supported two-patch outcome and interaction had already collapsed strongly relative to the paired one-patch state.

Across the declared cells, the number of cells with `K_fragment <= 4` was:

| patches | cells with K <= 4 / 12 |
|---:|---:|
| 1 | 0 |
| 2 | 0 |
| 3 | 0 |
| 4 | 3 |
| 6 | 9 |
| 8 | 12 |
| 12 | 12 |
| 16 | 12 |

Therefore `K=4` remains an analytical boundary for the specified one-state map, not a fitted finite fragmentation threshold. The finite projected system also depends on barrier position, source state, stochastic finite dynamics, and the declared trait/demographic closure.

## Full paired-gradient summary

| isolated patches | interaction retained | effective-size retained | high-trait mass retained | all three below paired n=1 |
|---:|---:|---:|---:|---:|
| 1 | 1.000000 | 1.000000 | 1.000000 | 0/1,037 |
| 2 | 0.001744 | 0.221311 | 0.282918 | 1,037/1,037 |
| 3 | 0.001539 | 0.148760 | 0.284800 | 1,037/1,037 |
| 4 | 0.001447 | 0.112705 | 0.301800 | 1,037/1,037 |
| 6 | 0.001339 | 0.074380 | 0.321053 | 1,037/1,037 |
| 8 | 0.001303 | 0.057851 | 0.321429 | 1,036/1,037 |
| 12 | 0.001260 | 0.041322 | 0.343290 | 1,037/1,037 |
| 16 | 0.001244 | 0.033058 | 0.393880 | 1,037/1,037 |

The single exception to the three-metric directional pattern at eight patches arose because one paired source had realised high-trait mass slightly above its own one-patch value while interaction and effective size were both much lower. It is retained rather than removed.

## Manuscript interpretation

The supplementary result supports a narrower and stronger statement than the original endpoint-only wording:

> The large one-versus-four fragmentation contrast was reproduced with fresh seeds and was already present after the first subdivision into two isolated patches. Interaction and local effective size then declined further with patch count, whereas realised high-trait mass showed a sharp initial loss followed by partial finite-state recovery.

It does **not** support:

- a universal smooth dose-response to fragment number;
- a claim that `K=4` is the observed finite collapse point;
- a claim that every functional state decreases monotonically with fragmentation;
- replacement of the historical H3 evidence ledger with this post-review campaign.

## Files

The publication artifact contains:

- `h3_fragmentation_gradient_records.csv` — 9,600 repeated-measures rows;
- `h3_fragmentation_gradient_cell_summary.csv` — 96 cell × patch-count summaries;
- `h3_fragmentation_gradient_pooled_summary.csv` — eight pooled patch-count summaries;
- `h3_fragmentation_gradient_metadata.json`;
- `figure_s_fragmentation_gradient.svg`;
- `MANIFEST.sha256`.
