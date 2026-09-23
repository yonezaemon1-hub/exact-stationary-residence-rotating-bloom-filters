# No.20 Claim Audit V0.4

Date: 2026-09-23

## Safe paper claim

The paper derives a stationary tagged-insertion residence-count law for an occupancy-triggered rotating Bloom-filter ring under an iid colliding-hash placement model. The contribution is the tagged residence law and its consequences, not occupancy hitting times, Bloom-filter recycling, two-phase buffering, or average false-positive/false-negative analysis.

## Closest prior work

- O'Neill (2022): negative occupancy / coupon-collector hitting-time distributions, convolution representations, generating functions, and moments.
- Dozier, Salamatian, Rubenstein (INFOCOM 2024): occupancy-triggered recycling Bloom filters, Markov/renewal analysis, expected cycle capacity, average false-positive rate, and two-phase variants.
- Dozier, Salamatian, Rubenstein (POMACS 2024): false-negative analysis for one- and two-phase recycling Bloom filters using Markov and renewal models.
- Yoon (2010) and Age-Partitioned Bloom Filters (2020): multi-buffer aging / sliding-window retention mechanisms.

## Mechanism convention audit

The sigma-bounded INFOCOM 2024 model applies the k hashes, sets the corresponding bits, and then resets when the post-insertion occupancy exceeds sigma. Its two-phase variant leaves the full active filter intact, clears the old frozen filter, and switches roles. This matches the paper's post-insertion threshold convention for a tagged item in the newly frozen filter.

## Claims to avoid

Do not claim:
- the occupancy threshold hitting distribution is new;
- occupancy-triggered recycling is new;
- two-phase / multi-buffer retention is new;
- this is the first analysis of recycling Bloom filters;
- false negatives caused by recycling are new.

## Residual statement

Targeted searches did not locate the same arbitrary-s stationary tagged-insertion residence PGF

G_Rs(z) = H(z)/E[L] * G_L(z)^(s-1)

together with the exact tagged residence moments, polynomial-geometric tail, and Poisson wall-clock transform for occupancy-triggered rings.

This is a targeted-search result, not a proof of global novelty.
