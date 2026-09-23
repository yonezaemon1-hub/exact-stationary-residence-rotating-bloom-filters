# Exact Stationary Residence-Time Laws for Occupancy-Triggered Rotating Bloom-Filter Rings

Author: Ryutaro Yonezu

Release candidate: `v1.0.0`

## Main result

For the stationary residence count `R_s` of a tagged insertion in an
`s`-filter occupancy-triggered rotating Bloom-filter ring,

`G_Rs(z) = H(z)/E[L] * G_L(z)^(s-1)`.

The package includes the reviewed V0.4 manuscript, exact standard-library
Python implementation, Windows real-machine evidence, claim audit, and
tail-asymptotic audits.

## Verification status

- exact kernel self-tests: PASS
- Windows real-machine numerical run: PASS
- 100,000-sample Monte Carlo check: PASS
- original tail-constant audit: PASS
- strengthened finite-n tail-ratio audit: PASS
- independent external review reported by the author: PASS
- publication claim boundary: occupancy hitting time and recycling mechanisms
  are treated as prior art; the candidate contribution is the stationary
  tagged-insertion residence law and its consequences.

## Frozen manuscript

`NO20_PAPER_V0_4.pdf` is the reviewed manuscript and is intentionally unchanged
by the strengthened audit.

Paper licensing can be declared separately on the archival record. Source code
in this package is MIT licensed.

## Planned repository

`yonezaemon1-hub/exact-stationary-residence-rotating-bloom-filters`

DOI fields are intentionally omitted until archival publication.
