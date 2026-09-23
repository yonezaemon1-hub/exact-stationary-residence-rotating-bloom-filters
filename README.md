# Exact Stationary Residence-Time Laws for Occupancy-Triggered Rotating Bloom-Filter Rings

Author: Ryutaro Yonezu

Release: `v1.0.0`

Software / Evidence DOI: `10.5281/zenodo.22927019`

Zenodo record: https://zenodo.org/records/22927019

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

Frozen manuscript SHA-256:

`3459e6caac9ca3f6a1fdc8c0b1dee8e0c51683db41e83d67bdf23d038338088e`

Release package SHA-256:

`857f7423a9a2c56fd37514bfb47935ae8efac7fbd45e0c76a05083e45e755b54`

Source code in this repository is MIT licensed.
