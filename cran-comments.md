## Resubmission (hdMTD 0.1.3)

### Summary of changes

- Added `plot.MTD()` method for objects of class "MTD", allowing graphical visualization of lag weights and contributions.
- Added `plot.MTDest()` method for "MTDest" objects, mirroring `plot.MTD()` and including EM iteration diagnostics (when available).
- Modified the tie-breaking rule in the Forward Selection (FS) procedure to ensure deterministic behavior.
- Updated documentation for `MTD-methods`, `MTDest-methods`, and `MTD-accessors` to remove redundant links and streamline method listings.
- `perfectSample()` now requires explicit sample size as an argument.
- Improved the error message in `logLik.MTD()` when a sample is not provided.

### R CMD check results (current)
0 errors | 0 warnings | 0 notes

---

## Previous submission notes (for reference)

version - hdMTD 0.1.2

### Summary of changes
- Added S3 class for EM fits: `MTDest` with `print()`, `summary()`, `coef()`, `logLik()` and `probs()`.
- Added S3 class for lag-selection fits: `hdMTD` with `print()` and `summary()`.
- Added public accessor functions for `MTD`/`MTDest`: `pj()`, `p0()`, `lambdas()`, `lags()` (ℤ⁻),
  `Lambda()`/`S()` (ℕ⁺), `states()`, and `transitP()` (for `MTD`). `S()` and `lags()` also apply
to `hdMTD`.
- Introduced `as.MTD()` to coerce an `MTDest` object to an `MTD` object.
- Clarified terminology: sample-based probs() has been renamed to empirical_probs(), while probs() is now an
 S3 generic for model objects.
- Added `MTD` methods: `print()`, `summary()`, `coef()`, `logLik()` and `probs()`.
- Replaced `any(is.na())` with `anyNA()` in `checkSample()`. 
- Cleaned package: removed unused datasets, marked internal helpers as @keywords internal, and updated examples
 to use simulated data.

---

## Previous submission notes (for reference)

version - hdMTD 0.1.1

- Relicensed the package from MIT to GPL-3.
- Removed an unintended `README.md` from the package source.
- Added `NEWS.md` documenting these changes.
- No code/API changes.

### R CMD check results (current)
0 errors | 0 warnings | 0 notes

version - hdMTD 0.1.0

R CMD check results
0 errors | 0 warnings | 1 note

## Acronym explanation

The following acronyms are used in the DESCRIPTION file:

- **MTD**: Mixture Transition Distribution, a class of models introduced by Raftery (1985)
for representing high-order Markov chains with fewer parameters.
- **hdMTD**: Name of the package, short for *high-dimensional Mixture Transition Distribution*.

## Proper nouns and author names flagged

In a previous submission, the following were flagged as possibly misspelled, but are
valid proper nouns:

- **MTD** (see above)
- **Ost** and **Takahashi**: authors of the reference included in the Description.
- **hdMTD**: the name of the package.

## Explanation of methods

- **Forward Stepwise and Cut**: This is a two-stage algorithm for lag selection in
  Mixture Transition Distribution models. The "Forward Stepwise" part incrementally adds
  relevant lags, and the "Cut" step prunes unnecessary ones. The method is described in 
  the referenced article by Ost & Takahashi (2023).

## URLs with status 403

In a previous submission, two URLs from the Australian Bureau of Meteorology (BOM) were
flagged with status 403 (Forbidden). These URLs are no longer referenced directly in the
documentation.

Instead, BOM is cited as the original data source, and a Kaggle-hosted version of the
dataset (where the data was actually acquired) is provided. A note has also been added
to clarify that the original BOM URL may require manual browser access.
