## Resubmission (hdMTD 0.1.2)

### Summary of changes
- Added S3 class for EM fits: `MTDest` with `print()`, `summary()`, `coef()`, and `logLik()`.
- Added public accessor functions for `MTD`/`MTDest`: `pj()`, `p0()`, `lambdas()`, `lags()` (ℤ⁻),
  `Lambda()`/`S()` (ℕ⁺), `states()`, and `transitP()` (for `MTD`).
- Added `MTD` methods: `print()`, `summary()`, `coef()`.
- Introduced `as.MTD()` to coerce an `MTDest` object to an `MTD` object.
- Replaced `any(is.na())` with `anyNA()` in `checkSample()`. 
- Marked internal helpers with `@keywords internal` so they no longer appear in `help(package="hdMTD")`.
- Removed unused datasets and updated examples to use simulated data.

### R CMD check results (current)
0 errors | 0 warnings | 0 notes

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
