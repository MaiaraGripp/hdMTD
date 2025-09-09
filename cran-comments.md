## Resubmission (hdMTD 0.1.2)

- `MTDest()`: now returns an S3 object ("MTDest") with iteration diagnostics (optional) and logLik.
  - Solved inconsistencies:
    a) When M was used to stop updates, the returned vector of delta log-likelihoods 
       (now called deltaLogLik, formerly distLogL) had size = (number_of_updates + 1), 
       because it included a last element smaller than M. Now, we've introduced the 
       variable `lastComputedDelta`, ensuring that length(deltaLogLik) = number_of_updates, 
       and we return `lastComputedDelta` separately. If `lastComputedDelta` differs from 
       the last element in `deltaLogLik`, it means `lastComputedDelta < M` and the convergence 
       criterion was met.
    b) When the max number of iterations was reached (i.e., the `nIter` criterion was used), 
       `oscillations` was being calculated using a model that was not updated with the final 
       estimated parameters (fixed).
    c) Improved performance: moved `indexa` out of the loop; avoided naming variables inside loops; 
       pre-allocated memory for `deltaLogLik`.

- Internal helper function `checkSample()` is now hidden from `help(package = "hdMTD")` 
by marking it with `@keywords internal`.
- Replaced `any(is.na(X))` with `anyNA(X)` in `checkSample()`.
- Removed unused datasets (`raindata`, `sleepscoring`, `testChains`), which were not required for package functionality.
- Updated examples to generate data using `perfectSample()` instead of relying on `testChains`.
- Changes documented at `NEWS.md`.

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
