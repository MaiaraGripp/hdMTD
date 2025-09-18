# hdMTD 0.1.2

## New
* Accessor functions for "MTD" and "MTDest": `pj()`, `p0()`, `lambdas()`,
  `lags()`, `Lambda()`/`S()`, `states()`, and `transitP()` (for `MTD`). These
  expose model components without relying on internal list structure. See `?MTD-accessors`.
* Methods for "MTD" and "MTDest" objects: added `print()`, `summary()`, and `coef()` for compact inspection of
  lag sets, state space, mixture weights, and a preview of the global transition matrix `P`. See `?MTD-methods` and `?MTDest-methods`.
* Coercion: new `as.MTD()` to rebuild an "MTD" object from an "MTDest" fit.
* `logLik()` method for "MTDest": allows direct extraction of the fitted log-likelihood.
* `probs()` S3 generic with methods for "MTD" and "MTDest". Returns one-step-ahead predictive probabilities
  either for specific contexts (`context=`) or from sample rows (`newdata=`). If neither is supplied, it returns
  the full global transition matrix (`transitP(object)` for `MTD`; `transitP(as.MTD(object))` for `MTDest`).

## Changes
* Renamed the sample-based estimator `probs(X, S, ...)` to `empirical_probs(X, S, ...)` to avoid ambiguity:
  `empirical_probs()` estimates transition probabilities from data, while `probs()` returns predictive probabilities
  from model/fit objects.

## Improvements
* `MTDest()` now returns an S3 object ("MTDest") with `print()`, `summary()`, `coef()`, and `logLik()` methods.
  This lets users inspect fits consistently and extract log-likelihoods.
* `hdMTD()` now returns an S3 object ("hdMTD") with `print()` and `summary()` methods for consistent inspection
  of high-dimensional lag selection results.

## Bug fixes
* Replaced `any(is.na(X))` with `anyNA(X)` in `checkSample()` for efficiency and clarity.

## Package cleanup
* Removed unused datasets (`raindata`, `sleepscoring`, `testChains`).
* Updated examples to use simulated data (via `perfectSample()`) instead of the removed `testChains` dataset.
* Internal helpers marked `@keywords internal` so they no longer appear in `help(package="hdMTD")`.

# hdMTD 0.1.1

* Relicensed the package from MIT to GPL-3.
* Removed an unintended `README.md` file from the package source.

