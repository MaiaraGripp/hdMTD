# hdMTD (development version)

# hdMTD 0.1.2

## New
* Accessor functions for "MTD": `pj()`, `p0()`, `lambdas()`,
  `lags()`, `Lambda()`, `states()`, and `transitP()`. See `?MTD-accessors`.
* Accessor functions for "MTDest": `pj()`, `p0()`, `lambdas()`,
  `lags()`, `S()` and `states()`. See `?MTD-accessors`.
* Accessor functions for "hdMTD": `S()` and `lags()`. See `?MTD-accessors`.
* Methods for "MTD" and "MTDest" objects: added `print()`, `summary()`, `coef()`, `logLik()`
and `probs()`. For compact inspection of lag sets, state space, mixture weights and more.
 See `?MTD-methods` and `?MTDest-methods`.
* Methods for "hdMTD" objects: added `print()` and `summary()` for compact inspection of
  lag selection results. See `?hdMTD-methods`.
* Coercion: new `as.MTD()` to rebuild an "MTD" object from an "MTDest" fit.

## Changes
* `probs()` is now a S3 generic with methods for "MTD" and "MTDest". Returns one-step-ahead predictive probabilities
  either for specific contexts (`context=`) or from sample rows (`newdata=`). If neither is supplied, it returns
  the full global transition matrix (`transitP(object)` for `MTD`; `transitP(as.MTD(object))` for `MTDest`).
* Renamed the sample-based estimator `probs(X, S, ...)` to `empirical_probs(X, S, ...)` to avoid ambiguity:
  `empirical_probs()` estimates transition probabilities from data, while `probs()` returns predictive probabilities
  from model/fit objects.

## Fixes
* Replaced `any(is.na(X))` with `anyNA(X)` in `checkSample()` for efficiency and clarity.

## Package cleanup
* Removed unused datasets (`raindata`, `sleepscoring`, `testChains`).
* Updated examples to use simulated data (via `perfectSample()`) instead of the removed `testChains` dataset.
* Internal helpers marked `@keywords internal` so they no longer appear in `help(package="hdMTD")`.

# hdMTD 0.1.1

* Relicensed the package from MIT to GPL-3.
* Removed an unintended `README.md` file from the package source.

