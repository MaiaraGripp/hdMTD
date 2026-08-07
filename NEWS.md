# hdMTD 0.1.5

## Documentation
* Added computational-cost warnings to the documentation of `hdMTD_BIC()` and `hdMTD_CUT()`, emphasizing that unrestricted candidate sets may substantially increase running time and memory requirements.
* Recommended using a reduced candidate lag set `S`, when possible, for example based on the output of `hdMTD_FS()`.
* Clarified in the documentation of functions that accept chain samples that missing values (`NA`) are not allowed.

## Vignettes

- Added an official package vignette with replication materials for the examples and analyses presented in Section 5 of the article describing the `hdMTD` package.
- Included pre-computed results for computationally intensive analyses so that the vignette can be built without rerunning long computations.

## Tests
* Added an automated test confirming that inference functions reject samples containing missing values with a clear error message.

# hdMTD 0.1.4

## Fixes and improvements
* Simplified the S3 class hierarchy: fitted objects of class `MTDest` now inherit from `MTD`.
* Removed redundant class checks from S3 methods, relying on method dispatch.
* Eliminated duplicated S3 methods (`coef.MTDest`, `probs.MTDest`), which now work via inheritance.
* Streamlined plotting code: overlapping plot types in `plot.MTDest` are delegated to `plot.MTD` using `NextMethod()`.
* Simplified `summary()` methods by printing directly from the summary output, removing auxiliary summary-print classes.
* Reorganized internal helper functions (e.g., moved `PI()` and `sx()` from `utils.R` to dedicated files).
* Improved and consolidated user-facing documentation, explicitly listing available S3 methods and accessors in the constructors’ help pages.

# hdMTD 0.1.3
  
## Fixes
* Modified the tie-breaking rule in the FS (Forward Selection) procedure to ensure deterministic behavior.
* Updated `MTD-methods`, `MTDest-methods` and `MTD-accessors` documentation to remove redundant links in the help system and streamline method listings.
* Sample size is now a required argument in `perfectSample()`.
* Improved the error message in `logLik.MTD()` when a sample is not provided.

## New Features
* Added a `plot.MTD()` method for visualizing MTD models, including bar plots of lag contributions and mixture weights, as well as directed weighted graphs (via igraph) representing each lag-specific transition matrix.
* Added a `plot.MTDest()` method for fitted `MTDest` objects, which mirrors `plot.MTD()` but also includes EM iteration diagnostics (log-likelihood variation per update) when available.

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

