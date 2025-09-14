# hdMTD 0.1.2

## New
* **Accessor functions** for `MTD` and `MTDest`: `pj()`, `p0()`, `lambdas()`,
  `lags()`, `Lambda()`/`S()`, `states()`, and `transitP()` (for `MTD`).
  These expose model components without relying on internal list structure. See `?MTD-accessors`. 
* **Methods for `MTD` and `MTDest` objects**: added `print()`, `summary()`, and `coef()` for compact inspection of
  lag sets, state space, mixture weights, and a preview of the global transition matrix `P`. See `?MTD-methods` and `MTDest-methods`.
* **Coercion**: new `as.MTD()` to rebuild an `MTD` object from an `MTDest` fit.  

## Improvements
* **`MTDest()` now returns an S3 object** ("MTDest") with `print()`, `summary()`, `coef()`,
  and `logLik()` methods. This lets users inspect fits consistently and extract log-likelihoods. 

## Bug fixes
* Replaced `any(is.na(X))` with `anyNA(X)` in `checkSample()` for efficiency and clarity.

## Package cleanup
* Removed unused datasets (`raindata`, `sleepscoring`, `testChains`).
* Updated examples to use simulated data (via `perfectSample()`) instead of the removed `testChains` dataset.
* Internal helpers marked `@keywords internal` so they no longer appear in `help(package="hdMTD")`.

# hdMTD 0.1.1

* Relicensed the package from MIT to GPL-3.
* Removed an unintended `README.md` file from the package source.

