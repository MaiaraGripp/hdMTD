# hdMTD 0.1.2

## Bug fixes
* Internal helper function `checkSample()` is now hidden from `help(package = "hdMTD")`
  by marking it with `@keywords internal`.
* Replaced `any(is.na(X))` by `anyNA(X)` in `checkSample()` for efficiency and clarity.

## Package cleanup
* Removed unused datasets (`raindata`, `sleepscoring`, `testChains`), leaving only `tempdata` as the included dataset.
* Updated all examples to use simulated data (via `perfectSample()`) instead of the removed `testChains` dataset.


# hdMTD 0.1.1

* Relicensed the package from MIT to GPL-3.
* Removed an unintended `README.md` file from the package source.

