# hdMTD (development version)

## Bug fixes
* Internal helper function `checkSample()` is now hidden from `help(package = "hdMTD")`
  by marking it with `@keywords internal`.
* Replaced `any(is.na(X))` by `anyNA(X)` in `checkSample()` for efficiency and clarity.


# hdMTD 0.1.1

* Relicensed the package from MIT to GPL-3.
* Removed an unintended `README.md` file from the package source.

