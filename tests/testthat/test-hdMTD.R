test_that("hdMTD works correctly with different methods", {
  # Create a sample MTD chain
  set.seed(1)
  X <- perfectSample(MTDmodel(c(1, 5), c(0, 1)), 5000)
  d <- 5

  ## helpers
  expect_hlagset <- function(obj, d, max_len = Inf) {
    expect_s3_class(obj, "hdMTD")

    S_hat   <- S(obj)

    expect_true(all(S_hat %in% 1:d))
    expect_lte(length(S_hat), max_len)

  }

  # Test the "FS" method
  fit_fs <- hdMTD(X = X, d = d, method = "FS", l = 2)
  expect_hlagset(fit_fs, d, max_len = 2)

  result_fs <- hdMTD(X = X, d = d, method = "FS", l = 2)
  expect_type(result_fs, "integer")
  expect_equal(length(result_fs), 2)
  expect_true(all(result_fs%in%1:d))
  # Add more specific checks based on your expectations for the "FS" method

  # Test the "FSC" method
  fit_fsc <- hdMTD(X = X, d = d, method = "FSC", alpha = 0.001, xi = 1, l = 3)
  expect_hlagset(fit_fsc, d, max_len = 3)
  # Add more specific checks based on your expectations for the "FSC" method

  # Test the "BIC" method
  fit_bic <- hdMTD(X = X, d = d, method = "BIC", xi = 1, minl = 3, maxl = 3)
  expect_hlagset(fit_bic, d, max_len = 3)
  expect_equal(length(S(fit_bic)), 3)

  fit_bic_extras <- hdMTD(X = X, d = d, method = "BIC",
                          xi = 1, minl = 3, maxl = 3,
                          byl = TRUE, BICvalue = TRUE)
  ex <- attr(fit_bic_extras, "extras")
  expect_true(is.list(ex))
  expect_true("BIC_out" %in% names(ex))
  expect_true(is.numeric(ex$BIC_out))
  expect_true(length(ex$BIC_out) >= 1)

  # Test the "CUT" method
  fit_cut <- hdMTD(X = X, d = d, method = "CUT",
                   S = c(1, 5, 3, 2), alpha = 0.001, xi = 1)
  expect_hlagset(fit_cut, d, max_len = 4)
  expect_true(all(S(fit_cut) %in% c(1, 2, 3, 5)))
})
