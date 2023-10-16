test_that("hdMTD works correctly with different methods", {
  # Create a sample MTD chain
  set.seed(1)
  X <- perfectSample(MTDmodel(c(1, 5), c(0, 1)), 5000)
  d <- 5

  # Test the "FS" method
  result_fs <- hdMTD(X = X, d = d, method = "FS", l = 2)
  expect_type(result_fs, "integer")
  expect_equal(length(result_fs), 2)
  expect_true(all(result_fs%in%1:d))
  # Add more specific checks based on your expectations for the "FS" method

  # Test the "FSC" method
  result_fsc <- hdMTD(X = X, d = d, method = "FSC", alpha = 0.001, xi = 1, l = 3)
  expect_type(result_fsc, "integer")
  expect_lt(length(result_fsc), 3.1)
  expect_true(all(result_fsc%in%1:d))
  # Add more specific checks based on your expectations for the "FSC" method

  # Test the "BIC" method
  result_bic <- hdMTD(X = X, d = d, method = "BIC", xi = 1, minl = 3, maxl = 3)
  expect_type(result_bic, "double")
  expect_lt(length(result_bic), 3.1)
  expect_true(all(result_bic%in%1:d))
  # Add more specific checks based on your expectations for the "BIC" method

  # Test the "CUT" method
  result_cut <- hdMTD(X = X, d = d, method = "CUT", S = c(1, 5, 3, 2), alpha = 0.001, xi = 1)
  expect_type(result_cut, "double")
  expect_lt(length(result_cut), 5.1)
  expect_true(all(result_cut%in%1:d))
  # Add more specific checks based on your expectations for the "CUT" method
})
