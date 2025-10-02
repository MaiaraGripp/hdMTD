test_that("hdMTD_CUT function works as expected", {
  # Create a sample MTD chain
  set.seed(1)
  X <- perfectSample(MTDmodel(c(1, 4), c(0, 1)), 1000)
  d <- 4

  # Test with the default parameters
  result_1 <- hdMTD_CUT(X = X, d = d)
  expect_true(is.numeric(result_1))
  expect_true(all(result_1 %in% 1:d))

  # Test with custom parameters
  result_2 <- hdMTD_CUT(X = X, d = d, S = c(1, 4), alpha = 0.02, mu = 1, xi = 0.4)
  expect_true(is.numeric(result_2))
  expect_true(all(result_2 %in% c(1, 4)))

  # Test with the default parameters
  result_3 <- hdMTD_CUT(X = X, d = d, alpha = 1)
  expect_true(is.numeric(result_1))
  expect_true(all(result_1 %in% 1:d))

  # Test with warn enabled
  expect_warning(hdMTD_CUT(X = X, d = d, S = c(1, 4), alpha = 0.0065, warn = TRUE))
})
