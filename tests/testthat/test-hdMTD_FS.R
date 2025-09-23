test_that("hdMTD_FS function works as expected", {
  # Create a sample MTD chain
  set.seed(1)
  X <- perfectSample(MTDmodel(c(2, 4), c(0, 1), lam0 = 0.05), 2000)
  d <- 5
  l <- 2

  # Test with the default parameters
  result_1 <- hdMTD_FS(X = X, d = d, l = l)
  expect_true(is.numeric(result_1))
  expect_true(all(result_1 %in% 1:d))

  # Test with custom parameters and elbowTest
  result_2 <- hdMTD_FS(X = X, d = d, l = l, elbowTest = TRUE)
  expect_true(is.numeric(result_2))
  expect_true(all(result_2 %in% 1:d))
  expect_lt(length(result_2), l+0.01)

  # Test with warn enabled
  expect_warning(hdMTD_FS(X = X, d = d, l = 3, warn = TRUE))
})
