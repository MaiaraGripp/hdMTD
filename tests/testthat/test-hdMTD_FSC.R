test_that("hdMTD_FSC function works as expected", {
  # Create a sample MTD chain
  set.seed(1)
  X <- perfectSample(MTDmodel(c(1, 3), c(0, 1)), 2000)
  d <- 4
  l <- 3
  alpha <- 0.02

  # Test with default parameters
  result_1 <- hdMTD_FSC(X = X, d = d, l = l, alpha = alpha)
  expect_true(is.numeric(result_1))
  expect_true(all(result_1 %in% 1:d))

  # Test with custom parameters
  result_2 <- hdMTD_FSC(X = X, d = d, l = 2, alpha = 0.001)
  expect_true(is.numeric(result_2))
  expect_true(all(result_2 %in% 1:d))
})
