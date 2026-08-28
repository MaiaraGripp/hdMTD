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

  # Test cut_fraction
  cut_fraction <- 0.3
  nCUT <- ceiling(cut_fraction * length(X))

  X_CUT <- X[seq_len(nCUT)]
  X_FS <- X[(nCUT + 1):length(X)]

  S_manual <- hdMTD_FS(X_FS, d = 10, l = l)

  expected <- hdMTD_CUT(X_CUT, d = 10, S = S_manual)

  result <- hdMTD_FSC(X, d = 10, l = l, cut_fraction = cut_fraction)

  expect_equal(result, expected)
})
