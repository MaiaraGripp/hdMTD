test_that("MTDmodel function works as expected", {
  # Test with default parameters
  result_1 <- MTDmodel(Lambda = c(1, 3), A = c(4, 8, 12))
  expect_true(is.list(result_1))
  expect_true(all(names(result_1) %in% c("P", "lambdas", "pj", "p0", "Lambda", "A", "single_matrix", "call")))

  expect_error(MTDmodel(Lambda = c(1, 3), A = c(0, 1), lam0 = 0.01, lamj = c(0.99)))
  expect_error(MTDmodel(Lambda = c(1, 3), A = c(0, 1), lam0 = 0.01, lamj = c(0.7, 0.3)))
  expect_error(MTDmodel(Lambda = c(1, 3), A = c(0, 1), lam0 = 0.01, lamj = c(0.69, 0.3), p0 = 1))

  expect_equal(probs(result_1), transitP(result_1))

})
