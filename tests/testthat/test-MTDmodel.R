test_that("MTDmodel function works as expected", {
  # Test with default parameters
  result_1 <- MTDmodel(Lambda = c(1, 3), A = c(4, 8, 12))
  expect_true(is.list(result_1))
  expect_true(all(names(result_1) %in% c("P", "lambdas", "p_j", "p0", "Lambda", "A")))

  # Test with custom parameters
  result_2 <- MTDmodel(Lambda = c(2, 4, 9), A = c(0, 1), w0 = 0.05, w_j = c(0.35, 0.2, 0.4),
                       p_j = list(matrix(c(0.5, 0.7, 0.5, 0.3), ncol = 2)), p0 = c(0.2, 0.8), single_matrix = TRUE)
  expect_true(is.list(result_2))
  expect_true(all(names(result_2) %in% c("P", "lambdas", "p_j", "p0", "Lambda", "A")))

  # Test with indep_part = FALSE
  result_3 <- MTDmodel(Lambda = c(2, 4, 9), A = c(0, 1), w0 = 0,
                       p_j = list(matrix(c(0.5, 0.7, 0.5, 0.3), ncol = 2)), single_matrix = TRUE, indep_part = FALSE)
  expect_true(is.list(result_3))
  expect_true(all(names(result_3) %in% c("P", "lambdas", "p_j", "p0", "Lambda", "A")))
})
