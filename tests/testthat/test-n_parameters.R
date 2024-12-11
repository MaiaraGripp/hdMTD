test_that("n_parameters function works as expected", {
  # Test with default parameters
  result_1 <- n_parameters(Lambda = c(1, 3), A = c(4, 8, 12))
  expect_equal(result_1, 16)

  # Test with single_matrix = TRUE
  result_2 <- n_parameters(Lambda = c(2, 4, 9), A = c(0, 1), single_matrix = TRUE)
  expect_equal(result_2, 6)

  # Test with indep_part = FALSE
  result_3 <- n_parameters(Lambda = c(2, 4, 9), A = c(0, 1), indep_part = FALSE)
  expect_equal(result_3, 8)

  # Test with single_matrix = TRUE and indep_part = FALSE
  result_4 <- n_parameters(Lambda = c(2, 4, 9), A = c(0, 1), single_matrix = TRUE, indep_part = FALSE)
  expect_equal(result_4, 4)

  result_5 <- n_parameters(Lambda = c(2, 4, 9), A = c(0, 1), indep_part = FALSE, zeta=1)
  expect_equal(result_5, 4)

  result_6 <- n_parameters(Lambda = c(2, 4, 9), A = c(0,1,2), indep_part = TRUE)
  expect_equal(result_6, 23)

  result_7 <- n_parameters(Lambda = c(2, 4, 9, 20), A = c(0,1,2), indep_part = TRUE, zeta=2)
  expect_equal(result_7, 18)
})
