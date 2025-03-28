test_that("checkSample works correctly", {
  # Create a valid sample for testing
  valid_sample <- sample(c(0,1),100,T)

  # Test for a valid sample
  expect_equal(checkSample(valid_sample), valid_sample)

  # Test for a valid data frame sample
  valid_df_sample <- data.frame(X = valid_sample)
  expect_equal(checkSample(valid_df_sample), valid_sample)

  # Test for insufficient sample size
  insufficient_sample <- 1
  expect_error(checkSample(insufficient_sample))

  # Test for non-numeric sample
  non_numeric_sample <- c(1, 2, "a", 4, 5)
  expect_error(checkSample(non_numeric_sample))

  # Test for sample with NA values
  sample_with_na <- c(1, 2, 3, NA, 5)
  expect_error(checkSample(sample_with_na))

  # Test for sample with all elements the same
  sample_all_same <- c(1, 1, 1, 1, 1)
  expect_error(checkSample(sample_all_same))

})
