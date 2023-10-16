test_that("dTV_pj works correctly", {
  # Create a valid stochastic matrix p_j for testing
  p_j <- matrix(c(0.1, 0.2, 0.7, 0.4, 0.5, 0.1, 0.5, 0.3, 0.2), nrow = 3)
  rows <- matrix(c(1, 2), ncol = 2)
  which_max <- FALSE

  # Test for a valid p_j
  result <- dTV_pj(p_j, rows, which_max)

  # Check the result
  expect_equal(result, sum(abs(p_j[1,]-p_j[2,]))/2)

  # Test for which_max = TRUE
  which_max <- TRUE
  rows <- matrix(c(1,1,2,2,3,3),ncol=2)
  result <- dTV_pj(p_j, rows, which_max)

  # Check the result
  expect_equal(result[[1]], sum(abs(p_j[1,]-p_j[3,]))/2)
  expect_equal(result[[2]], "1x3")

  # Add more test cases as needed
})
