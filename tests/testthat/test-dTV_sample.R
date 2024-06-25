test_that("dTV_sample works correctly", {
  # Create sample data for testing
  X <- sample(c(0,1),200,replace = TRUE)
  S <- c(1, 4)
  j <- 3
  A <- c(0,1,2)
  lenA <- length(A)
  cTab <- countsTab(X,d=4)
  base <- freqTab(S,j,A,cTab,complete = TRUE)
  A_pairs <- matrix(c(0,1,0,2,1,2),byrow=T,ncol=2)
  x_S <- c(1, 0)

  # Test for a valid input
  result <- dTV_sample(S=S, j=j, lenA=lenA, base=base, A_pairs=A_pairs, x_S=x_S)

  # Check the result
  expect_equal(length(result), nrow(A_pairs))
  # Add more specific checks based on your expectations

  # Add more test cases as needed
})
