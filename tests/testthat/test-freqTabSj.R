test_that("freqTabSj works correctly", {
  # Create sample data for testing
  Lambda <- c(1,3,5)
  S <- c(1, 5)
  j <- 3
  d <- 5
  A <- c(0,1)
  N <- 100
  X <- sample(A,N,T)
  cTab <- countsTab(X,d)
  fTab <- freqTab(S,j,A,cTab)

  # Test for a valid input
  result <- freqTabSj(S, j, fTab, lenX=length(X), d)

  # Check the result
  expect_equal(class(result)[1], "tbl_df")
  expect_equal(ncol(result),length(c(S,j))+1)
  expect_equal(sum(result$Nx_Sj),N-d)
  # Add more specific checks based on your expectations

  # Add more test cases as needed
})
