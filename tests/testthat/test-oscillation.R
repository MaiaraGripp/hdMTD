test_that("oscillation function works as expected", {

  # Test with an MTD model object
  set.seed(1)
  MTD_model <- MTDmodel(Lambda = c(1, 3), A = c(4, 8, 12))
  result_1 <- oscillation(MTD_model)
  expect_equal(round(as.vector(result_1),3), c(0.018, 0.108))

  # Test with a chain sample, providing S and A
  Lambda <- c(1,10)
  A <- c(0,1)
  MTD <- MTDmodel(Lambda,A)
  X <-perfectSample(MTD,N=1000)
  result_2 <- oscillation(X, S = c(1, 10))
  expect_equal(round(as.vector(result_2),3), c(0.056,0.227))

})
