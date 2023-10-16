test_that("sx function calculates sx correctly", {

  set.seed(1)
  Lambda <- c(1, 3, 5)
  S <- c(1,5)
  A <- c(0,1)
  lenA <- length(A)
  x_S <- c(1,0)
  mu <- 0.1
  alpha <- 0.2
  xi <- 0.3
  MTD <- MTDmodel(Lambda,A)
  X <- perfectSample(MTD,N=500)

  cTab <- countsTab(X,d=5)
  freqTab <- freqTab(S=S,j=3,A=A,countsTab=cTab)

  expected_sx <- c(0.1236801, 0.1454461)


  result_sx <- sx(S=S, freqTab, lenA, x_S, mu, alpha, xi)

  expect_true(all.equal(result_sx, expected_sx, tolerance = 1e-5))
})
