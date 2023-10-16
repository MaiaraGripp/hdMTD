test_that("Test for the PI function", {
  Lambda <- c(1,3,5)
  A <- c(0,1)
  MTD <- MTDmodel(Lambda,A)
  X <- perfectSample(MTD,N=1000)
  lenX <- length(X)
  d <- 3
  S <- c(1,5)
  cTab <- countsTab(X,5)
  fTab <- freqTab(S=S,j=3,A=A,countsTab=cTab)
  fTabSj <- freqTabSj(S=Lambda,j=5,freqTab = fTab,lenX=lenX,d=d)
  x_S <- c(0,1)
  # Replace with your test values
  result <- PI(S=S, fTabSj, x_S, lenX, d)

  expect_type(result,"double")
  expect_length(result, length(A))
})
