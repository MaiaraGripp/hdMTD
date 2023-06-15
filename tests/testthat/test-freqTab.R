#should not work
test_that("countsTab works", {
  cTab <- countsTab(X=c(0,0,1,1,0,1,1,1,0,0,0),d=2)
  expect_error(freqTab(S=NULL,j=NULL,A=c(0,1),cTab))
})
test_that("countsTab works", {
  cTab <- countsTab(X=c(0,0,1,1,0,1,1,1,0,0,0),d=2)
  expect_error(freqTab(S=2,j=NULL,A=c(0,2),cTab))
})
test_that("countsTab works", {
  cTab <- countsTab(X=c(0,0,1,1,0,1,1,1,0,0,0),d=3)
  expect_error(freqTab(S=2,j=4,A=c(0,1),cTab))
})
test_that("countsTab works", {
  cTab <- countsTab(X=c(0,0,1,1,0,1,1,1,0,0,0),d=3)
  colnames(cTab) <- c("a","b","c","d","f")
  expect_error(freqTab(S=2,j=NULL,A=c(0,1),cTab))
})


#Should Work
test_that("freqTab works", {
  cTab <- countsTab(X=c(0,0,1,1,0,1,1,1,0,0,0),d=2)
  fTab <- freqTab(S=2,j=NULL,A=c(0,1),cTab)

  expect_equal(fTab$Nxa_Sj, c(1,3,3,2))
  expect_equal(fTab$Nx_Sj, c(4,4,5,5))
  expect_equal(fTab$qax_Sj, c(1,3,3,2)/c(4,4,5,5))
})
