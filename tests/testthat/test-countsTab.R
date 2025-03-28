
test_that("countsTab works", {
  expect_error(countsTab(c(1,1,1,0,1,0,0,0,1),10))
})
test_that("countsTab works", {
  expect_error(countsTab(c(1,1,1,0,1,0,0,0,1),15.5))
})

#Should work
A <- c(1,8,9,11)
n <- 5000
testVar <- sample(A, n, replace = TRUE)
d = 5
test_that("countsTab works", {
  expect_true(dim(countsTab(testVar,d))[1] <= length(A)^(d + 1) &
                dim(countsTab(testVar,d))[2] == d + 2)
  expect_equal(sum(countsTab(testVar,d)$Nxa), n - d)
})

