#should not work
test_that("shapeSample works", {
  expect_error(shapeSample(c(1,1,1,0,1,0,0,0,1),10))
})
test_that("shapeSample works", {
  expect_error(shapeSample(c(1,1,1,0,1,0,0,0,1),1))
})
test_that("shapeSample works", {
  expect_error(shapeSample(c(1,1,1,0,1,0,0,0,1),15.5))
})

A <- c(1,8,9,11)
n <- 100
testVar <- rbind(sample(A,n,replace = TRUE),sample(A,n,replace = TRUE))
d=5
test_that("shapeSample works", {
  expect_error(shapeSample(testVar,d))
})

testVar <- sample(A,n,replace = TRUE)
testVar[1]=NA
test_that("shapeSample works", {
  expect_error(shapeSample(testVar,d))
})

#Should work
A <- c(1,3,5)
n <- 10000
testVar <- sample(A,n,replace = TRUE)
d <- 4
test_that("shapeSample works", {
  expect_true(dim(shapeSample(testVar,d))[1] <= length(A)^(d+1) &
                dim(shapeSample(testVar,d))[2] == d+1+1)
  expect_equal(sum(shapeSample(testVar,d)$Nxa),n-d)
})
A <- c(1,8,9,11)
n <- 5000
testVar <- sample(A,n,replace = TRUE)
d=5
test_that("shapeSample works", {
  expect_true(dim(shapeSample(testVar,d))[1] <= length(A)^(d+1) &
                dim(shapeSample(testVar,d))[2] == d+1+1)
  expect_equal(sum(shapeSample(testVar,d)$Nxa),n-d)
})

