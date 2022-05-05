test_that("checkSample works", {
  expect_error(checkSample(c("a","b")))
})
test_that("checkSample works", {
  expect_error(checkSample(c(1,1,1)))
})
test_that("checkSample works", {
  expect_error(checkSample(c("a","b")))
})
test_that("checkSample works", {
  expect_error(checkSample(matrix(c(1,2,3,4),2,2)))
})
test_that("checkSample works", {
  expect_error(checkSample(c(1,"NA",2)))
})
test_that("checkSample works", {
  expect_silent(checkSample(c(1,2,1)))
})
