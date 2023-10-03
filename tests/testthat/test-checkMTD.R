test_that("checkMTD works MTD", {
  expect_error(checkMTD(c(1,2)))
})
test_that("checkMTD works p0", {
  obj <- MTDmodel(c(1,2,3),c(0,1))
  obj$p0 <- NULL
  expect_error(checkMTD(obj))
  obj$p0 <- c(2,2)
  expect_error(checkMTD(obj))
  obj$p0 <- 0
  expect_error(checkMTD(obj))
  obj$p0 <- c(0,0)
  expect_no_error(checkMTD(obj))
})
test_that("checkMTD works lambdas", {
  obj <- MTDmodel(c(1,2,3),c(0,1))
  obj$lambdas <- NULL
  expect_error(checkMTD(obj))
  obj$lambdas <- c(0.5,0.5)
  expect_error(checkMTD(obj))
  obj$lambdas <- c(0.5,0.5,0.5,0.5)
  expect_error(checkMTD(obj))
})
test_that("checkMTD works Lambda", {
  obj <- MTDmodel(c(1,2,3),c(0,1))
  obj$Lambda <- NULL
  expect_error(checkMTD(obj))
})
test_that("checkMTD works p_j", {
  obj <- MTDmodel(c(1,2,3),c(0,1))
  obj$p_j <- c(1,2,3)
  expect_error(checkMTD(obj))
})
