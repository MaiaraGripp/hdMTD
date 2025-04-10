test_that("checkMTD works correctly", {
  # Create a valid MTD object for testing
  valid_MTD <- list(
    Lambda = c(1, 3),
    A = c(2, 3, 5),
    p0 = c(0.2, 0.5, 0.3),
    lambdas = c(0.01, 0.49, 0.5),
    pj = list(
      matrix(c(0.1, 0.2, 0.7, 0.4, 0.5, 0.1, 0.5, 0.3, 0.2), nrow = 3),
      matrix(c(0.1, 0.2, 0.7, 0.4, 0.5, 0.1, 0.5, 0.3, 0.2), nrow = 3)
      )
  )
  class(valid_MTD) <- "MTD"

  # Test for a valid MTD
  expect_silent(checkMTD(valid_MTD))

  valid_MTD2 <- valid_MTD
  valid_MTD2$p0 <- 0
  # Test for a valid MTD
  expect_silent(checkMTD(valid_MTD2))

  # Test for MTD with invalid Lambda
  invalid_MTD_Lambda <- valid_MTD
  invalid_MTD_Lambda$Lambda <- c(1, -2, 3)
  expect_error(checkMTD(invalid_MTD_Lambda))
  invalid_MTD_Lambda$Lambda <- c(1, 2, 3)
  expect_error(checkMTD(invalid_MTD_Lambda))

  # Test for MTD with invalid p0
  invalid_MTD_p0 <- valid_MTD
  invalid_MTD_p0$p0 <- c(0.2, 0.3, 0.6)
  expect_error(checkMTD(invalid_MTD_p0))
  invalid_MTD_p0$p0 <- c(0.4, 0.6)
  expect_error(checkMTD(invalid_MTD_p0))
  invalid_MTD_p0$p0 <- c(-1.6, 0.3,0.3)
  expect_error(checkMTD(invalid_MTD_p0))
  invalid_MTD_p0$p0 <- c(0, 0)
  expect_error(checkMTD(invalid_MTD_p0))
  # Add more tests as needed
})

test_that("checkMTD works p0", {
  obj <- MTDmodel(c(1,2,3),c(0,1))
  obj$p0 <- NULL
  expect_error(checkMTD(obj))
  obj$p0 <- c(2,2)
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
test_that("checkMTD works pj", {
  obj <- MTDmodel(c(1,2,3),c(0,1))
  obj$pj <- c(1,2,3)
  expect_error(checkMTD(obj))
})
