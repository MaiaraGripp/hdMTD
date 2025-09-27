test_that("MTDest function works as expected", {
  # Create a sample MTD chain
  set.seed(1)
  Lambda <- c(1, 10)
  A <- c(0, 1)
  MTD <- MTDmodel(Lambda, A, 0.01)
  X <- perfectSample(MTD, 2000)
  S <- c(1, 10)

  init <- list('p0' = MTD$p0, 'lambdas' = MTD$lambdas, 'pj' = MTD$pj)
  init$p0 <- c(0.5, 0.5)

  # Test with default parameters
  result_1 <- MTDest(X = X, S = S, M = 0.05, init = init)
  expect_true(class(result_1) == "MTDest")
  expect_true(all(names(coef(result_1)) %in% c("lambdas", "pj", "p0")))
  expect_equal(sum(p0(result_1)),1)
  expect_equal(sum(lambdas(result_1)),1)

  init$p0 <- c(0)
  expect_warning(expect_error(MTDest(X = X, S = S, M = 0.05, init = init)))

  init$lambdas[2] <- init$lambdas[2] + init$lambdas[1]
  init$lambdas[1] <- 0
  # Test with default parameters
  result_1 <- MTDest(X = X, S = S, M = 0.05, init = init)
  expect_true(class(result_1) == "MTDest")
  expect_true(all(names(coef(result_1)) %in% c("lambdas", "pj", "p0")))
  expect_true(all(c(p0(result_1), lambdas(result_1)[1])==c(0,0,0)))
  # Test with custom parameters
  init2 <- list('p0' = c(0.4, 0.6),
               'lambdas' = c(0.05, 0.45, 0.5),
               'pj' = list(
                 matrix(c(0.2, 0.8, 0.45, 0.55), byrow = TRUE, ncol = 2),
                 matrix(c(0.25, 0.75, 0.3, 0.7), byrow = TRUE, ncol = 2))
               )
  result_2 <- MTDest(X = X, S = S, init = init2, iter = TRUE)
  expect_true(class(as.MTD(result_2)) == "MTD")
})
