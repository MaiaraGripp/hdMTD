test_that("hdMTD_BIC function works as expected", {
  # Create a sample MTD chain
  set.seed(1)
  X <- perfectSample(MTDmodel(c(1, 5), c(0, 1)), 5000)
  d <- 5

  # Test when minl is equal to maxl
  result_1 <- hdMTD_BIC(X = X, d = d, minl = 2, maxl = 2, xi = 0.5)
  expect_true(is.numeric(result_1))
  expect_true(all(result_1 %in% 1:d))

  # Test when minl is less than maxl
  result_2 <- hdMTD_BIC(X = X, d = d, minl = 1, maxl = 3, xi = 0.1, BICvalue = TRUE)
  expect_true(is.numeric(result_2))
  expect_true(all(names(result_2) %in% 1:d))

  # Test with byl = TRUE
  result_3 <- hdMTD_BIC(X = X, d = d, minl = 1, maxl = 3, xi = 0.2, byl = TRUE)
  expect_true(is.character(result_3))
  expect_length(result_3, 4)
  for (i in 1:(length(result_3)-1)) {
  expect_true(all( as.integer(unlist(strsplit(result_3[i],","))) %in% 1:d) )
  }

  # Test with warning = TRUE
  expect_warning(hdMTD_BIC(X = X, d = d, minl = 1, maxl = 2, warning = TRUE))
})
