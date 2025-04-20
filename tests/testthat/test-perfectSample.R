test_that("perfectSample function works as expected", {
  # Test with an MTD model object
  MTD_model <- MTDmodel(Lambda = c(1, 4), A = c(0, 2))
  result_1 <- perfectSample(MTD_model, N = 200)
  expect_length(result_1, 200)

})
test_that("perfectSample function is not set for non-MTD objects.", {
  # Test with a non-MTD object
  non_MTD_object <- c(0, 2, 0, 2, 0)
  expect_error(perfectSample(non_MTD_object, N = 200))
})
test_that("perfectSample function does not work for MTD objects that lack an independent distribution.", {
  # Test with a non-MTD object
  MTD_model <- MTDmodel(Lambda = c(1, 4), A = c(0, 2), indep_part = FALSE)
  expect_error(perfectSample(MTD_model, N = 200))
})


