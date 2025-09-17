test_that("empirical_probs function estimates probabilities correctly", {

  # Relevant lag set
  Lambda <- c(1, 3)
  # State space
  A <- c(1, 2)
  #MTD object
  MTD <- MTDmodel(Lambda, A)
  # Create a sample MTD Markov Chain
  sample_chain <-perfectSample(MTD, N = 500)

  # Calculate probabilities
  result <- empirical_probs(sample_chain, Lambda)

  # Check if the result matches the expected output
  expect_true(all(colnames(result)==c("past_{ -3,-1 }","a","p(a|past)")))
  expect_equal(result[,1], c('11','11','12','12','21','21','22','22'))
  expect_equal(result[,2], c( 1, 2, 1, 2, 1, 2, 1, 2))
  expect_type(result[,3],"double")
  expect_true(all(result[,3]>=0))
  expect_true(all(result[,3]<=1))
})
