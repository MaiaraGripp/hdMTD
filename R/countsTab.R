#' Counts sequences of length d+1 in a sample
#'
#' Creates a tibble containing all unique sequences of length \code{d+1} found in
#' the sample, along with their absolute frequencies.
#'
#' @details The function generates a tibble with \code{d+2} columns. In the first
#' \code{d+1} columns, each row displays a unique sequence of size \code{d+1}
#' observed in the sample. The last column, called \code{Nxa}, contains the
#' number of times each of these sequences appeared in the sample.
#'
#' The number of rows in the output varies between \eqn{1} and \eqn{|A|^{d+1}},
#' where \eqn{|A|} is the number of unique states in \code{X}, since it depends
#' on the number of unique sequences that appear in the sample.
#'
#' @param X A numeric vector, a single-column data frame, or a list with a
#' sample from a Markov chain. The first element must be the most recent observation.
#' @param d A positive integer specifying the number of elements in each sequence,
#' which will be \code{d+1}. Typically, \code{d} represents the chain order or
#' serves as an upper limit for it.
#'
#' @return  A tibble with all observed sequences of length \code{d+1} and their
#' absolute frequencies.
#' @export
#' @importFrom dplyr %>%
#'
#' @examples
#' countsTab(c(1,2,2,1,2,1,1,2,1,2),3)
#'
#' # Using test data.
#' countsTab(testChains[, 1],2)
#'
countsTab <-function(X,d){
  # Validate and process the input sample (X)
  X <- checkSample(X)
  # Checks if d is a positive integer
  if (!is.numeric(d) || d %% 1 != 0 || d <= 0)
      stop("d must be a positive integer.")
  # Ensures the sample size is sufficient
  if (length(X) <= d )
      stop("The sample size must be greater than d.")

  n <- length(X)
  d1 <- d + 1
  X <- rev(X)  # Reverses the sample so the output is read from left to right


  if(n - d >= d1){
      # If there are d1 or more sequences of size d1 in the sample
      X_list <- vector("list", d1) # List to store d1 matrices with d1 columns
      for (i in 1:d1) {
          aux <- (n - (i - 1)) %% d1
          X_list[[i]] <- matrix(X[i:(n - aux)], ncol = d1, byrow = TRUE)
      }
  } else {
      # If there are fewer sequences, store the available ones
      X_list <- vector("list", n-d)
      for (i in 1:(n - d)) {
          X_list[[i]] <- X[i:(i + d)]
      }
  }

  # Bind all matrices/sequences of X_list in a single matrix
  XTab <- do.call(rbind, X_list)
  colnames(XTab) <- c(paste0("x", seq(d, 1)), "a")

  # Convert to tibble and count occurrences
  XTab <- dplyr::as_tibble(XTab) %>%
    dplyr::group_by(across(everything())) %>%
    dplyr::summarise(Nxa = dplyr::n(), .groups = "drop")
  XTab
}

# dplyr seems to be faster than data.table when d is large

