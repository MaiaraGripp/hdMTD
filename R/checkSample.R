#' Checks a sample
#'
#' Checks if a sample is a suitable argument for some functions within the package.
#'
#' @param X  A numeric vector, a single-column data frame, a list, or a matrix
#' with a single row or a single column.
#'
#' @return Returns the sample as a vector or identifies any possible sample problems.
#'
checkSample <- function(X){
  # Allow matrices with a single column or a single row
  if (is.matrix(X)) {
      if (ncol(X) == 1 || nrow(X) == 1) {
          X <- as.vector(X)  # Convert to vector
      } else {
          stop("X must be a vector, a single-column data frame, a list, or a matrix
               with a single row or a single column.")
      }
  }

  # Convert single-column data frame to vector
  if (is.data.frame(X)) {
      if (ncol(X) != 1)
          stop("X must be a single chain; multiple columns are not accepted.")
      if (nrow(X) <= 1)
          stop("Insufficient sample size.")
      X <- as.vector(X[, 1])
  }

  # Convert list to vector if necessary
  if (is.list(X)) {
      X <- unlist(X)
  }
  # Basic validation checks
  if (length(X) <= 1)
      stop("Insufficient sample size.")
  if (!is.numeric(X))
      stop("X must be a numeric dataset.")
  if (any(is.na(X)))
      stop("NA values are not allowed in the sample.")
  if (length(unique(X)) == 1)
      stop("The sample must contain at least two distinct values.")

  X
}

