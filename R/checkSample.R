#' Checks a sample
#'
#' Checks if a sample is a suitable argument for some functions within
#'the package.
#'
#' @param X A vector or a single-column data frame with a sample of a MTD chain.
#'
#' @return Returns the sample as a vector or identifies any possible sample problems.
#'
checkSample <- function(X){
  if(is.data.frame(X)){
    if(ncol(X)!=1)stop("X must be a single chain so multiple columns are not accepted.")
    if(nrow(X)<=1)stop("Insufficient sample size.")
    X <- X[,1]
  }
  if (length(X)<=1)stop("Insufficient sample size.")
  if ( !is.numeric(X) ) { stop("X must be a numeric dataset.") }
  if ( length(ncol(X)) != 0 ) { stop("X must have only 1 dimension.") }
  if ( any(is.na(X)) ) { stop("NA values are not allowed in the sample.") }
  if ( length(unique(X)) == 1 ) { stop("All elements in the sample are the same.") }
  X
}

