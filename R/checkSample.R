#' Check a Markov chain sample
#'
#' Checks if a Markov chain sample is suitable for use
#'
#' @param X A vector with a Markov Chain.
#'
#' @return Returns adjusted chain sample and possible sample errors.
#'
checkSample <- function(X){
  if(is.data.frame(X)){
    if(ncol(X)!=1)stop("X must be a single Markov Chain so multiple columns are not accepted.")
    if(nrow(X)<=1)stop("Insufficient sample size.")
    X <- X[,1]
  }
  if ( !is.numeric(X) ) { stop("X must be a numeric dataset.") }
  if ( length(ncol(X)) != 0 ) { stop("X must have only 1 dimension.") }
  if ( any(is.na(X)) ) { stop("NA values are not allowed in the sample.") }
  if ( length(unique(X)) == 1 ) { stop("All elements in the sample are the same.") }
  X
}

