#' Estimated transition probabilities
#'
#' Computes the Maximum likelihood estimators (MLE) for an MTD Markov chain with
#' relevant lag set \code{S}.
#'
#' @param X A vector or single-column data frame containing a chain sample.
#' @param S A numeric vector of positive integers representing a set of relevant lags.
#' @param matrixform Logical. If \code{TRUE}, the output is formatted as a stochastic transition matrix.
#' @param A A numeric vector representing the state space, consisting of distinct integers.
#' If not provided, this function will set \code{A = unique(X)}.
#' @param warning Logical. If \code{TRUE}, the function warns the user when the state space is automatically set as \code{A = unique(X)}.
#'
#' @return The MLE estimates for the given chain sample, conditioned on the relevant lags in \code{S}.
#' If \code{matrixform = TRUE}, the result is a stochastic matrix.
#' @export
#'
#' @examples
#' X <- testChains[, 3]
#' probs(X, S = c(1, 30))
#' probs(X, S = c(1, 15, 30))
#'
probs <- function( X, S, matrixform=FALSE, A=NULL, warning=FALSE ){

  X <- checkSample(X)
  check_probs_inputs( X, S, matrixform, A, warning )

  if(length(A) == 0) { A <- sort(unique(X)) } else { A <- sort(A) }

  S <- sort(S,decreasing = TRUE)
  lenS <- length(S)

  base <- countsTab(X,max(S))
  base <- freqTab(S=S,A=A,countsTab=base,complete=TRUE)

  probs <- data.frame(
    apply(base[, 1:lenS], 1, paste0, collapse = ""),
    base[, lenS + 1],
    base$qax_Sj
  )
  names(probs) <- c(paste("past_{",paste0(-S,collapse = ","),"}"),"a","p(a|past)")

  if(matrixform){
    Pest <- probs$`p(a|past)`
    dim(Pest) <- c(length(A),length(A)^lenS)
    Pest <- t(Pest)
    colnames(Pest) <- A
    rownames(Pest) <- unique(probs[,1])
    probs <- Pest
  }
  probs
}
