#' Estimated transition probabilities
#'
#' Maximum likelihood estimators (MLE) for a MTD Markov chain with relevant lag set \code{S}.
#'
#' @param X A vector or single-column data frame containing a chain sample.
#' @param S A numeric vector of positive integers. Typically, \code{S} represents a set of relevant lags.
#' @param matrixform Logical. If \code{TRUE}, the probability estimates outputted by this function
#' come in the form of a stochastic matrix.
#' @param A A vector with positive integers representing the state space. If not informed,
#' this function will set \code{A=unique(X)}.
#' @param warning Logical. If \code{TRUE}, informs the user if the state space was set as \code{A=unique(X)}..
#'
#' @return The MLE for a given chain sample with relevant lags set \code{S}.
#' @export
#'
#' @examples
#' X <- testChains[,3]
#' probs(X,S=c(1,30))
#' probs(X,S=c(1,15,30))
#'
probs <- function(X,S,matrixform=FALSE,A=NULL,warning=FALSE){
## Check X, S, A, matrixform
  X <- checkSample(X)
  if(length(S) < 1  ||
     !is.numeric(S) ||
     any( S%%1 != 0) ){stop("S must be a vector of at least 2 integer numbers.")
  }
  if(length(A)==0){
    if(warning==TRUE){
      warning("States space A not informed. Code will set A <- sort(unique(X)).")
    }
    A <- sort(unique(X))
  }else{A <- sort(A)}
  if( length(A)<=1   ||
      any(A%%1 !=0)   )stop("States space A must be a numeric vector with at least two integers.")
  if ( !all( unique(X) %in% A ) ) {
    stop("Verify the state space to ensure it includes all states that appear in the sample.")
  }
  if(!is.logical(matrixform)){stop("matrixform should be a logical parameter that determines the return format of the function.")}

  S <- sort(S,decreasing = TRUE)
  lenS <- length(S)
  base <- countsTab(X,S[1])
  base <- freqTab(S=S,A=A,countsTab=base,complete=TRUE)
  probs <- cbind(apply(base[,1:lenS],1,paste0,collapse=""),base[,lenS+1])
  probs <- cbind(probs,base$qax_Sj)
  names(probs) <- c(paste("past_{",paste0(-S,collapse = ","),"}"),"a_0","p(a|past)")

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
