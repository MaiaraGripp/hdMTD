#' Estimated transition probabilities
#'
#' Maximum likelihood estimators for a MTD Markov chain with relevant lag set S.
#'
#' @param X A MTD Markov Chain.
#' @param S The estimated relevant lags set.
#' @param A The state space. "A" only needs to be informed if X does not already contain all elements of "A"
#' @param warning If TRUE may return warnings.
#'
#' @return The MLE for a given MTD Markov Chain sample with relevant lags set S.
#' @export
#'
#' @examples
#' X <- perfectSample(MTDmodel(Lambda=c(1,7),A=c(1,2)),N=500)
#' probs(X,S=c(1,7))
#' probs(X,S=c(1,7,9))
#'
probs <- function(X,S,A=NULL,warning=FALSE){
  X <- checkSample(X)
  if(length(S) < 2  ||
     !is.numeric(S) ||
     any( S%%1 != 0) ){stop("S must be a vector of at least 2 integer numbers.")
  }
  if(length(A)==0){
    if(warning==TRUE){
      warning("States space A not informed. Code will set A <- sort(unique(X)).")
    }
    A <- unique(X)
  }
  if( !is.numeric(A) ||
      length(A)<=1   ||
      length(dim(A))!=0 )stop("States space A must be a numeric vector with at least two values.")
  if ( !all( unique(X) %in% A ) ) {
    stop("Check the states space, it include all states that occur in the sample.")
  }

  S <- sort(S,decreasing = TRUE)
  lenS <- length(S)
  base <- countsTab(X,S[1])
  base <- freqTab(S=S,A=A,countsTab=base,complete=TRUE)
  probs <- cbind(apply(base[,1:lenS],1,paste0,collapse=""),base[,lenS+1])
  probs <- cbind(probs,base$qax_Sj)
  names(probs) <- c(paste("past_{",paste0(-S,collapse = ","),"}"),"a_0","p(a|past)")
  probs
}
