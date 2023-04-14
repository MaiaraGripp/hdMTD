#' A function for inference in MTD Markov chains with BIC.
#'
#' Estimates a relevant lag set \eqn{\Lambda} of MTD models using BIC.
#'
#' @param X A Markov chain.
#' @param A The states space.
#' @param d An upper threshold for the chains order.
#' @param S A set of relevant lags, if empty S=\eqn{1,2,\dots, d}.
#' @param l An upper bound for the number of elements to be returned in the estimated relevant lag set.
#' @param xi The BIC constant. Defaulted to 1/2. Smaller xi `(near 0)` reduces the impact of overparameterization.
#' @param warning If TRUE may return warnings.
#'
#' @details See function [sparseMarkov::BIC_l].
#'
#' @return Returns an estimation, using BIC, of the relevant lag set with size \eqn{1,2, \dots l} and the relevant lag set with smallest BIC independent of size.
#' @export
sparseMarkov_BIC <- function(X,A=NULL,d,S=1:d,l,xi=1/2,warning=FALSE){
  #Checking inputs
  if(length(S) < 2  ||
     !is.numeric(S) ||
     any( S%%1 != 0) ){stop("S must be a vector of at least 2 integer numbers.")}

  if( !is.numeric(d) || d<2 || (d %% 1)!=0  ||  d<max(S)){
    stop("The order d must be an integer number greater than 2 or the greatest elements in S.")
  }
  if(length(A)==0){
    if(warning==TRUE){
      warning("States space A not informed. Code will set A <- sort(unique(X)).")
    }
    A <- unique(X)
  }else{
    if( !is.numeric(A) ||
        length(A)<=1   ||
        length(dim(A))!=0 )stop("States space A must be a numeric vector with at least two values.")
  }
  while ( is.na(l) || !is.numeric(l) || l%%1 != 0 || l>d || l>length(S) ) {
    cat("l value is not valid. l should be a positive integer lower or equal to d or the number of elements in S.")
    l <- readline(prompt = "Please enter a valid l : ")
    l <- suppressWarnings(as.numeric(l))
  }
  while ( is.na(xi) || !is.numeric(xi) || xi <= 0 ) {
    cat("BIC constant xi value is not valid. xi should be a positive number.")
    xi <- readline(prompt = "Please enter a valid xi: ")
    xi <- suppressWarnings(as.numeric(xi))
  }
  A <- sort(A)

    pML <- BIC_l(X=X,A=A,S=sort(S),l=l,xi=xi)
    smallest <- names(unlist(pML))[order(unlist(pML))][1]
    estLambda <- sapply(sapply(pML, orderedNames), dplyr::first )
    estLambda <- append(estLambda,smallest)
    names(estLambda) <- c(paste0("l=",1:l),"smallest")

  estLambda
}

orderedNames <- function(x){
  names(x[order(x)])
}
