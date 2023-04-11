#' Bayesian Information Criterion (BIC) of a MTD Markov Chain.
#'
#' Calculates the log maximum likelihood for each possible set of lags (with elements from \eqn{1 to d})
#'  of sizes \eqn{1,2, \dots, l} and penalizes by the number of parameters.
#'  \deqn{ -log ML + n_parameters x N x c }. See function [sparseMarkov::n_parameters()].
#'
#' @param X Markov chain.
#' @param A The states space.
#' @param S A set that includes all relevant lags. If only the upper bound d for the chains order in known set S=1:d.
#' @param l A upper threshold for the number of elements in the relevant lag sets.
#' @param xi The BIC constant. Defaulted to 1/2. Smaller xi (near 0) reduces the impact of overparameterization.
#' @importFrom utils combn
#'
#' @return The BIC of a Markov Chain for each possible set of lags with size 1,2...,l.
#' @export
#'
BIC_l <- function(X,A,S,l,xi=1/2){
  #checks
  if(length(S) < 2  ||
     !is.numeric(S) ||
     any( S%%1 != 0) ){stop("S must be a vector of at least 2 integer numbers.")
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
  while ( is.na(l) || !is.numeric(l) || l%%1 != 0 || l>length(S) ) {
    cat("l value is not valid. l should be a positive integer lower or equal to the number of elements in S.")
    l <- readline(prompt = "Please enter a valid l : ")
    l <- suppressWarnings(as.numeric(l))
  }
  while ( is.na(xi) || !is.numeric(xi) || xi <= 0 ) {
    cat("BIC constant xi value is not valid. xi should be a positive number.")
    xi <- readline(prompt = "Please enter a valid xi: ")
    xi <- suppressWarnings(as.numeric(xi))
  }
  #\.
  #gatherin inputs
  A <- sort(A)
  S <- sort(S)
  d <- max(S)
  lenS <- length(S)
  base <- shapeSample(X,d)
  #\.

  if(l==1){
    tryCombs <- matrix( c( rep(0,lenS), rep(n_parameters(1,A)*log(length(X))*xi,lenS) ),
                        byrow = T,nrow = 2 )
    colnames(tryCombs) <- S
    rownames(tryCombs) <- c("log_ML","penalty")
    for (k in S) {
      b <- base_Sja(k,j=NULL,A,base,complete = FALSE)
      tryCombs[1,k] <- -sum(b$Nxa_Sj*log(b$qax_Sj))
    }
    penalizedML <- apply(tryCombs,2,sum)

  }else{#if l>1
    tryCombs <- list()
    for (i in 1:l) {
      ncombs <- choose(lenS,i)
      tryCombs[[i]] <- matrix( rep(0,2*ncombs) ,byrow = T,nrow = 2 )
      tryCombs[[i]][2,] <- n_parameters(Lambda=(1:i),A)*log(length(X))*xi
      aux <- apply(t(combn(S,i)), 1, paste0,collapse=",")
      for ( k in 1:ncombs ) {
        G <- as.numeric(unlist(strsplit(aux[k], ",")))
        b <- base_Sja(G,j=NULL,A,base,complete = FALSE)
        tryCombs[[i]][1,k] <- -sum(b$Nxa_Sj*log(b$qax_Sj))
      }
      colnames(tryCombs[[i]]) <- aux
      rownames(tryCombs[[i]]) <- c("log_ML","penalty")
    }
    penalizedML <- sapply(tryCombs, apply, 2,sum)
  }
  penalizedML
}
