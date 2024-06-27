#' Perfectly samples an MTD Markov chain
#'
#' Samples an MTD Markov Chain from the stationary distribution.
#'
#' @param MTD An MTD object, see [MTDmodel()] for properly generating a MTD object.
#' @param N The sample size. If \code{NULL} sample size will be set to 1000.
#' @return Returns a sample from an MTD model (the first element is the most recent).
#'
#' @details This perfect sample algorithm requires that the MTD model has
#' an independent distribution and a positive weight (i.e., \code{MTD$lambdas["lam0"]>0} which
#' means \eqn{\lambda_0>0}).
#'
#' @export perfectSample
#'
#' @examples
#' perfectSample(MTDmodel(Lambda=c(1,4),A=c(0,2)),N=200 )
#' perfectSample(MTDmodel(Lambda=c(2,5),A=c(1,2,3)),N=1000 )
#'
perfectSample <- function(MTD,N=NULL){

  UseMethod("perfectSample")
}

#' @export
perfectSample.MTD <- function(MTD,N=NULL){
## Check inputs
  if(length(N)==0){
    warning("Sample size N not informed. N will be set to 1000.")
    N <- 1000
  }
  checkMTD(MTD)
  if( sum(MTD$p0)==0 || MTD$lambdas[1]==0){stop("The provided MTD model does not have an independent distribution or a positive weight for it, but the algorithm that this function uses for perfect sampling requires both.")}
  if( N<=max(MTD$Lambda) ||
      !is.numeric(N) ||
      length(N)!=1 ||
      N%%1!=0 )stop("Sample size N needs a integer number greater than max(Lambda)")


  a_k <- cumsum(MTD$lambdas)
  acum_p0 <- cumsum(MTD$p0)
  chain <- NULL
  repeat{ # Runs the perfect sample algorithm until we obtain a sample that contains a length d sequence
# of observations with no gaps between them.

    chain <- c(NA,chain) # initiate first chain element as NA
    Yt <- NULL # initiate vector with sampling positions
    Yt[1] <- 1
    Kt <- 1
    cont <- 1

    x <- MTD$A[which(acum_p0>runif(1))[1]] #x is the independent sampling result (x~p0)
    while (Kt>0) { # Kt>0 implies Yt>Yt-1 (Kt=0 iff Yt=Yt-1).
      u <- runif(1)
      Kt <- ifelse( (which(a_k>u)[1]-1)==0 , 0 , MTD$Lambda[which(a_k>u)[1]-1] ) # tells us from which distribution to sample.
      cont <- cont+1
      Yt[cont] <- Yt[cont-1]+Kt # Yt is the position of the sample that is relevant for sampling chain[Yt-1]

      if( !is.na(chain[Yt[cont]]) ){ # verifies if there is already something sampled chain[Yt], if TRUE:
        x <- chain[Yt[cont]] # x is replaced with whatever was in chain[Yt]
        Yt[cont+1] <- Yt[cont]
        break # stops the loop because we've found a sampled past.
      }
      # If there was nothing previously sampled (in other loops) this loop only stops in Kt=0
    }
    lenYt <- length(Yt)-1
    chain[Yt[lenYt]] <- x #sampled from p0 if Kt=0, or taken from chain if there was something sampled
    if(lenYt > 1){ # if we did not get to sample from p0 right from the start we need to sample all pasts in Yt
      for (i in lenYt:2) {
        Kt <- Yt[i]-Yt[i-1]
        p_Kt <- MTD$pj[[which(MTD$Lambda==Kt)]]
        p_KtRow <- which(MTD$A==chain[Yt[i]])# !A needs to be sorted just like lines in pj
        acum_p_KtRow <- cumsum(p_Kt[p_KtRow,])
        chain[Yt[i-1]] <- MTD$A[which(acum_p_KtRow>runif(1))[1]]
      }
    }
    if( all( !is.na(chain[1:max(MTD$Lambda)]) ) ){ break } #stops when we get a whole sequence of length d without gaps
  }
  chain <- chain[1:max(MTD$Lambda)]

  for( i in 1:(N-max(MTD$Lambda)) ){ #completes the remaining sample size with initial sample coming "chain".
    chain <- c(NA,chain)
    u <- runif(1)
    Kt <- ifelse( (which(a_k>u)[1]-1)==0, 0, MTD$Lambda[which(a_k>u)[1]-1] )
    chain[1] <-
      if(Kt==0){
        MTD$A[which(acum_p0>runif(1))[1]]
      }else{
        p_Kt <- MTD$pj[[which(MTD$Lambda==Kt)]]
        p_KtRow <- which(MTD$A==chain[1+Kt])
        acum_p_KtRow <- cumsum(p_Kt[p_KtRow,])
        chain[1] <- MTD$A[which(acum_p_KtRow>runif(1))[1]]
      }
  }
  chain #X1...XN
}

#' @export
perfectSample.default <- function(MTD,N){
  print("The implemented perfect sample algorithm works only for MTD class objetcs.")
}
