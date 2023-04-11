#' Perfectly samples a MTD Markov chain.
#'
#' Can sample a MTD Markov Chain with stationary distribution.
#'
#' @title Perfect Sample.
#' @param MTD a canonical MTD object.
#' @param N the sample size. If NULL sample size will be set to 1000.
#' @return Returns a sample from o MTD model.
#'
#' @export perfectSample
perfectSample <- function(MTD,N=NULL){

  UseMethod("perfectSample")
}

#' @export
perfectSample.MTD <- function(MTD,N=NULL){
  if(length(N)==0){
    warning("Sample size N not informed. N will be set to 1000.")
    N <- 1000
  }
  checkMTD(MTD) #checks if MTD object is correctly set
  if( N<=max(MTD$Lambda) ||
      !is.numeric(N) ||
      length(N)!=1 ||
      N%%1!=0 )stop("Sample size N needs a integer number greater than max(Lambda)")

  a_k <- cumsum(MTD$lambdas)
  acum_p0 <- cumsum(MTD$p0)
  chain <- NULL
  repeat{
    chain <- c(NA,chain)
    Yt <- NULL #times chain
    Yt[1] <- 1
    Kt <- 1
    cont <- 1

    x <- MTD$A[which(acum_p0>runif(1))[1]] #the sample value to be used when sampling independently
    while (Kt>0) { #when this happens the last two Yt entrances will be identical
      u <- runif(1)
      Kt <- ifelse( (which(a_k>u)[1]-1)==0 , 0 , MTD$Lambda[which(a_k>u)[1]-1] )
      cont <- cont+1
      Yt[cont] <- Yt[cont-1]+Kt
      if( !is.na(chain[Yt[cont]]) ){
        x <- chain[Yt[cont]]
        Yt[cont+1] <- Yt[cont]
        break #breaks if chain has already been sampled for this Yt
      }
    }
    lenYt <- length(Yt)-1
    chain[Yt[lenYt]] <- x #sampled from p0 if Kt=0 or repeated from chain if there was already a symbol
    if(lenYt > 1){
      for (i in lenYt:2) {
        Kt <- Yt[i]-Yt[i-1]
        p_Kt <- MTD$p_j[[which(MTD$Lambda==Kt)]]
        p_KtRow <- which(MTD$A==chain[Yt[i]])#carefull a needs to be sorted
        #cant use rownames cause if p_j is inputed they can be missing
        acum_p_KtRow <- cumsum(p_Kt[p_KtRow,])
        chain[Yt[i-1]] <- MTD$A[which(acum_p_KtRow>runif(1))[1]]
      }
    }
    if( all( !is.na(chain[1:max(MTD$Lambda)]) ) ){ break }
  }
  chain <- chain[1:max(MTD$Lambda)]

  for( i in 1:(N-max(MTD$Lambda)) ){
    chain <- c(NA,chain)
    u <- runif(1)
    Kt <- ifelse( (which(a_k>u)[1]-1)==0, 0, MTD$Lambda[which(a_k>u)[1]-1] )
    chain[1] <-
      if(Kt==0){
        MTD$A[which(acum_p0>runif(1))[1]]
      }else{
        p_Kt <- MTD$p_j[[which(MTD$Lambda==Kt)]]
        p_KtRow <- which(MTD$A==chain[1+Kt])
        acum_p_KtRow <- cumsum(p_Kt[p_KtRow,])
        chain[1] <- MTD$A[which(acum_p_KtRow>runif(1))[1]]
      }
  }
  chain
}

#' @export
perfectSample.default <- function(MTD,N){
  print("The implemented perfect sample algorithm works only for MTD class objetcs.")
}

#export
#print.tildeMTD_01 <- function(x){
 # y <- unlist(MTD$perfectSample)
 # print(y)
#}

#export
#summary.tildeMTD_01 <- function(x){
 # list("P"=MTD$P,"lambdas"=MTD$lambdas,"p0"=MTD$p0,"Lambda"=MTD$Lambda,"A"=MTD$A,"signal"=MTD$signal,"perfectSample"=MTD$perfectSample)
#}

