#' Check a MTD object.
#'
#' Checks if a MTD object has all the needed parameters.
#'
#' @param MTD An MTD object
#'
checkMTD <- function(MTD){
  #list + MTD class
  if(!is.list(MTD)    ||
     class(MTD)!="MTD" )stop("MTD must be an object of class MTD which is a list with especific parameters. The use of MTDmodel() function in order to create MTD is recommended.")
  #Lambda
  if( any(MTD$Lambda<=0)   ||
     !all(MTD$Lambda%%1==0)||
     !is.vector(MTD$Lambda) )stop("Lambda must be a numeric vector with positive integers.")
  if(all(sort(MTD$Lambda)!=MTD$Lambda))stop("The Lambda set must be ordered from smallest to greatest, be carefull with matching the order of weights w_j accordingly.")
  lenL <- length(MTD$Lambda)

  #A
  if( length(MTD$A)<=1   ||
      !is.vector(MTD$A)  ||
      any(MTD$A%%1!=0)    )stop("States space A must be a numeric vector with at least two integers.")
  if(all(sort(MTD$A)!=MTD$A))stop("The states space A must be ordered from smallest to largest.")
  lenA <- length(MTD$A)

  #p0
  if( !is.numeric(MTD$p0) ||
      !is.vector(MTD$p0)  ||
      !all(MTD$p0>=0)      )stop("p0 must be a numeric nonnegative vector.")
  if( !length(MTD$p0) %in% c(1,lenA) )stop("p0 must be either 0 or a vector of size ",lenA)
  if( round(sum(MTD$p0),3)!=1 & sum(MTD$p0)!=0 )stop("Either each element in p0 is 0 or they must sum up to 1.")

  #lambdas
  if( !is.numeric(MTD$lambdas)     ||
      round(sum(MTD$lambdas),3)!=1 ||
      !all(MTD$lambdas>=0)         ||
      length(MTD$lambdas)!=(lenL+1) )stop("lambdas must be a vector of size ",lenL+1, "(the same size of the vector of relevant lags Lambda + 1), with nonnegative numbers that must add up to 1. The first element of the lambdas vector is the weight for the independent distribution p0, if your MTD model doesn't have and independent distribution set lambdas[1]=0.")

  if(!is.list(MTD$pj)                ||
     length(MTD$pj)!=lenL            ||
     !all(sapply(MTD$pj, is.matrix)) ||
     !all(sapply(MTD$pj,dim)==c(lenA,lenA))
     )stop(paste0("pj must be a list with ", lenL," stochastic matrices ", lenA, "x",lenA))
  aux <- do.call(rbind,MTD$pj)
  if(!is.numeric(aux)                      ||
     !all(round(apply(aux, 1, sum),3)==1)  ||
     !all(aux>=0) )stop(paste0("pj must be a list with ", lenL," stochastic matrices ", lenA, "x",lenA,". In other words, each matrix row must sum up to 1."))
}
