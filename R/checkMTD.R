#' Checks if a MTD object has all the needed parameters.
#'
#' @param MTD An MTD object
#'
checkMTD <- function(MTD){
  if(!is.list(MTD)|| class(MTD)!="MTD")stop("MTD must be a object of class MTD which is a list of especific parameters. You can resort to MTDmodel function to create such object.")
  #testing Lambda and A
  if(!is.numeric(MTD$Lambda) ||
     any(MTD$Lambda<=0) ||
     !all(MTD$Lambda%%1==0)||
     !is.vector(MTD$Lambda))stop("Lambda must be a numeric vector with positive integers.")
  if(!is.numeric(MTD$A) ||
     any(MTD$A<0) ||
     !is.vector(MTD$A))stop("A must be a numeric vector with nonnegative numbers.")

  #sorting parameters
  if(all(sort(MTD$Lambda)!=MTD$Lambda))stop("Lambda set must be ordered from smallest to greater, be carefull with matching order of weights.")
  Lambda <- sort(MTD$Lambda)
  if(all(sort(MTD$A)!=MTD$A))stop("Stace space A must be ordered from smallest to greater.")
  MTD$A <- sort(MTD$A)

  lenL <- length(MTD$Lambda)
  lenA <- length(MTD$A)
  if( !is.numeric(MTD$p0) ||
      round(sum(MTD$p0),3)!=1 ||
      !all(MTD$p0>0) ||
      length(MTD$p0)!=lenA )stop("p0 must be a vector of size ",lenA, " numeric, nonnegative and must add up to 1.")

  if( !is.numeric(MTD$lambdas) ||
       round(sum(MTD$lambdas),3)!=1 ||
       !all(MTD$lambdas>0) ||
       length(MTD$lambdas)!=(lenL+1) )stop("p0 must be a vector of size ",lenL+1, " numeric, nonnegative and must add up to 1.")

   if(length(MTD$p_j)!=lenL)stop("p_j must be a list with ", lenL," stochastic matrices.")
  aux <- do.call(rbind,MTD$p_j)
  if( !all(round(apply(aux, 1, sum),3)==1)  ||
      !all(aux>0) ||
      !is.numeric(aux) ||
      ncol(aux)!=lenA ||
      any(sapply(MTD$p_j,dim)!=lenA)) stop(paste0("p_j must be a list of stochastic matrices ", lenA, "x",lenA))
}
