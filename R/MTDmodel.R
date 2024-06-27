#' Creates an MTD model
#'
#' Given a set of parameters it generates an MTD model as an object of class MTD.
#'
#' @param Lambda A numeric vector of positive integers representing the relevant lag set.
#' The elements will be sorted from smallest to greatest. The smallest number represents the latest
#'  (most recent) time in the past, and the greatest number represents the earliest time in the past.
#' @param A A vector with positive integers representing the state space.
#' @param lam0 The weight of the independent distribution, must be a value in `[0,1)`.
#' @param lamj A vector of weights for the distributions in the argument \code{pj}. Values must be
#' in the range `[0, 1)`. The first element in \code{lamj} must be the weight of the first
#' element of the list \code{pj} and so on.
#' @param pj A list with \code{length(Lambda)} stochastic matrices.
#' @param p0 A vector with the independent distribution part of an MTD model. If not informed
#' and argument \code{indep_part=TRUE}, the distribution will be sampled from a uniform distribution.
#' If \code{indep_part=FALSE}, then there is no independent distribution and \code{p0} entries will
#' be set to zero. If you enter \code{p0=0}, then \code{indep_part} will be set to \code{FALSE}.
#' @param single_matrix Logical. If \code{TRUE}, all matrices in list \code{pj} are equal.
#' @param indep_part Logical. If \code{FALSE} there is no independent distribution and \code{p0}
#' is set to zero.
#'
#' @details This MTD object can be used by functions such as [oscillation()], which retrieves the
#'  model's oscillation, and [perfectSample()], which will perfectly sample an MTD Markov chain
#'  with the model's parameters.
#'
#' @return An MTD model as an MTD object.
#' @importFrom stats runif
#' @export
#'
#' @examples
#' MTDmodel(Lambda=c(1,3),A=c(4,8,12))
#'
#' MTDmodel(Lambda=c(2,4,9),A=c(0,1),lam0=0.05,lamj=c(0.35,0.2,0.4),
#' pj=list(matrix(c(0.5,0.7,0.5,0.3),ncol=2)),p0=c(0.2,0.8),single_matrix=TRUE)
#'
#' MTDmodel(Lambda=c(2,4,9),A=c(0,1),lam0=0.05,
#' pj=list(matrix(c(0.5,0.7,0.5,0.3),ncol=2)),single_matrix=TRUE,indep_part=FALSE)
MTDmodel <- function(Lambda,
                     A,
                     lam0 = NULL,
                     lamj = NULL,
                     pj = NULL,
                     p0 = NULL,
                     single_matrix = FALSE,
                     indep_part = TRUE)
{
## Check Lambda and A
  if(!is.numeric(Lambda) ||
     any(Lambda<=0) ||
     !all(Lambda%%1==0)||
     !is.vector(Lambda))stop("Lambda must be a numeric vector with positive integers.")
  if( length(A)<=1   ||
      any(A%%1 !=0)   )stop("States space A must be a numeric vector with at least two integers.")
  # Sorting
  if( all(sort(Lambda)!=Lambda) & length(lamj)!=0 )warning("The Lambda set will be ordered from smallest to greatest, be carefull with matching order of lamj accordingly.")
  Lambda <- sort(Lambda)
  A <- sort(A)
## Gathering some inputs
  lenA <- length(A)
  lenL <- length(Lambda)
  lenAL <- lenA^lenL

## Check indep_part and p0 (independent distribution)
  if(!is.logical(indep_part)) stop("Argument indep_part must be TRUE or FALSE.")
  if( length(p0)!=0 & all(p0==0) & indep_part==TRUE ){
    indep_part <- FALSE
    warning("Since p0=0 indep_part will be set to FALSE")
  }
  if(!indep_part){
    p0 <- rep(0,lenA)
    if(length(lam0)!=0 && lam0!=0 )warning("Since indep_part=FALSE lam0 will be set to 0 and p0 will be a vector of 0.")
    lam0 <- 0
  }
  if(length(p0)!=0){ #p0!=NULL
    if( !is.numeric(p0) ||
        !all(p0>=0) ||
        length(p0)!=lenA)stop("p0 must be 0 or a vector of length ",lenA, " consisting of nonnegative numbers.")
    if( round(sum(p0),3)!=1 & sum(p0)!=0 )stop("Either each element in p0 is 0 or they must sum up to 1.")
  }else{
    p0 <- stats::runif(lenA)
    p0 <- p0/sum(p0)
  }
  names(p0) <- c(paste0("p0(",A,")"))

## init lambdas and check lam0
  lambdas <- numeric(lenL+1)
  if(length(lam0)!=0){
    if( length(lam0) > 1 ||
        !is.numeric(lam0) ||
        lam0 < 0 ||
        lam0 >= 1) stop("lam0 must be a number in the range [0,1) .")
    lambdas[1] <- lam0
  }
## check lamj
  if(length(lamj)!=0){
    if( !is.numeric(lamj) ||
        !lenL==length(lamj) ||
        !all(lamj>0) ) stop(paste0("lamj must be a vector of length ",lenL, ", numeric and nonnegative."))
    if(length(lam0)==0){
      lambdas <- c(1-sum(lamj),lamj)
    }else{
      lambdas <- c(lambdas[1],lamj)
    }
  }else{ #if lamj=NULL, sample:
    if(length(lam0)==0){
      lambdas <- runif(lenL+1)
      lambdas <- lambdas/sum(lambdas)
    }else{
      aux <- runif(lenL)
      aux <- (1-lambdas[1])*(aux/sum(aux))
      lambdas <- c(lambdas[1],aux)
    }
  }
  if(round(sum(lambdas),6)!=1) stop("The sum of weights lam0 + [lamj] must sum up to 1. If indep_part=FALSE, lam0 will be set to 0, so be careful in case lamj was inputted.")
  names(lambdas) <- c("lam0",paste0("lam-", Lambda) )

## check single_matrix and pj
  if(!is.logical(single_matrix))stop("Argument single_matrix must be TRUE OR FALSE.")
  if(single_matrix){
    if( !(length(pj)%in%c(0,1)) )stop("If single_matrix=TRUE pj must be a list with a single matrix or NULL.")
  }else{
    if( !(length(pj)%in%c(0,lenL)) )stop(paste0("If single_matrix=FALSE pj must be NULL or a list with ",lenL, " matrices."))
  }
  if(length(pj)!=0){
    if(!is.list(pj)) stop("pj must be a list.")
    aux <- do.call(rbind,pj)
    if( !all(round(apply(aux, 1, sum),2)==1)  ||
        !all(aux>=0) ||
        !is.numeric(aux) ||
        ncol(aux)!=lenA ||
        any(sapply(pj,dim)!=lenA)) stop(paste0("pj must be a list of stochastic matrices ", lenA, "x",lenA))
  }else{ #if pj=NULL
    pj <- list()
    for (j in 1:lenL) {
      R <- matrix(runif(lenA^2),ncol=lenA,nrow = lenA)
      R <- R/apply(R,1,sum)
      pj[[j]] <- R #each [[]] is a matrix pj, lenA x lenA, forall j in Lambda
      colnames(pj[[j]])=rownames(pj[[j]]) <- A
    }
  }
  if(single_matrix && lenL>=2){
    for (j in 2:lenL) {
      pj[[j]] <- pj[[1]]
    }
  }
  names(pj) <- paste0("p-",Lambda)
## Calculating P
  #initiating
  subx <- try(expand.grid(rep(list(seq(1:lenA)),lenL)),silent = TRUE) #all possible sequences x_{Lambda}
  if(class(subx)=="try-error"){stop(paste0("For length(Lambda)=",lenL," the dataset with all pasts sequences (x of length(Lambda)) with elements of A is too large."))}
  subx <- subx[,order(lenL:1)]
  P <- matrix(0,ncol = lenA,nrow = lenAL)

  if (lenL==1) {
    for (i in 1:lenAL) { # runs in all pasts of length(Lambda) with elements of A (the lines of P)
      P[i,] <- lambdas%*%rbind(p0,pj[[1]][i,])
    }
    rownames(P) <- A
  }else{
    for (i in 1:lenAL) { # runs in all pasts of length(Lambda)  with elements of A (the lines of P)
      aux <- matrix(0,ncol=lenA,nrow = lenL)
      for (j in 1:lenL) {
        aux[j,] <- pj[[j]][subx[i,(lenL+1-j)],] # the lines in aux are each from a different pj
      }
      P[i,] <- lambdas%*%rbind(p0,aux) #\sum_{j\in Lambda and 0} lambda_j*pj(.|subx[i,j])
    }
  }
  colnames(P) <- A
  if(lenL>1){
    subx <- as.matrix(expand.grid(rep(list(A),lenL))) #calculate subx again so the actual values of A can be the names of P
    subx <- subx[,order(lenL:1)]
    rownames(P) <- apply(subx, 1, paste0,collapse="")
  }

  MTD <- list(P=P,
              lambdas=lambdas,
              pj=pj,
              p0=p0,
              Lambda=Lambda,
              A=A)

  class(MTD) <- "MTD"
  MTD
}

#' @export
print.MTD <- function(x, ...){
  ind <- seq_along(x)
  if ( all(x$p0==0) ) {
    ind <- ind[-4]
  }
  print(x[ind])
  return(invisible(x))
}
