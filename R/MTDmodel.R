#' Creates an MTD model
#'
#' Given a set of parameters generates a Mixture transition distribution model as an object of class MTD.
#'
#' @param Lambda The relevant lag set should consist of positive integers and will be sorted from smallest
#'  to greatest. The smallest number represents the latest (most recent) time in the past, and the greatest
#'   number represents the earliest time in the past.
#' @param A The states space.
#' @param w0 The weight of the independent distribution, must be a value in `[0,1)`.
#' @param w_j A vector of weights for the distributions p_j. Values must be in the range `[0, 1)`. The order of weights must match the order
#'  of the sorted Lambda set.
#' @param p_j A list with lenL matrices of size `lenA x lenA`. Here, lenL represents the number of elements in Lambda, and lenA represents the number
#'  of elements in A.
#' @param p0 A vector with the independent distribution part of an MTD model. If not informed and indep_part=TRUE, the distribution will be sampled
#' from a uniform distribution. If indep_part=FALSE, the distribution will be set to zero.
#' @param single_matrix If TRUE all p_j matrices are equal.
#' @param indep_part If FALSE the independent distribution is set to zero.
#'
#' @details This MTD object can be used by functions such as [oscillation()], which retrieves the model's oscillation, and [perfectSample()],
#'  which will perfectly sample a MTD Markov chain with the model's parameters
#'
#' @return An MTD model as a MTD object.
#' @importFrom stats runif
#' @export
#'
#' @examples
#' MTDmodel(Lambda=c(1,3),A=c(4,8,12))
#'
#' MTDmodel(Lambda=c(2,4,9),A=c(0,1),w0=0.05,w_j=c(0.35,0.2,0.4),
#' p_j=list(matrix(c(0.5,0.7,0.5,0.3),ncol=2)),p0=c(0.2,0.8),single_matrix=TRUE)
#'
#' MTDmodel(Lambda=c(2,4,9),A=c(0,1),w0=0.05,
#' p_j=list(matrix(c(0.5,0.7,0.5,0.3),ncol=2)),single_matrix=TRUE,indep_part=FALSE)
MTDmodel <- function(Lambda,
                     A,
                     w0 = NULL,
                     w_j = NULL,
                     p_j = NULL,
                     p0 = NULL,
                     single_matrix = FALSE,
                     indep_part = TRUE)
{
  #testing Lambda and A
  if(!is.numeric(Lambda) ||
     any(Lambda<=0) ||
     !all(Lambda%%1==0)||
     !is.vector(Lambda))stop("Lambda must be a numeric vector with positive integers.")
  if( length(A)<=1   ||
      any(A%%1 !=0)   )stop("States space A must be a numeric vector with at least two integers.")

  #sorting parameters
  if( all(sort(Lambda)!=Lambda) & length(w_j)!=0 )warning("The Lambda set will be ordered from smallest to greatest, be carefull with matching order of w_j accordingly.")
  Lambda <- sort(Lambda)
  A <- sort(A)

  #getting info from parameters
  lenA <- length(A)
  lenL <- length(Lambda)
  lenAL <- lenA^lenL
  tilde_A <- seq(1:lenA)

  #init variables
  lambdas <- numeric(lenL+1)

  #test for p0 (independent distribution)
  if(length(p0)!=0){
    if( !is.numeric(p0) ||
        round(sum(p0),3)!=1 ||
        !all(p0>0) ||
        length(p0)!=lenA )stop("p0 must be a vector of size ",lenA, " consisting of numeric, nonnegative values that sum up to 1.")
  }else{
    p0 <- stats::runif(lenA)
    p0 <- p0/sum(p0)
  }
  names(p0) <- c(paste0("p_0(",A,")"))

  #rewrite p0 and set w0=0 if indep_part=FALSE
  if(!is.logical(indep_part)) stop("Argument indep_part must be TRUE OR FALSE.")
  if(!indep_part){
      p0 <- rep(0,lenA)
      if(length(w0)!=0 && w0!=0 )warning("Since indep_part=FALSE w0 and p0 will be set to 0.")
      w0 <- 0
    }


  #tests for lambdas
  if(length(w0)!=0){
    if( length(w0) > 1 ||
        !is.numeric(w0) ||
        w0 < 0 ||
        w0 >= 1) stop("w0 must be a number in the range [0,1) .")
    #if w0 passes test...
    lambdas[1] <- w0
  }

  if(length(w_j)!=0){
    if( !is.numeric(w_j) ||
        !lenL==length(w_j) ||
        !all(w_j>0)) stop(paste0("w_j must be a vector of size ",lenL, ", numeric and nonnegative."))
    #if w_j passes tests...
    if(length(w0)==0){
      lambdas <- c(1-sum(w_j),w_j)
    }else{
      lambdas <- c(lambdas[1],w_j)
    }
  }else{ #if w_j=NULL...
    if(length(w0)==0){
      lambdas <- runif(lenL+1)
      lambdas <- lambdas/sum(lambdas)
    }else{
      aux <- runif(lenL)
      aux <- (1-lambdas[1])*(aux/sum(aux))
      lambdas <- c(lambdas[1],aux)
    }
  }

  if(round(sum(lambdas),6)!=1) stop("The sum of weights w0 + [w_j] must add up to 1. If indep_part=FALSE, w0 will be set to 0, so be careful in case w_j was inputted.")

  names(lambdas) <- c("lam0",paste0("lam-", Lambda) )

  #tests for p_j
  if(!is.logical(single_matrix))stop("Argument single_matrix must be TRUE OR FALSE.")
  if(single_matrix){
    if( !(length(p_j)%in%c(0,1)) )stop("If single_matrix=TRUE p_j must be a list with a single matrix or NULL.")
  }else{
    if( !(length(p_j)%in%c(0,lenL)) )stop(paste0("If single_matrix=FALSE p_j must be NULL or a list with ",lenL, " matrices."))
  }
  if(length(p_j)!=0){
    if(!is.list(p_j)) stop("p_j must be a list.")
    aux <- do.call(rbind,p_j)
    if( !all(round(apply(aux, 1, sum),2)==1)  ||
        !all(aux>0) ||
        !is.numeric(aux) ||
        ncol(aux)!=lenA ||
        any(sapply(p_j,dim)!=lenA)) stop(paste0("p_j must be a list of stochastic matrices ", lenA, "x",lenA))
  }else{ #if p_j=NULL
    p_j <- list()
    for (j in 1:lenL) {
      R <- matrix(runif(lenA^2),ncol=lenA,nrow = lenA)
      R <- R/apply(R,1,sum)
      p_j[[j]] <- R #each [[]] is a matrix p_j, lenA x lenA, forall j in Lambda
      colnames(p_j[[j]])=rownames(p_j[[j]]) <- A
    }
  }
  if(single_matrix && lenL>=2){
    for (j in 2:lenL) {
      p_j[[j]] <- p_j[[1]]
    }
  }
  names(p_j) <- paste0("p_-",Lambda)

  subx <- try(expand.grid(rep(list(tilde_A),lenL)),silent = TRUE) #all possible sequences x_{Lambda}
  if(class(subx)=="try-error"){stop(paste0("For lenL=",lenL," the dataset with all pasts sequences (x of size lenL) with elements of A is too large."))}
  subx <- subx[,order(lenL:1)]

  P <- matrix(0,ncol = lenA,nrow = lenAL)
  if (lenL==1) {
    for (i in 1:lenAL) {
      P[i,] <- lambdas%*%rbind(p0,p_j[[1]][i,])
    }
    rownames(P) <- A
  }else{
    for (i in 1:lenAL) {
      aux <- matrix(0,ncol=lenA,nrow = lenL)
      for (j in 1:lenL) {
        aux[j,] <- p_j[[j]][subx[i,(lenL+1-j)],] # the lines in aux are each from a different p_j
      }
      P[i,] <- lambdas%*%rbind(p0,aux) #\sum_{j\in Lambda and 0} lambda_j*p_j(.|subx[i,j])
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
              p_j=p_j,
              p0=p0,
              Lambda=Lambda,
              A=A)
  class(MTD) <- "MTD"
  MTD
}

#' @export
print.MTD <- function(x, ...){
  print(x[seq_along(x)])
  return(invisible(x))
}
