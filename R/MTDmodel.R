
#' Creates an MTD model
#'
#' Given a set of parameters it generates an MTD model as an object of class MTD.
#'
#' @param Lambda A numeric vector of positive integers representing the relevant lag set.
#' The elements will be sorted from smallest to greatest. The smallest number represents the latest
#'  (most recent) time in the past, and the greatest number represents the earliest time in the past.
#' @param A A vector with nonnegative integers representing the state space.
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
#' @importFrom methods is
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
  # The user may set indep_part to FALSE by providing p0=0 or p0 as a vector of 0
  if( length(p0) != 0 & all(p0 == 0) & indep_part == TRUE ){
    indep_part <- FALSE
    warning("Since p0=0 indep_part will be set to FALSE")
  }

  # Validates MTDmodel inputs
  check_MTDmodel_inputs(Lambda, A, lam0, lamj, pj, p0, single_matrix, indep_part)

  # Warns if Lambda needs to be sorted and lamj was provided since their order must match
  if( any( Lambda != sort(Lambda) ) & length(lamj) != 0 ) {
    warning("The Lambda set will be ordered from smallest to largest, the user should match the order of lamj accordingly.")
  }
  # Sorting
  Lambda <- sort(Lambda)
  A <- sort(A)

  lenA <- length(A)
  lenL <- length(Lambda)
  lenAL <- lenA^lenL

  # If indep_part is FALSE sets p0 <- c(0,...,0)
  if( !indep_part ){
    p0 <- rep(0,lenA)
    # If the user provides lam0!=0 but indep_part=FALSE, warns that lam0 <- 0
    if( length( lam0 ) != 0 && lam0 != 0 ) {
      warning("Since indep_part=FALSE lam0 will be set to 0 and p0 will be a vector of 0.")
    }
    lam0 <- 0
  }
  # If p0 is not provided samples a distribution from uniforms
  if ( is.null(p0) ) {
    p0 <- stats::runif(lenA)
    p0 <- p0/sum(p0)
  }
  names(p0) <- c(paste0("p0(",A,")"))

  lambdas <- numeric(lenL+1) # init vector for MTD weights
  if( !is.null(lam0) ){
    lambdas[1] <- lam0
  }
  if( !is.null(lamj) ){
    if( is.null(lam0) ){
      lambdas <- c(1-sum(lamj),lamj)
    } else {
      lambdas <- c(lambdas[1],lamj)
    }
  } else {
    if( is.null(lam0) ){
      lambdas <- runif(lenL+1)
      lambdas <- lambdas/sum(lambdas)
    } else {
      aux <- runif(lenL)
      aux <- (1-lambdas[1])*(aux/sum(aux))
      lambdas <- c(lambdas[1],aux)
    }
  }
  # check if weights add to 1
  if( round(sum(lambdas), 5) != 1) stop("The sum of weights lam0 + sum(lamj) must be equal to 1.")
  names(lambdas) <- c("lam0",paste0("lam-", Lambda) )

  # generate pj matrices if pj is not provided
  if( is.null(pj) ){
    pj <- list()
    if ( !single_matrix ) { # generates lenL matrices
      for ( j in seq_len(lenL) ) {
        R <- matrix(runif(lenA^2), ncol = lenA, nrow = lenA)
        R <- R/apply(R, 1, sum)
        pj[[j]] <- R
        colnames(pj[[j]]) <- rownames(pj[[j]]) <- A
      }
    } else { # generates 1 matrix
      R <- matrix(runif(lenA^2), ncol = lenA, nrow = lenA)
      R <- R/apply(R, 1, sum)
      pj[[1]] <- R
      colnames(pj[[1]]) <- rownames(pj[[1]]) <- A
    }
  }
  if( single_matrix && lenL >= 2 ){ # repeats the informed or generated matrix for all lags
    for (j in 2:lenL) {
      pj[[j]] <- pj[[1]]
    }
  }
  names(pj) <- paste0("p-",Lambda)

  # Calculating P (the transition matrix)
  subx <- try(expand.grid(rep(list(seq(1:lenA)),lenL)),silent = TRUE) #all possible sequences x_{Lambda}
  if( inherits(subx,"try-error") ){
    stop(paste0("For length(Lambda)=",lenL," the dataset with all pasts sequences (x of length(Lambda)) with elements of A is too large."))
    }
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

