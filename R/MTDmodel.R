#' Creates a Mixture Transition Distribution (MTD) Model
#'
#' Generates an MTD model as an object of class \code{MTD} given a set of parameters.
#'
#' @param Lambda A numeric vector of positive integers representing the relevant lag set.
#' The elements will be sorted from smallest to greatest. The smallest number represents the latest
#'  (most recent) time in the past, and the largest number represents the earliest time in the past.
#' @param A A vector with nonnegative integers representing the state space.
#' @param lam0 A numeric value in `[0,1)`, representing the weight of the independent distribution.
#' @param lamj A numeric vector of weights for the transition probability matrices in \code{pj}.
#' Values must be in the range `[0, 1)`, and their sum with `lam0` must be equal to 1.
#' The first element in \code{lamj} must be the weight for the first element in \code{Lambda} and so on.
#' @param pj A list with \code{length(Lambda)} stochastic matrices, each of size `length(A) x length(A)`.
#' The first matrix in \code{pj} must refer to the first element in \code{Lambda} and so on.
#' @param p0 A probability vector for the independent component of the MTD model. If \code{NULL}
#' and \code{indep_part=TRUE}, the distribution will be sampled from a uniform distribution.
#' If \code{indep_part=FALSE}, then there is no independent distribution and \code{p0} entries will
#' be set to zero. If you enter \code{p0=0}, \code{indep_part} is set to \code{FALSE}.
#' @param single_matrix Logical. If \code{TRUE}, all matrices in list \code{pj} are identical.
#' @param indep_part Logical. If \code{FALSE}, the model does not include an independent distribution
#' and \code{p0} is set to zero.
#'
#' @details The resulting MTD object can be used by functions such as [oscillation()], which retrieves the
#'  model's oscillation, and [perfectSample()], which will sample an MTD Markov chain from its invariant
#'  distribution.
#'
#' @return A list of class \code{MTD} containing:
#' \describe{
#'   \item{`P`}{The transition probability matrix of the MTD model.}
#'   \item{`lambdas`}{A vector with MTD weights (`lam0` and `lamj`).}
#'   \item{`pj`}{A list of stochastic matrices defining conditional transition probabilities.}
#'   \item{`p0`}{The independent probability distribution.}
#'   \item{`Lambda`}{The vector of relevant lags.}
#'   \item{`A`}{The state space.}
#' }
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

  # Validates MTDmodel inputs
  check_MTDmodel_inputs(Lambda, A, lam0, lamj, pj, p0, single_matrix, indep_part)

  # Warns the user if lamj or pj were provided but Lambda needs to be sorted, since their elements order must match.
  if( any( Lambda != sort(Lambda) ) & length(lamj) != 0 ) {
    if( !is.null(lamj) || !is.null(pj)  ){
      warning("The Lambda set will be ordered from smallest to largest, the user should match the order of lamj accordingly.")
    }
  }
  # Sorting Lambda and A
  Lambda <- sort(Lambda)
  A <- sort(A)

  lenA <- length(A)
  lenL <- length(Lambda)
  lenAL <- lenA^lenL

  # If the user provides p0 to a zero vector, automatically set indep_part to FALSE
  if( !is.null(p0) && all(p0 == 0) && indep_part ){
    indep_part <- FALSE
    warning("Since p0=0, indep_part is set to FALSE.")
  }

  # If indep_part is FALSE, enforce p0 as a zero vector and set lam0 = 0
  if ( !indep_part ) {
    p0 <- rep(0, lenA)
    if( !is.null(lam0) && lam0 != 0 ) {
      warning("Since indep_part=FALSE, lam0 is set to 0 and p0 is a zero vector.")
    }
    lam0 <- 0
  }

  # If p0 is not provided, sample from a uniform distribution
  if ( is.null(p0) ) {
    p0 <- stats::runif(lenA)
    p0 <- p0 / sum(p0)
  }
  names(p0) <- paste0("p0(", A, ")")

  # Constructing MTD weight vector (lambdas)
  lambdas <- numeric(lenL+1)
  if( !is.null(lamj) ){
    if( is.null(lam0) ){
      lam0 <- 1 - sum(lamj)
    }
      lambdas <- c(lam0, lamj)
  } else {
    # Generate random weights ensuring they sum to 1
    if( is.null(lam0) ){
      lambdas <- runif(lenL+1)
      lambdas <- lambdas/sum(lambdas)
    } else {
      aux <- runif(lenL)
      aux <- (1-lam0)*(aux/sum(aux))
      lambdas <- c(lam0,aux)
    }
  }

  # Check if weights add to 1
  if(round(sum(lambdas), 5) != 1) {
    stop("The sum of weights lam0 + sum(lamj) must be equal to 1.")
  }
  names(lambdas) <- c("lam0",paste0("lam-", Lambda))

  # Generate stochastic matrices (pj) if not provided
  if( is.null(pj) ) {
    pj <- list()
    if ( !single_matrix ) { # Generates lenL matrices
      for ( j in seq_len(lenL) ) {
        R <- matrix(stats::runif(lenA^2), ncol = lenA, nrow = lenA)
        R <- R/rowSums(R)
        pj[[j]] <- R
        colnames(pj[[j]]) <- rownames(pj[[j]]) <- A
      }
    } else { # Generates 1 matrix
      R <- matrix(stats::runif(lenA^2), ncol = lenA, nrow = lenA)
      R <- R/rowSums(R)
      pj[[1]] <- R
      colnames(pj[[1]]) <- rownames(pj[[1]]) <- A
    }
  }
  if( single_matrix && lenL >= 2 ){
    # Repeats the informed or generated matrix for all lags
    for (j in 2:lenL) {
      pj[[j]] <- pj[[1]]
    }
  }
  names(pj) <- paste0("p-",Lambda)

# Calculating P (the transition matrix)
  # Generate all possible size lenL sequences with digits from 1 to lenA
  subx <- try(expand.grid(rep(list(seq_len(lenA)), lenL)), silent = TRUE)
  if( inherits(subx,"try-error") ){
    stop(paste0("For length(Lambda)=",lenL," the dataset with all pasts sequences (x of length(Lambda)) with elements of A is too large."))
    }
  subx <- subx[,order(lenL:1)]

  P <- matrix(0,ncol = lenA,nrow = lenAL)

  if (lenL==1) {
    for (i in seq_len(lenAL)) { # runs in all lines of P
      P[i, ] <- lambdas %*% rbind(p0, pj[[1]][i, ])
    }
    rownames(P) <- A
  }else{
    for (i in seq_len(lenAL)) { # runs in the lines of P (each refers to a sequence of size lenL with elements from A)
      aux <- matrix(0, ncol = lenA, nrow = lenL)
      for (j in seq_len(lenL)) {
        aux[j, ] <- pj[[j]][subx[i, (lenL + 1 - j)], ] # the lines in aux are each from a different pj
      }
      P[i, ] <- lambdas %*% rbind(p0, aux)
    }
  }
  colnames(P) <- A
  if(lenL > 1){
    subx <- as.matrix(expand.grid(rep(list(A), lenL))) # Now subx carries the actual values of A
    subx <- subx[, order(lenL:1)]
    rownames(P) <- apply(subx, 1, paste0, collapse="")
  }

  MTD <- list(
    P = P,
    lambdas = lambdas,
    pj = pj,
    p0 = p0,
    Lambda = Lambda,
    A = A
    )
  class(MTD) <- "MTD"
  MTD
}

#' @export
print.MTD <- function(x, ...) {
  ind <- seq_along(x)
  if (all(x$p0 == 0)) {
    ind <- ind[-which(names(x) == "p0")]
  }
  print(x[ind])
  return(invisible(x))
}

