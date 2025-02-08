#' Checks a class MTD object
#'
#' Verifies if an object of class MTD is correctly structured, containing all
#'  the necessary parameters, and if these parameters satisfy their respective
#'  constraints.
#'
#' @details This package includes a function called [MTDmodel()], which outputs
#' a properly structured MTD object that does not require additional checks.
#' However, since the user can create the MTD object manually, this [checkMTD()]
#' is used within the functions that use MTD objects as inputs to prevent errors.
#' The correct structure and instances of a class MTD object are detailed in the
#' [MTDmodel()] function documentation.
#'
#' @param MTD An object of class MTD.
#'
#' @importFrom methods is
#'
checkMTD <- function(MTD){

  # Verifies if the object is a list with class MTD.
  if (!is.list(MTD))
    stop("MTD must be a list.")
  if (!is(MTD, "MTD"))
    stop("MTD must be an object of class MTD.")

  # Checks if Lambda is a numeric vector of positive integers in ascending order.
  if (any(MTD$Lambda <= 0) || !all(MTD$Lambda%%1 == 0) || !is.vector(MTD$Lambda))
    stop("Lambda must be a numeric vector of positive integers.")
  if (any(sort(MTD$Lambda) != MTD$Lambda))
    stop("Lambda must be sorted in ascending order.")

  lenL <- length(MTD$Lambda)

  # Checks if A is a numeric vector of integers (length ≥ 2), sorted in ascending order.
  if (length(MTD$A) <= 1 || !is.vector(MTD$A) || any(MTD$A%%1 != 0))
    stop("State space A must be a numeric vector containing at least two integers.")
  if (any(sort(MTD$A) != MTD$A))
    stop("State space A must be sorted in ascending order.")

  lenA <- length(MTD$A)

  # Checks if p0 is a numeric nonnegative vector of length 1 or length(A), summing to 1.
  if (!is.numeric(MTD$p0) || !is.vector(MTD$p0) || !all(MTD$p0 >= 0))
    stop("p0 must be a nonnegative numeric vector.")
  if (!length(MTD$p0) %in% c(1, lenA))
    stop(paste0("p0 must be either a scalar 0 or a numeric vector of length ", lenA, "."))
  if (round(sum(MTD$p0), 3) != 1 & sum(MTD$p0) != 0)
    stop("The elements in p0 must either sum to 1 or all be 0.")

  # Checks if lambdas is a numeric nonnegative vector of length (length(Lambda) + 1) that sums to 1.
    if (!is.numeric(MTD$lambdas) || round(sum(MTD$lambdas), 3) != 1 ||
        !all(MTD$lambdas >= 0) || length(MTD$lambdas) != (lenL + 1))
        stop(paste0(
        "lambdas must be a vector of length ", lenL + 1, " (the number
        of relevant lags in Lambda plus 1), consisting of nonnegative numbers
        that sum to 1. The first element of the lambdas vector is the weight for
        the independent distribution p0, if your MTD model does not include an
        independent distribution, set lambdas[1] to 0."
        ))

  # Checks if pj is a list with length(Lambda) elements, each containing a stochastic matrix of size length(A) x length(A).
  if(!is.list(MTD$pj) || length(MTD$pj) != lenL ||
     !all(sapply(MTD$pj, is.matrix)) || !all(sapply(MTD$pj,dim) == c(lenA,lenA)))
    stop(paste0("pj must be a list with ", lenL, " stochastic matrices ", lenA,"x",lenA))
  aux <- do.call(rbind, MTD$pj)
  if(!is.numeric(aux) || !all(round(apply(aux, 1, sum),3)==1) || !all(aux>=0))
    stop(paste0("pj must be a list with ", lenL, " stochastic matrices ", lenA, "x",lenA,".
    In other words, each matrix row must sum up to 1."))
}
