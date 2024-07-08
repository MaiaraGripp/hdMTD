#' Number of Parameters in an MTD Model
#'
#' Calculates the number of parameters in an MTD model with a relevant lag set Lambda and state space A.
#'
#' @param Lambda A numeric vector of positive integers representing the relevant lag set.
#' The elements will be sorted from smallest to greatest. The smallest number represents the latest
#'  (most recent) time in the past, and the greatest number represents the earliest time in the past.
#' @param A A vector with positive integers representing the state space.
#' @param single_matrix Logical. If \code{TRUE}, the chain sample is thought to come from an MTD model
#' where the stochastic matrices \eqn{p_j} are constant across all lags \eqn{j\in \Lambda}. So there
#' are fewer parameters in the penalization term
#' @param indep_part Logical. If \code{FALSE} there is no independent distribution and \eqn{\lambda_0=0} which
#' reduces the number of parameters in the penalization term.
#'
#' @return The number of parameters in an MTD model.
n_parameters <- function(Lambda, A, single_matrix = FALSE, indep_part = TRUE){
  lenA <- length(A)
  lenL <- length(Lambda)
  n_parameters <- lenL-1 #number of lambda params if lam0=0
  if (indep_part) { #then lambda0 > 0
    n_parameters <- n_parameters+lenA  # number of params with lambda0 (+1) and p_0 distribution (+lenA-1)
  }
  if (single_matrix) {
    n_parameters <- n_parameters+lenA*(lenA-1)
  }else{
    n_parameters <- n_parameters+lenA*(lenA-1)*lenL
  }
  n_parameters
}
