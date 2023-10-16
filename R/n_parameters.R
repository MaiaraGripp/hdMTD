#' Number of Parameters in an MTD Model
#'
#' Calculates the number of parameters in an MTD model with a relevant lag set Lambda and state space A.
#'
#' @param Lambda The relevant lag set.
#' @param A The state space.
#' @param single_matrix If TRUE all p_j matrices are equal.
#' @param indep_part If FALSE independent distribution is set to zero.
#'
#' @return The number of parameters in a MTD model.
n_parameters <- function(Lambda, A, single_matrix = FALSE, indep_part = TRUE){
  lenA <- length(A)
  lenL <- length(Lambda)
  n_parameters <- lenL-1 #number of lambda params if lambda0=0
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
