#' Calculates the number of parameters in a MTD model given the Lambda set.
#'
#' @param Lambda the relevant lag set.
#' @param A the states space.
#' @param single_matrix If TRUE all p_j matrix are equal.
#' @param indep_part If FALSE independent distribution is set to zero.
#'
#' @return The number of parameters of a MTD model.
n_parameters <- function(Lambda, A, single_matrix = FALSE, indep_part = TRUE){
  #single_matrix=TRUE => p_i(.|a)=p_j(.|a) for all i \neq j and a\in Alphabet.
  #indep_part=FALSE => no term independent from past.
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
