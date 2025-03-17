#' The total variation distance between distributions
#'
#' Calculates the total variation distance between distributions conditioned
#' in a given past sequence.
#'
#' @details This function computes the total variation distance between distributions
#' found in \code{base}, which is expected to be the output of the function [freqTab()].
#' Therefore, \code{base} must follow a specific structure (e.g., column names must
#' match, and a column named qax_Sj, containing transition distributions, must be
#' present). For more details on the output structure of [freqTab()], refer to its
#' documentation..
#'
#' If you provide the state space \code{A}, the function calculates:
#' \code{lenA <- length(A)} and \code{A_pairs <- t(utils::combn(A, 2))}.
#' Alternatively, you can input \code{lenA} and \code{A_pairs} directly and let
#' \code{A <- NULL}, which is useful in loops to improve efficiency.
#'
#' @param S A numeric vector of positive integers (or \code{NULL}) representing
#' a set of past lags. The distributions from which this function will calculate the total
#' variation distance are conditioned on a fixed sequence indexed by \code{S} ( the user must
#'  also input the sequence through the argument \code{x_S}).
#' @param j A positive integer representing a lag in the \eqn{complement} of \code{S}.
#'  The symbols indexed by \code{j} vary along the state space \code{A}, altering the distribution
#'  through this single lag, and the size of this change is what this function seeks to measure.
#' @param A  A vector of unique positive integers (state space) with at least two elements. \code{A}
#' represents the state space. You may leave \code{A=NULL} (default) if you provide the function
#' with the arguments \code{lenA} and \code{A_pairs} (see *Details* below).
#' @param base A data frame with sequences of elements from \code{A} and their transition
#' probabilities. \code{base} is meant to be an output from function [freqTab()],
#' and must be structured as such. The data frame must contain all required transitions
#' conditioned on \code{x_S} (i.e. \code{length(A)^2} rows with sequence \code{x_S}).
#' See *Details* section for further information.
#' @param lenA An integer \code{>= 2}, representing \code{length(A)}. Required if \code{A} is not
#' provided.
#' @param A_pairs A two-column matrix with all unique pairs of elements from \code{A}.
#' Required if \code{A} is not provided.
#' @param x_S A vector of length \code{length(S)} or \code{NULL}. If \code{S==NULL}, \code{x_S} will
#' be set to \code{NULL}. \code{x_S} represents a sequence of symbols from \code{A} indexed by
#' \code{S}. This sequence remains constant across the conditional distributions to be compared,
#' representing the fixed configuration of the past.
#'
#' @return Returns a vector of total variation distances, where each entry corresponds to the
#' distance between a pair of distributions conditioned on the same fixed past \code{x_S},
#' differing only in the symbol indexed by \code{j}, which varies across all distinct pairs
#' of elements in \code{A}.The output has length equal to the number of unique pairs
#' in \code{A_pairs}.
#'
#' @import dplyr
#' @import purrr
#' @export
#' @examples
#' #creating base argument through freqTab function.
#' pbase <- freqTab(S=c(1,4),j=2,A=c(1,2,3),countsTab = countsTab(testChains[,2],d=5))
#' dTV_sample(S=c(1,2),j=4,A=c(1,2,3),base=pbase,x_S=c(2,3))
#' pbase <- freqTab(S=NULL,j=1,A=c(1,2,3),countsTab = countsTab(testChains[,2],d=5))
#' dTV_sample(S=NULL,j=1,A=c(1,2,3),base=pbase)
#'
dTV_sample <- function(S, j, A = NULL, base, lenA = NULL, A_pairs = NULL, x_S){

  check_dTVsample_inputs(S, j, A, base, lenA, A_pairs, x_S)

  # If A is provided, compute lenA and A_pairs
  if( !is.null(A) ){
    A <- sort(A)
    lenA <- length(A)
    A_pairs <- t(utils::combn(A, 2))
  }
  nrowA_pairs <- nrow(A_pairs)

  S <- sort(S,decreasing = TRUE)
  lenS <- length(S)
  # If S is empty, set x_S to NULL
  if(lenS == 0) x_S <- NULL


  # Filters base according to x_S if S is not NULL
  if (!is.null(S)) {
    filtr_S <- paste0("x", S)
    base <- base %>%
      mutate(match = purrr::pmap_lgl(pick(all_of(filtr_S)), ~ all(c(...) == x_S))) %>%
      filter(match) %>%
      select(-match)
  }

  # Check if base contains the expected number of rows
  if(nrow(base)!=lenA^2) {
    stop("The base does not contain all required distributions conditioned on x_S.")
  }

  # Compute total variation distances
  filtr_j <- paste0("x",j)
  disTV <- matrix(0,ncol=nrowA_pairs)

  for (i in seq_len(nrowA_pairs)) {
    # Filter rows where column `xj` matches one of the elements in the current pair from A_pairs
    D <- base %>% dplyr::filter(.data[[filtr_j]] %in% A_pairs[i, ])
    # Compute total variation distance between the two distributions defined by the pair
    disTV[i] <- sum(abs(D$qax_Sj[1:lenA]-D$qax_Sj[(lenA+1):(2*lenA)]))/2
  }
  colnames(disTV) <- apply(A_pairs, 1, paste0, collapse="x")
  rownames(disTV) <- paste0(x_S,collapse = "")
  disTV
}


