#' A tibble showing the frequency of sample sequences
#'
#' Creates a tibble that displays the frequencies of sequences within the sample.
#'
#' @param S A numeric vector of positive integers or \code{NULL}. Represents a set
#'  of past lags that must be present within the columns of the \code{freqTab} argument. Typically,
#'  \code{S} is a subset of \code{1:d}.
#' @param j An integer or \code{NULL}. Typically represents a lag in the
#'  complement of \code{S} (in \code{1:d}). Both \code{S} and \code{j} cannot be \code{NULL}
#'  at the same time.
#' @param freqTab A tibble generated by the [freqTab()] function, must contain columns referring to lags
#' \code{S} and \code{j}.
#' @param lenX A positive integer representing the number of elements in the sample.
#' @param d A positive integer representing an upper bound for the chain order.
#'
#' @return A tibble with sequences indexed by the lags in \code{S} and \code{j}, their frequencies,
#'  and the maximum likelihood estimator (MLE) of the transition probabilities.
#' @importFrom dplyr %>%
#'
freqTabSj <- function(S,j,freqTab,lenX,d){
  Sj <- sort(c(S,j),decreasing = TRUE)
  if ( is.numeric(Sj) ) {
    freqTabSj <- freqTab %>%
                  dplyr::group_by_at(paste0("x",Sj)) %>%
                  dplyr::summarise(Nx_Sj=sum(Nxa_Sj), .groups="drop")
    freqTabSj
  }else{
    matrix(c(0,lenX-d),ncol=2)
  }
}
