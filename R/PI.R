#' Estimated stationary distribution.
#'
#' An estimated stationary distribution vector for any given sequence.
#'
#' @param S A numeric vector of positive integers. Typically, \code{S} is a subset of \code{1:d},
#' and represents a set of relevant lags.
#' @param freqTabSj A tibble with frequencies of sequences in a chain sample. Typically
#' \code{freqTabSj} is an output from [freqTabSj()].
#' @param x_S A vector of length \code{length(S)} or \code{NULL}. If \code{S=NULL}, \code{x_S} will
#'  be set to \code{NULL}. \code{x_S} represents a sequence of symbols from \code{A} indexed by
#'  \code{S}. This sequence remains constant across the stationary distributions to be calculated,
#'  representing a fixed configuration of the past, while the observations in any other lag displayed
#'  in \code{freqTabSj} that are not in \code{S} may vary.
#' @param lenX An integer positive number with the length of the sample vector.
#' @param d A positive integer representing an upper bound for the chain order.
#'
#' @return Returns a vector of estimated stationary distributions for a set of sequences. However, these
#' sequences might have a fixed part of elements \code{x_S} common among them and other elements varying
#' across the entire state space.
#' @importFrom dplyr %>%
#
#'
PI <- function(S,freqTabSj,x_S,lenX,d){
  if (is.numeric(S)) {
    filtr_S <- paste0("x",S)
    B <- freqTabSj
    B$test <- apply( B %>% dplyr::select_at(filtr_S),1,is_xS,x_S) #true for all B$x_S == x_S
    C <- dplyr::filter(B,test == TRUE) #only rows with x_S remain
  }else{
    C <- freqTabSj
  }
  inv <- matrix(C$Nx_Sj/(lenX-d),ncol = 1)
  colnames(inv) <- paste0(x_S,collapse = "")
  inv
}

is_xS <- function(x,y) {
  return( all( x == y ) )
}
