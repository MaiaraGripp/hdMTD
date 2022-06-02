#' Checks which row in data frame is equal to a vector
#'
#' @param x A data frame whose rows are sequencies indexed by S.
#' @param y A specific x_S sequence.
#'
#' @return A vector of true and false if a row in x is exactly equal to x_S.
is_xS <- function(x,y) {
  return( all( x == y ) ) #I used return() so R would let me use Roxygen
  }
