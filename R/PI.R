#' Estimated stationary distribution.
#'
#' Estimated stationary distribution vector for every a_j and x_S
#'
#' @param S A subset of 1:d
#' @param freqTabSj A table with frequencies of sequences in a Markov chain sample.
#' @param x_S A specific sequence indexed by S.
#' @param lenX The lenght of the sample vector.
#' @param d The order of the chain.
#'
#' @return Returns a vector of estimated stationary distributions of a sequence X_S
#' with the element in j varying in all the states space.
#' @importFrom dplyr %>%
#
#'
PI <- function(S,freqTabSj,x_S,lenX,d){
  if (is.numeric(S)) {
    filtr_S <- paste0("x",S)
    B <- freqTabSj
    B$test <- apply( B %>% dplyr::select_at(filtr_S),1,is_xS,x_S) #true forall B$x_S == x_S
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
