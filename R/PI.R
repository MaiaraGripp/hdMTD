#' Estimated stationary distribution vector for every a_j and a sequence x_S
#'
#' @param S A set of lags.
#' @param base A table with frequencies of a Markov chain sample.
#' @param x_S A specific sequence indexed by S.
#' @param lenX The lenght of the sample vector.
#' @param d The order of the chain.
#'
#' @return Returns a vector of estimated stationary distributions of a sequence X_S
#' with the element in j varying in all the states space.
#' @importFrom dplyr %>%
#
#'
PI <- function(S,base,x_S,lenX,d){
  if (is.numeric(S)) {
   #S <- sort(S,decreasing = TRUE) #S needs to be decreasing for filtering
    filtr_S <- paste0("x",S)
    B <- base
    B$test <- apply( B %>% dplyr::select_at(filtr_S),1,is_xS,x_S) #true forall B$x_S == x_S
    C <- dplyr::filter(B,test == TRUE) #only rows with x_S remain
  }else{
    C <- base
  }
  inv <- matrix(C$Nx_Sj/(lenX-d),ncol = 1)
  colnames(inv) <- paste0(x_S,collapse = "")
  inv
}
##(pi(xa_Sj),pi(xb_Sj),pi(xc_Sj),...)
