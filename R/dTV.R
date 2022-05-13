#' Total variation distance between distributions given x_S
#'
#' @param S A set of lags.
#' @param j A specific lag \eqn{j \in S^c}.
#' @param lenA The length of the states space A
#' @param base  A set with the estimated transition probabilities q(a|xb_Sj)
#' @param A_pairs A list with all pairs with elements of A
#' @param x_S A specific sequence indexed by S
#'
#' @return Returns a vector with the total variation distance between estimated distributions
#' given a x_S
#' @importFrom dplyr %>%
#' @details (dTv_xS p(.|a_j),p(.|b_j),dTv_xS p(.|a_j),p(.|c_j) ...)
#'
dTV <- function(S,j,lenA,base,A_pairs,x_S){

  nrowA_pairs <- nrow(A_pairs)
  if ( is.numeric(S) ) {
    S <- sort(S,decreasing = TRUE) #S needs to be decreasing for filtering
    filtr_S <- paste0("x",S)
    B <- base
    B$test <- apply(B %>% dplyr::select_at(filtr_S),1,is_xS,x_S)
    C <- dplyr::filter(B,test==TRUE)
  }else{
    C <- base
  }

  filtr_j <- paste0("x",j)
  disTV <- matrix(0,ncol=nrowA_pairs)
  for (i in 1:nrowA_pairs) {
    D <- C %>% dplyr::filter_at(filtr_j, ~. %in% A_pairs[i,])
    disTV[i] <- sum(abs(D$qax_Sj[1:lenA]-D$qax_Sj[(lenA+1):(2*lenA)]))/2
# sum_{\in a} |q(a|xb_Sj)-q(a|xc_Sj)|/2
  }
  colnames(disTV) <- apply(A_pairs, 1, paste0, collapse="x")
  rownames(disTV) <- paste0(x_S,collapse = "")
  disTV #vector size nrowA_pairs
}
