#' Total variation distance for p_j distribution.
#'
#' Calculates the total variation distance given a stochastic matrix p_j.
#'
#' @param p_j A stochastic matrix.
#' @param rows A matrix with two columns, the pair of numbers in each row represents the rows of p_j to be compared,
#'  i.e. the distributions in p_j to be compared.
#' @param which_max If TRUE, the function will also return the rows of the p_j matrices with the distributions
#' that have the greatest total variation distance.
#'
#' @return Returns the total variation distance of a set of distributions in p_j.
#'
dTV_pj <- function(p_j,rows,which_max=FALSE){
  dV <- numeric(nrow(rows))
  for (i in 1:nrow(rows)) {
    dV[i] <- sum(abs(p_j[rows[i,1],]-p_j[rows[i,2],]))/2
  }
  if(which_max==FALSE){
    max(dV)
  }else{
    pos <- which(dV==max(dV))
    list(max(dV),paste0(rows[pos,],collapse="x"))
  }
}
