#' Total variation distance for pj distribution.
#'
#' Calculates the total variation distance given a stochastic matrix pj.
#'
#' @param pj A stochastic matrix.
#' @param rows A matrix with two columns, the pair of numbers in each row represents the rows of pj that
#'  will be compared, that is, the distributions in pj to be compared.
#' @param which_max If TRUE, the function will also return the rows of the pj matrices with the distributions
#' that have the greatest total variation distance.
#'
#' @return Returns the total variation distance of a set of distributions in pj.
#'
dTV_pj <- function(pj,rows,which_max=FALSE){
  dV <- numeric(nrow(rows))
  for (i in 1:nrow(rows)) {
    dV[i] <- sum(abs(pj[rows[i,1],]-pj[rows[i,2],]))/2
  }
  if(which_max==FALSE){
    max(dV)
  }else{
    pos <- which(dV==max(dV))
    list(max(dV),paste0(rows[pos,],collapse="x"))
  }
}
