#' The maximum total variation distance
#'
#' Computes the total variation distance between any two rows of a stochastic matrix.
#' Can analyze several pairs of rows at once, and if so, will return the greatest
#' total variation distance calculated among the pairs.
#'
#' @param pj A stochastic matrix.
#' @param rows A numeric matrix with two columns. Each row must contain a pair of
#' numbers representing the indices of the rows in \code{pj} that will be compared,
#' i.e., the distributions.
#' @param which_max Logical. If TRUE, the function will also return the pair in
#' \code{rows} that represents the row indices of the matrix \code{pj} with the
#' greatest total variation distance.
#'
#' @return Returns the maximum total variation distance of a set of distributions in
#' \code{pj}. If \code{which_max=TRUE}, it returns a list containing the maximum total
#' variation distance along with the pair of row indices in \code{pj} representing the
#' distributions with the greatest distance.
#'
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
