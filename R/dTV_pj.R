#' Total distance variation for p_j distribution.
#'
#' Calculates the total distance variation given a stochastic matrix p_j
#'
#' @param p_j A stochastic matrix.
#' @param rows A matrix with two columns, the pair of numbers in each row represents the rows of p_j to be compared,
#'  i.e. the distributions in p_j to be compared.
#'
#' @return Returns de total variation distance of a set of distributions in p_j.
#'
dTV_pj <- function(p_j,rows){
  dV <- numeric(nrow(rows))
  for (i in 1:nrow(rows)) {
    dV[i] <- sum(abs(p_j[rows[i,1],]-p_j[rows[i,2],]))/2
  }
  max(dV)
}
