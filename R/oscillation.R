#' Calculates the oscillation of a MTD model
#'
#' @param MTD an MTD object.
#'
#' @return Returns the oscillations of a MTD object, that is, for each
#' relevant lag j, returns lambda_j times the maximum of the total variation
#'  distance between the distributions in matrix p_j.
#' @export
oscillation <- function(MTD){
checkMTD(MTD)
lenA <- length(MTD$A) #number of rows/cols in each p_j
rows <- t(combn(lenA,2))
x <- MTD$lambdas[-1]*sapply(MTD$p_j,dTV_pj,rows)
names(x) <- paste0("-",Lambda)
x
}
