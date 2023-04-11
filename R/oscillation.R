#' Oscillations of a MTD model
#'
#' Calculates the oscillations of a MTD model
#'
#'
#' @param MTD an MTD object.
#'
#' @details The oscillations of a MTD model
#' (\eqn{\delta_j} for \eqn{j in \Lambda}), are the product of the weight \eqn{\lambda_j} times the maximum of the total variation distances between the distributions in the matrix p_j.
#' These values are important because they measure how much a relevant lag j influences the model.
#'
#' @return Returns the oscillations for each relevant lag of a MTD object.
#' @export
oscillation <- function(MTD){
checkMTD(MTD)
lenA <- length(MTD$A) #number of rows/cols in each p_j
rows <- t(combn(lenA,2))
x <- MTD$lambdas[-1]*sapply(MTD$p_j,dTV_pj,rows)
names(x) <- paste0("-",MTD$Lambda)
x
}
