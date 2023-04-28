#' Oscillations of a MTD model
#'
#' Calculates the oscillations of a MTD model object.
#'
#' @param MTD an MTD object.
#'
#' @details The oscillations of a MTD model
#' (\eqn{\delta_j} for \eqn{j in \Lambda}), are the product of the weight \eqn{\lambda_j} times the maximum of the total variation distances between the distributions in the matrix p_j.
#' These values are important because they measure how much a relevant lag j influences the model.
#'
#' @return Returns the oscillations for each relevant lag of a MTD object.
#' @export oscillation
oscillation <- function(MTD){
  UseMethod("oscillation")
}

#' @export
oscillation.MTD <- function(MTD){
checkMTD(MTD)
lenA <- length(MTD$A) #number of rows/cols in each p_j
rows <- t(combn(lenA,2))
x <- MTD$lambdas[-1]*sapply(MTD$p_j,dTV_pj,rows)
names(x) <- paste0("-",MTD$Lambda)
class(x) <- "MTDoscillation"
x
}

#' @export
oscillation.default <- function(MTD){
  print("The implemented method can only calculate oscillations for MTD objects for now.")
}

#' @export
print.MTDoscillation <- function(x, ...){
 cat("Calculating \U03B4\U2096 = \U03BB\U2096*max{b,c in A: dTV[p\U2096(.|b),p\U2096(.|c)]}, \n")
 cat("for each k in Lambda: \n")
 cat("\n")
 #y <- as.numeric(names(x))
 x <- oscillation(MTD)
 print(rbind(x))
 return(invisible(x))
}
