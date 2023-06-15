#' Oscillations of a MTD model
#'
#' Calculates the oscillations of a MTD model object.
#'
#' @param x Must be a MTD object.
#'
#' @details The oscillations of a MTD model
#' (\eqn{\delta_k} for \eqn{k in \Lambda}), are the product of the weight \eqn{\lambda_k} times the maximum of the total variation distances between the distributions in the matrix p_k.
#' These values are important because they measure how much a relevant lag k influences the model.
#'
#' @return Returns the oscillations for each relevant lag of a MTD object.
#' @export oscillation
#' @examples oscillation( MTDmodel(Lambda=c(1,4),A=c(2,3) ) )
#' oscillation( MTDmodel(Lambda=c(1,4),A=c(2,3),w0=0.01 ) )
#' oscillation( MTDmodel(Lambda=c(1,4),A=c(2,3),w0=0.01,w_j=c(0.49,0.5)) )
#' oscillation( MTDmodel(Lambda=c(1,4),A=c(2,3),w0=0.01,w_j=c(0.49,0.5),
#' p_j=list(matrix(c(0.1,0.9,0.9,0.1),ncol=2)), single_matrix=TRUE))
#'
oscillation <- function(x){ UseMethod("oscillation") }

#' @export
oscillation.MTD <- function(x){
  checkMTD(x)
  lenA <- length(x$A) #number of rows/cols in each p_k
  rows <- t(combn(lenA,2))
  y <- x$lambdas[-1]*sapply(x$p_j,dTV_pj,rows)
  names(y) <- paste0("-",x$Lambda)
  class(y) <- "MTDoscillation"
  y
}

#' @export
oscillation.default <- function(x){
  print("The implemented method can only calculate oscillations for MTD objects for now.")
}


#' @export
print.MTDoscillation <- function(x, ...){
 class(x) <- NULL
 cat("Calculating \U03B4\U2096 = \U03BB\U2096*max{b,c in A: dTV[p\U2096(.|b),p\U2096(.|c)]}, \n")
 cat("for each k in Lambda: \n")
 cat("\n")
 print(x)
 return(invisible(x))
}

