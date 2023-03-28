#' Estimates a relevant lag set Lambda with BIC
#'
#' @param X A Markov chain.
#' @param d A upper threshold for the chains order.
#' @param l A upper threshold for the number of elements in the relevant lag set.
#' @param c The BIC constant. Defaulted to 1/2. Smaller c `(near 0)` reduces the impact of overparameterization.
#'
#' @return Return a estimator of the relevant lag set using BIC.
#' @export
sparseMarkov_BIC <- function(X,d,l,c=1/2){
  pML <- BIC_l(X=X,d=d,l=l,c=c)
    smallest <- names(unlist(pML))[order(unlist(pML))][1]
    estLambda <- sapply(sapply(pML, orderedNames), dplyr::first )
    estLambda <- append(estLambda,smallest)
    names(estLambda) <- c(paste0("l=",1:l),"smallest")
  estLambda
}

orderedNames <- function(x){
  names(x[order(x)])
}
