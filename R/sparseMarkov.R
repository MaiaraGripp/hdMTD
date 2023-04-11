#' A function for inference in MTD Markov chains
#'
#' This function returns an estimation of the relevant lag set \eqn{\Lambda}
#' of a MTD model. The method for this estimation must be set in the `method` parameter which is defaulted to "FSC".
#' See [Ost and Takahashi](https://arxiv.org/abs/2202.08007).
#'
#' @param X A MTD Markov chain sample.
#' @param A The states space.
#' @param d The order or an upper bound for the order.
#' @param method One of the following methods for estimation:
#'   ## "FSC":
#'      The default method "Forward Stepwise and Cut". See [Ost and Takahashi](https://arxiv.org/abs/2202.08007).
#'
#'   ## "FS":
#'      Applies only the first part of the "FSC" method, i.e. the "Forward Stepwise" part.
#'
#'   ## "CUT":
#'      Applies only the second part of the "FSC" method, i.e. the "CUT" part.
#'
#'   ## "BIC":
#'      A method that selects a relevant lag set using Bayesian Information Criterion
#'
#'
#' @param l A stop point for the "FS" and "BIC" methods.
#' @param alpha "CUT" parameter. Defaulted to 0.05.
#' @param mu "CUT" parameter. Defaulted to 1.
#' @param xi "CUT" parameter if method="FSC" or "CUT". BIC constant if method="BIC".Defaulted to 0.5.
#' @param warning If TRUE may return warnings.
#'
#' @return The estimated set of relevant lags.
#' @export
#'
sparseMarkov <- function(X,A=NULL,d,method="FSC",l=NULL, alpha=0.05, mu=1, xi=0.5, warning=FALSE){

  #check methods
  if(method!="CUT" & length(l)==0)stop("Parameter l can't be NULL for the choosen method.")
  if(method!="CUT"){
    if(!is.numeric(l) ||
       length(l)!=1 ||
       l%%1!=0 ||
       l<=0 ||
       l>=d )stop("l must be a integer number greater than 0 and equal to or lower than d.")
  }
  method_list <- c("FSC","FS","CUT","BIC")
  if( !(toupper(method) %in% method_list) ) {
    stop("Unknown method.")
  }else{
      method <- toupper(method)
    }
  #\.
  # Gathering inputs from sample
  lenX <- length(X)
  #\.
  if ( method == "FSC" ) {
    # Checing restrictions for FSC
    if ( lenX <= 2*(d+1)) {
      stop("The FSC method splits data in two,
           therefore the sample size must be greater than 2*(d+1).")
    }
    #\.
    # Gathering inputs for FSC from sample
    m <- lenX %/% 2
    Xm <- X[1:m]
    Xn <- X[(m+1):lenX]
    n <- length(Xn)

    S <- sparseMarkov_FS(Xm,l=l,A=A,d=d,warning=warning)
    S <- sparseMarkov_CUT(Xn, A=A, d=d, alpha=alpha, mu=mu, xi=xi, S=S)
  }
  if ( method == "FS" ) {
    S <- sparseMarkov_FS(X=X,l=l,A=A,d=d,warning=warning)
  }
  if ( method == "CUT" ) {
    S <- sparseMarkov_CUT(X=X, A=A, d=d, alpha=alpha, mu=mu, xi=xi)
  }
  if (method == "BIC" ){
    S <- sparseMarkov_BIC(X=X,A=A,d=d,l=l,xi=xi)[l]
    S <- as.numeric(unlist(strsplit(S, ",")))
  }
  return(S)
}

