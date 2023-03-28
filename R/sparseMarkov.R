#' A function for inference in sparse Markov chains
#'
#' @param X A markov chain sample.
#' @param d The order or an upper bound for the order.
#' @param method A method for estimation. The default method is the Forward Stepwise and Cut ("FSC"="FS"+"CUT"). Alternatively, the methods "FS" and "CUT" can be used separately.
#' @param l A stop point for the "FS" and "BIC" methods.
#' @param alpha "CUT" parameter. Defaulted to 0.05.
#' @param mu "CUT" parameter. Defaulted to 1.
#' @param xi "CUT" parameter if method="FSC" or "CUT". BIC constant if method="BIC".Defaulted to 0.5.
#'
#' @return The set of relevant lags.
#' @export
#'
sparseMarkov <- function(X,d,method="FSC",l=NULL, alpha=0.05, mu=1, xi=0.5){
  #l <- list(...)
  if(method!="CUT" & length(l)==0)stop("Parameter l can't be NULL for the choosen method.")
  if(method!="CUT"){
    if(!is.numeric(l) ||
       length(l)!=1 ||
       l%%1!=0 ||
       l<=0 ||
       l>=d )stop("l must be a integer number greater than 0 and equal to or lower than d.")
  }
  # Checking restrictions
  checkSample(X)
  if( !is.numeric(d) | d<2 | (d %% 1)!=0 ){
    stop("The order d must be an integer number greater than 2.")
  }
  method_list <- c("FSC","FS","CUT","BIC")
  if( !(toupper(method) %in% method_list) ) {
    stop("Unknown method.")
  }else{
      method <- toupper(method)
    }
  #\.
  # Gathering inputs from sample
  A <- sort(unique(X))
  lenA <- length(A)
  lenX <- length(X)

  #A_pairs <- t(combn(A,2))
  #nrowA_pairs <- nrow(A_pairs)

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

    S <- sparseMarkov_FS(Xm,l=l,A=A,d=d)
    S <- sparseMarkov_CUT(Xn, A=A, d=d, alpha=alpha, mu=mu, xi=xi, S=S)
  }
  if ( method == "FS" ) {
    S <- sparseMarkov_FS(X=X,l=l,A=A,d=d)
  }
  if ( method == "CUT" ) {
    S <- sparseMarkov_CUT(X=X, A=A, d=d, alpha=alpha, mu=mu, xi=xi)
  }
  if (method == "BIC" ){
    S <- sparseMarkov_BIC(X=X,d=d,l=l,c=xi)
  }
  return(S)
}


