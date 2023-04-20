#' A function for inference in MTD Markov chains with FSC method.
#'
#' This function estimates the relevant lag set \eqn{\Lambda}
#' of a MTD model through the FSC algorithm.
#'
#' @param X A MTD Markov chain sample.
#' @param A The States space.
#' @param d An upper threshold for the chains order.
#' @param l Stop point for FS algorithm.
#' @param alpha A parameter of CUT.
#' @param mu A parameter of CUT.
#' @param xi A parameter of CUT.
#' @param warning If TRUE may return warnings.
#'
#' @details The "Forward Stepwise and Cut" (FSC)is an algorithm for inference in
#' Mixture Transition Ditribution (MTD) models.
#' It consists in the application of the "Forward Stepwise" (FS) step followed by the CUT algorithm.
#' This method was developed by [Ost and Takahashi](https://arxiv.org/abs/2202.08007) and is specially useful for High order MTD Markov chains.
#'
#' @return Returns a estimated set of relevant lags.
#' @export
sparseMarkov_FSC <- function(X,A=NULL,d,l=NULL, alpha=0.05, mu=1, xi=0.5, warning=FALSE){
  # Checking inputs
    # Sample
  X <- checkSample(X)
    # A
  if(length(A)==0){
    if(warning==TRUE){
      warning("States space A not informed. Code will set A <- sort(unique(X)).")
    }
    A <- unique(X)
  }
  if( !is.numeric(A) ||
      length(A)<=1   ||
      length(dim(A))!=0 )stop("States space A must be a numeric vector with at least two values.")
  if ( !all( unique(X) %in% A ) ) {
    stop("Check the states space, it must have all states that occur in the sample.")
  }
    # d
  if( !is.numeric(d) || d<2 || (d %% 1)!=0 ){
    stop("The order d must be an integer number greater than 2.")
  }
    # l
  while ( is.na(l) || !is.numeric(l) || l%%1 != 0 || l>d || l>length(S) ) {
    cat("l value is not valid. l should be a positive integer lower or equal to d or the number of elements in S.")
    l <- readline(prompt = "Please enter a valid l : ")
    l <- suppressWarnings(as.numeric(l))
  }
    # alpha
  while ( is.na(alpha) || !is.numeric(alpha) || alpha <= 0 ) {
    cat("alpha value is not valid for CUT step. alpha should be a positive number.")
    alpha <- readline(prompt = "Please enter a valid alpha: ")
    alpha <- suppressWarnings(as.numeric(alpha))
  }
    # mu
  while ( is.na(mu) || !is.numeric(mu) || mu <= 0 ) {
    cat("mu value is not valid for CUT step. mu should be a positive number.")
    mu <- readline(prompt = "Please enter a valid mu: ")
    mu <- suppressWarnings(as.numeric(mu))
  }
    # xi
  while ( is.na(xi) || !is.numeric(xi) || xi <= 0 ) {
    cat("xi value is not valid for CUT step. xi should be a positive number.")
    xi <- readline(prompt = "Please enter a valid xi: ")
    xi <- suppressWarnings(as.numeric(xi))
  }
  if(!logical(warning)){stop("warning must be TRUE or FALSE.")}


  lenX <- length(X)
  if ( lenX <= 2*(d+1)) {
    stop("The FSC method splits data in two,
           therefore the sample size must be greater than 2*(d+1).")
  }
  m <- lenX %/% 2
  Xm <- X[1:m]
  Xn <- X[(m+1):lenX]
  n <- length(Xn)
  S <- sparseMarkov_FS(Xm,A=A,d=d,l=l,warning=warning)
  S <- sparseMarkov_CUT(Xn,A=A,d=d,S=S,alpha=alpha,mu=mu,xi=xi)
  s
}
