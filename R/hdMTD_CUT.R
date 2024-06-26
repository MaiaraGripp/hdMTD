#' The Cut method for inference in MTD models
#'
#' A function that estimates the set of relevant lags of an MTD model using the CUT method.
#'
#' @param X A vector or single-column data frame containing a chain sample.
#' @param d A positive integer representing an upper bound for the chain order.
#' @param S A numeric vector of positive integers from which this function will select
#' a set of relevant lags. Typically, \code{S} is a subset of \code{1:d}. If \code{S}
#' is not provided, by default \code{S=1:d}.
#' @param A A vector with positive integers representing the state space. If not informed,
#' this function will set \code{A=unique(X)}.
#' @param alpha A positive real number, \code{alpha}, appears in a threshold used in the CUT
#'  step to determine if two distributions are different enough. The larger the \code{alpha},
#'  he greater the distance required to consider that there is a difference between a set
#' of distributions.
#' @param mu A positive real number between 0 and 3. \code{mu}is also a component of the same
#' threshold as \code{alpha}.
#' @param xi A positive real number, \code{xi} is also a component of the same threshold as
#'  \code{alpha} and \code{mu}.
#' @param warning Logical. If \code{TRUE}, informs the user if the state space was set as
#' \code{A=unique(X)}.
#' @param ... Other parameters. This is used to accommodate any additional argument passed
#' to [hdMTD_CUT()] through the [hdMTD()] function.
#'
#' @details The "Forward Stepwise and Cut" (FSC) is an algorithm for inference in
#' Mixture Transition Distribution (MTD) models. It consists
#' in the application of the "Forward Stepwise" (FS) step followed by the CUT algorithm.
#' This method and its steps where developed by [Ost and Takahashi](https://arxiv.org/abs/2202.08007)
#' and are specially useful for inference in high-order MTD Markov chains. This specific function
#' will only apply the CUT step of the algorithm and return an estimated relevant lag set.
#'
#' @references
#' Ost, G. and Takahashi, D. (2022), ‘Sparse markov models for high-dimensional inference’.
#' [ arXiv:2202.08007](https://arxiv.org/abs/2202.08007)
#'
#' @return Returns a set of relevant lags estimated using the CUT algorithm.
#' @export
#' @examples
#' X <- testChains[,3]
#' hdMTD_CUT(X,4,alpha=0.02,mu=1,xi=0.4)
#' hdMTD_CUT(X,d=6,S=c(1,4,6),alpha=0.0065)
#'
hdMTD_CUT <- function(X, d, S=1:d, alpha=0.05, mu=1, xi=0.5, A=NULL, warning=FALSE,...){

## Checking S, d, X, A
  if(length(S) < 2  ||
     !is.numeric(S) ||
     any( S%%1 != 0) ){stop("Parameter 'S' must be a numeric vector containing at least 2 integer values.")}
  if( !is.numeric(d) || d<2 || (d %% 1)!=0 || d<max(S)){
    stop("The 'd' parameter must be an integer greater than 2 and at least as large as the greatest value in 'S'.")
  }
  X <- checkSample(X)
  if(length(A)==0){
    if(warning==TRUE){
      warning("States space A not informed. Code will set A <- sort(unique(X)).")
    }
    A <- sort(unique(X))
  }else{A <- sort(A)}
  if( length(A)<=1   ||
      any(A%%1 !=0)   )stop("States space A must be a numeric vector with at least two integers.")
  if ( !all( unique(X) %in% A ) ) {
    stop("Check the states space, it must include all states that occur in the sample.")
  }
## ask to change constantes alpha, mu, xi if needed
  while ( is.na(alpha) || !is.numeric(alpha) || alpha <= 0 ) {
    cat("The alpha value is not valid for CUT step. alpha should be a positive number.")
    alpha <- readline(prompt = "Please enter a valid alpha: ")
    alpha <- suppressWarnings(as.numeric(alpha))
  }
  while ( is.na(mu) || !is.numeric(mu) || mu <= 0 || mu >=3) {
    cat("The mu value is not valid for CUT step. mu should be a positive number less than or equal to 3.")
    mu <- readline(prompt = "Please enter a valid mu: ")
    mu <- suppressWarnings(as.numeric(mu))
  }
  while ( is.na(xi) || !is.numeric(xi) || xi <= 0 ) {
    cat("The xi value is not valid for CUT step. xi should be a positive number.")
    xi <- readline(prompt = "Please enter a valid xi: ")
    xi <- suppressWarnings(as.numeric(xi))
  }

## Gathering inputs
  lenA <- length(A)
  lenS <- length(S)
  subx <- as.matrix(expand.grid(rep(list(A),lenS-1))[,(lenS-1):1],ncol=lenS-1)
  nrow_subx <- nrow(subx)
  dec_S <- sort(S,decreasing = TRUE)
  A_pairs <- t(utils::combn(A,2))
  A_pairsPos <- t(utils::combn(1:lenA,2))
  nrowA_pairs <- nrow(A_pairs)
  base <- countsTab(X=X,d=d)
  b_Sja <- freqTab(S=S,A=A,countsTab=base)


  dTV_txy <- numeric(lenS)
  for (z in 1:lenS) { #It runs through all lags in S, singling out one in each loop.
    j <- dec_S[z]
    Sminusj <- dec_S[ -which( dec_S == j ) ]
    Q <- matrix(0,ncol=nrowA_pairs,nrow = nrow_subx)
    R <- matrix(0,ncol = lenA, nrow = nrow_subx)
    for (k in 1:nrow_subx) { #runs through all sequences x of size |S|-1 ( i.e. x_{S\j})
      Q[k,] <- dTV_sample(S=Sminusj,j=j,lenA=lenA,base=b_Sja,A_pairs=A_pairs,x_S=subx[k,])
      R[k,] <- sx(S=Sminusj,freqTab=b_Sja,lenA=lenA,x_S=subx[k,],mu=mu,alpha=alpha,xi=xi)
      #see auxiliary function sx() below.
    }
    # Q: a matrix for all total var distances when changing the element in past j
    # R: a matrix with thresholds for the entries in Q
    colnames(Q) <- apply(A_pairs, 1, paste0, collapse="x")
    rownames(Q) <- apply(subx, 1, paste0, collapse="")
    assign(paste0("dtv_j",j),Q)
    colnames(R) <- A
    rownames(R) <- apply(subx, 1, paste0,collapse="")
    assign(paste0("sx_Sj",j),R)

    dTV_txy[z] <- max(Q-txy(R=R,A_pairs=A_pairs,A_pairsPos=A_pairsPos))
    #see auxiliary function txy() below
  }
  S <- sort(S,decreasing = TRUE)
  S <- S[which(dTV_txy>0)]
  S
}

############################################################
# Local auxiliary functions for calculating the CUT threshold:
############################################################
  # sx(): Returns a threshold with respect to some sequence x_S
sx <- function(S,freqTab,lenA,x_S,mu,alpha,xi){
  filtr_S <- paste0("x",S)
  B <- freqTab
  B$test <- apply(B %>%
                    dplyr::select_at(filtr_S),1,is_xS,x_S)
  C <- dplyr::filter(B,test==TRUE)
  Nx <- C$Nx_Sj[seq(1,nrow(C),by=lenA)]
  result <- numeric(lenA)
  for (k in 1:lenA) {
    sums <- 0
    for (i in 1:lenA) {
      sums = sums +
        sqrt(( C$qax_Sj[i+(k-1)*lenA] +
                 alpha/Nx[k])*mu/(2*mu-exp(mu)+1))
    }
    result[k]=sqrt(0.5*alpha*(1+xi)/Nx[k])*sums+alpha*lenA/(6*Nx[k])
  }
  result
}

# tyx(): txy = sx + sy is the threshold that determines if the difference
# between a distribution conditioned on x_S and another conditioned on y_S is indeed relevant.
txy <- function(R,A_pairs,A_pairsPos){
  tn <- matrix(0,nrow=nrow(R),ncol=nrow(A_pairs))
  for (s in 1:nrow(A_pairs)) {
    tn[,s] <- apply(R[,A_pairsPos[s,]],1,sum)
  }
  colnames(tn)=apply(A_pairs, 1, paste0,collapse="x")
  rownames(tn)=rownames(R)
  tn
}



