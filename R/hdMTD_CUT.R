#' The Cut method.
#'
#' A function for inference in MTD Markov chains with CUT method. It applies CUT algorithm to estimate a relevant lag set \eqn{\Lambda} of a MTD model.
#'
#' @param X A MTD chain sample.
#' @param d An upper bound for the chains order.
#' @param S A subset of 1:d that contains the relevant lag set, if empty S=\eqn{1,2,\dots, d}.
#' @param A The states space. "A" only needs to be informed if X does not already contain all elements of "A".
#' @param alpha A parameter of the CUT algorithm. Alpha is constant of a threshold used in the CUT
#'step to determine if two distributions are different enough. The larger the alpha, the
#'greater the distance required between the distributions (to be considered different).
#' @param mu A parameter of the CUT algorithm. mu is also a component of the same threshold as alpha.
#' @param xi A parameter of the CUT algorithm. xi is also a component of the same threshold as alpha and mu.
#' @param warning If TRUE may return warnings.
#' @param ... Used to accommodate any extra arguments passed by the [hdMTD()] function.
#'
#'
#' @details The "Forward Stepwise and Cut" (FSC)is an algorithm for inference in
#' Mixture Transition Distribution (MTD) models.
#' It consists in the application of the "Forward Stepwise" (FS) step followed by the CUT algorithm.
#' This method was developed by [Ost and Takahashi](https://arxiv.org/abs/2202.08007) and is specially useful for high-order MTD Markov chains.
#' This function will only apply the CUT step of the algorithm.
#'
#' @return Returns a estimated set of relevant lags.
#' @export
#' @examples
#' X <- perfectSample(MTDmodel(Lambda=c(1,4),A=c(0,1)),N=1000)
#' hdMTD_CUT(X,4,alpha=0.02,mu=1,xi=0.4)
#' hdMTD_CUT(X,d=6,S=c(1,4,6),alpha=0.0065)
#'
hdMTD_CUT <- function(X, d, S=1:d, alpha=0.05, mu=1, xi=0.5, A=NULL, warning=FALSE,...){

  #Checking inputs
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
  }
  if( length(A)<=1   ||
      any(A%%1 !=0)   )stop("States space A must be a numeric vector with at least two integers.")
  if ( !all( unique(X) %in% A ) ) {
    stop("Check the states space, it must include all states that occur in the sample.")
  }
  while ( is.na(alpha) || !is.numeric(alpha) || alpha <= 0 ) {
    cat("The alpha value is not valid for CUT step. alpha should be a positive number.")
    alpha <- readline(prompt = "Please enter a valid alpha: ")
    alpha <- suppressWarnings(as.numeric(alpha))
  }
  while ( is.na(mu) || !is.numeric(mu) || mu <= 0 ) {
    cat("The mu value is not valid for CUT step. mu should be a positive number.")
    mu <- readline(prompt = "Please enter a valid mu: ")
    mu <- suppressWarnings(as.numeric(mu))
  }
  while ( is.na(xi) || !is.numeric(xi) || xi <= 0 ) {
    cat("The xi value is not valid for CUT step. xi should be a positive number.")
    xi <- readline(prompt = "Please enter a valid xi: ")
    xi <- suppressWarnings(as.numeric(xi))
  }
  #\.
  # Gathering inputs
  A <- sort(A)
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
  #\.

  dTV_txy <- numeric(lenS)
  for (z in 1:lenS) {
    j <- dec_S[z]
    Sminusj <- dec_S[ -which( dec_S == j ) ]

    Q <- matrix(0,ncol=nrowA_pairs,nrow = nrow_subx)
    R <- matrix(0,ncol = lenA, nrow = nrow_subx)
    for (k in 1:nrow_subx) { #runs in all x_S
      Q[k,] <- dTV_sample(S=Sminusj,j=j,lenA=lenA,base=b_Sja,A_pairs=A_pairs,x_S=subx[k,])
      R[k,] <- sx(S=Sminusj,freqTab=b_Sja,lenA=lenA,x_S=subx[k,],mu=mu,alpha=alpha,xi=xi)
    }
    colnames(Q) <- apply(A_pairs, 1, paste0, collapse="x")
    rownames(Q) <- apply(subx, 1, paste0, collapse="")
    assign(paste0("dtv_j",j),Q)
    colnames(R) <- A
    rownames(R) <- apply(subx, 1, paste0,collapse="")
    assign(paste0("sx_Sj",j),R)

    dTV_txy[z] <- max(Q-txy(R=R,A_pairs=A_pairs,A_pairsPos=A_pairsPos))
  }
  S <- sort(S,decreasing = TRUE)
  S <- S[which(dTV_txy>0)]
  S
}

#local auxiliary function
txy <- function(R,A_pairs,A_pairsPos){
  tn <- matrix(0,nrow=nrow(R),ncol=nrow(A_pairs))
  for (s in 1:nrow(A_pairs)) {
    tn[,s] <- apply(R[,A_pairsPos[s,]],1,sum)
  }
  colnames(tn)=apply(A_pairs, 1, paste0,collapse="x")
  rownames(tn)=rownames(R)
  tn
}
#\.


