#' A function for inference in MTD Markov chains with CUT method.
#'
#' Applies Cut algorithm to estimate a relevant lag set \eqn{\Lambda} of a MTD model.
#'
#' @param X A mixture transition distribution (MTD) chain sample.
#' @param A The states space.
#' @param d An upper threshold for the chains order..
#' @param S A set of relevant lags, if empty S=\eqn{1,2,\dots, d}.
#' @param alpha A parameter of CUT.
#' @param mu A parameter of CUT.
#' @param xi A parameter of CUT.
#' @param warning If TRUE may return warnings.
#'
#'
#' @details The "Forward Stepwise and Cut" (FSC)is an algorithm for inference in
#' Mixture Transition Ditribution (MTD) models.
#' It consists in the application of the "Forward Stepwise" (FS) step followed by the CUT algorithm.
#' This method was developed by [Ost and Takahashi](https://arxiv.org/abs/2202.08007) and is specially useful for High order MTD Markov chains.
#' This function will apply only the CUT step of the algorithm.
#'
#' @return Returns a estimated set of relevant lags.
#' @export
#'
sparseMarkov_CUT <- function(X, A, d, S=1:d, alpha=0.05, mu=1, xi=0.5, warning=FALSE){

  #Checking inputs
  if(length(S) < 2  ||
     !is.numeric(S) ||
     any( S%%1 != 0) ){stop("S must be a vector of at least 2 integer numbers.")}
  if( !is.numeric(d) || d<2 || (d %% 1)!=0 || d<max(S)){
    stop("The order d must be an integer number greater than 2 or the greatest element in S.")
  }
  checkSample(X)
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
  while ( is.na(l) || !is.numeric(l) || l%%1 != 0 || l>d || l>length(S) ) {
    cat("l value is not valid. l should be a positive integer lower or equal to d or the number of elements in S.")
    l <- readline(prompt = "Please enter a valid l : ")
    l <- suppressWarnings(as.numeric(l))
  }
  while ( is.na(alpha) || !is.numeric(alpha) || alpha <= 0 ) {
    cat("alpha value is not valid for CUT step. alpha should be a positive number.")
    alpha <- readline(prompt = "Please enter a valid alpha: ")
    alpha <- suppressWarnings(as.numeric(alpha))
  }
  while ( is.na(mu) || !is.numeric(mu) || mu <= 0 ) {
    cat("mu value is not valid for CUT step. mu should be a positive number.")
    mu <- readline(prompt = "Please enter a valid mu: ")
    mu <- suppressWarnings(as.numeric(mu))
  }
  while ( is.na(xi) || !is.numeric(xi) || xi <= 0 ) {
    cat("xi value is not valid for CUT step. xi should be a positive number.")
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
  base <- shapeSample(X=X,d=d)
  b_Sja <- base_Sja(S=S,A=A,base=base)
  #\.

  dTV_txy <- numeric(lenS)
  for (z in 1:lenS) {
    j <- dec_S[z]
    Sminusj <- dec_S[ -which( dec_S == j ) ]

    Q <- matrix(0,ncol=nrowA_pairs,nrow = nrow_subx)
    R <- matrix(0,ncol = lenA, nrow = nrow_subx)
    for (k in 1:nrow_subx) { #runs in all x_S
      Q[k,] <- dTV_sample(S=Sminusj,j=j,lenA=lenA,base=b_Sja,A_pairs=A_pairs,x_S=subx[k,])
      R[k,] <- sx(S=Sminusj,base=b_Sja,lenA=lenA,x_S=subx[k,],mu=mu,alpha=alpha,xi=xi)
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
#tn=|s(1111)+s(2111)|s(1111)+s(3111)|s(2111)+s(3111)|
#   |s(1211)+s(2211)|s(1211)+s(3211)|s(2211)+s(3211)|...
