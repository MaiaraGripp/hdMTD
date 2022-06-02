#' Applies Cut algorithm
#'
#' @param X A mixture transition distribution chain sample.
#' @param A The states space.
#' @param d The chains order.
#' @param alpha A parameter of CUT.
#' @param mu A parameter of CUT.
#' @param xi A parameter of CUT.
#' @param S A set of relevant lags.
#'
#' @return Returns a set that contains relevant lags.
#' @export
#'
sparseMarkov_CUT <- function(X, A, d, alpha, mu, xi, S=1:d){

  #Checking inputs
  if ( is.numeric(S) ) {
    if ( length(S) == 1 ) { stop("S must have at least 2 entries")}
    if ( any( S%%1 != 0 ) ) { stop("All elements in S must be integers")}
  }else{
    stop("S must be a numeric vector with at least 2 entries")
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

  lenA <- length(A)
  lenS <- length(S)
  subx <- as.matrix(expand.grid(rep(list(A),lenS-1))[,(lenS-1):1],ncol=lenS-1)
  nrow_subx <- nrow(subx)
  dec_S <- sort(S,decreasing = TRUE)
  A_pairs <- t(utils::combn(A,2))
  nrowA_pairs <- nrow(A_pairs)
  base <- shapeSample(X=X,d=d)
  b_Sja <- base_Sja(S=S,A=A,base=base)

  dTV_txy <- numeric(lenS)
  for (z in 1:lenS) {
    j <- dec_S[z]
    Sminusj <- dec_S[ -which( dec_S == j ) ]

    Q <- matrix(0,ncol=nrowA_pairs,nrow = nrow_subx)
    R <- matrix(0,ncol = lenA, nrow = nrow_subx)
    for (k in 1:nrow_subx) { #runs in all x_S
      Q[k,] <- dTV(S=Sminusj,j=j,lenA=lenA,base=b_Sja,A_pairs=A_pairs,x_S=subx[k,])
      R[k,] <- sx(S=Sminusj,base=b_Sja,lenA=lenA,x_S=subx[k,],mu=mu,alpha=alpha,xi=xi)
    }
    colnames(Q) <- apply(A_pairs, 1, paste0, collapse="x")
    rownames(Q) <- apply(subx, 1, paste0, collapse="")
    assign(paste0("dtv_j",j),Q)
    colnames(R) <- A
    rownames(R) <- apply(subx, 1, paste0,collapse="")
    assign(paste0("sx_Sj",j),R)

    dTV_txy[z] <- max(Q-txy(R=R,A_pairs=A_pairs))
  }
  S <- sort(S,decreasing = TRUE)
  S <- S[which(dTV_txy>0)]
  S
}
