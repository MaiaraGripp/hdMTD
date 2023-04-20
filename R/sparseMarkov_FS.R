#' A function for inference in MTD Markov chains with FS method.
#'
#' Applies Forward Stepwise (FS) algorithm to estimate a relevant lag set for MTD models.
#'
#' @param X A mixture transition distribution (MTD) chain sample.
#' @param A The States space.
#' @param d An upper threshold for the chains order.
#' @param l Stop point for FS algorithm.
#' @param warning If True may return warnings
#'
#' @importFrom utils combn
#'
#' @details The "Forward Stepwise and Cut" (FSC)is an algorithm for inference in
#' Mixture Transition Ditribution (MTD) models.
#' It consists in the application of the "Forward Stepwise" (FS) step followed by the CUT algorithm.
#' This method was developed by [Ost and Takahashi](https://arxiv.org/abs/2202.08007) and is specially useful for High order MTD Markov chains.
#' This function will apply only the FS step of the algorithm.
#'
#' @return Returns a estimated S set of relevant lags with size l .
#' @export
#'
sparseMarkov_FS <- function(X,A=NULL,d,l,warning=FALSE){
  # Cheking inputs
  while ( is.na(l) || !is.numeric(l) || l%%1 != 0 || l>d ) {
    cat("l value is not valid for FS step. l should be a positive integer lower or equal to d.")
    l <- readline(prompt = "Please enter a valid l : ")
    l <- suppressWarnings(as.numeric(l))
  }
  X <- checkSample(X)
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
  if( !is.numeric(d) || d<2 || (d %% 1)!=0 ){
    stop("The order d must be an integer number greater than 2.")
  }
  #\.
  #Gathering inputs
  A <- sort(A)
  lenA <- length(A)
  lenX <- length(X)
  base <- shapeSample(X=X,d=d)
  A_pairs <- t(utils::combn(A,2))
  A_pairsPos <- t(utils::combn(1:lenA,2))
  nrowA_pairs <- nrow(A_pairs)
  #\.

  #maxnuj <- numeric(l)
  S <- NULL
  lenS <- 0
  while ( lenS < l ) {

    if( is.numeric(S) ){
      subx <- as.matrix(expand.grid(rep(list(A),lenS))[,lenS:1],
                        ncol=lenS)
      nrow_subx=nrow(subx)
      Sc=(1:d)[-S]
    }else{
      subx=matrix(0,ncol = 1) # needs to be in this format
      nrow_subx=1
      Sc=1:d
    }
    lenSc <- length(Sc)
    Sc <- sort(Sc,decreasing = TRUE)
    dec_S <- sort(S,decreasing = TRUE) # functions dTV_sample and PI need a decreasing S

    nuj <- numeric(lenSc)
    for (z in 1:lenSc) {
      j <- Sc[z]

      b_Sja <- base_Sja(S=S,j=j,A=A,base=base)
      b_Sj <- base_Sj(S=S,j=j,b_Sja,lenX=lenX,d=d)
      b_S <- base_Sj(S=S,j=NULL,b_Sja,lenX=lenX,d=d)#if S=NULL b_S<-matrix(c(0,lenX-d),ncol=2)
      ncolb_S <- ncol(b_S)

      if( lenS > 0) {
        PositNx_S <- which(b_S$Nx_Sj>0)
      }else{
        PositNx_S <- 1 }
#PositNx_S = {index of all sequences x_S : N(x_S)>0}

      lenPositNx_S <- length(PositNx_S)
      for (k in 1:lenPositNx_S) { #runs in all sequences x_S: N(x_S)>0
        t <- PositNx_S[k]
        cont <- 0

        PIs <- PI(S=dec_S,base=b_Sj,x_S=subx[t,],lenX=lenX,
                  d=d)
#(pi(xa_Sj),pi(xb_Sj),pi(xc_Sj),...)
        dTVs <- dTV_sample(S=dec_S,j=j,lenA=lenA,base=b_Sja,
                           A_pairs=A_pairs,x_S=subx[t,]) #receives b_Sja
#(dTv_xS[p(.|a_j),p(.|b_j)],dTv_xS[p(.|a_j),p(.|c_j)]...)

        for (y in 1:nrowA_pairs) {# runs in pairs (b,c) such that b\in A, c\in A and b\neq c
          cont <- cont + prod(PIs[A_pairsPos[y,]])*dTVs[y]
        }

        PI_xS <- as.numeric(b_S[t,ncolb_S]/(lenX-d))#will be 1 when S=NULL.
        nuj[z] <- nuj[z]+cont/PI_xS
      }
    }
    s <- Sc[which(nuj==max(nuj))]
    S <- c(S,s)
    lenS <- length(S)
    #maxnuj[lenS] <- max(nuj)
  }
  #cbind(S,maxnuj)###
  S
}

