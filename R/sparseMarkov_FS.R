#' Aplies Forward Stepwise (FS) algorithm
#'
#' @param X A Markov Chain sample with mixture transition distribution.
#' @param l Stop point for FS algorithm.
#' @param A States space.
#' @param d Chains order.
#'
#' @return Returns a set S of size l that contains relevant lags.
#' @export
#'
sparseMarkov_FS <- function(X,l,A,d){
  # Cheking l


  while ( is.na(l) || !is.numeric(l) || l%%1 != 0 || l>d ) {
    cat("l value is not valid for FS step. l should be a positive integer number lower or equal to d.")
    l <- readline(prompt = "Please enter a valid l : ")
    l <- suppressWarnings(as.numeric(l))
  }
  #\.
  # Cheking A
  if ( !all( sort(unique(X)) %in% sort(A) ) ) {
    stop("Check the states space, it must have all possible configurations of X.")
  }
  #\.
  #Gathering inputs from inputs
  lenX <- length(X)
  base <- shapeSample(X=X,d=d)
  A <- sort(A)
  lenA <- length(A)
  A_pairs <- t(utils::combn(A,2))
  nrowA_pairs <- nrow(A_pairs)
  #\.

  S <- NULL
  lenS <- 0
  while ( lenS < l ) {

    if( is.numeric(S) ){
      subx <- as.matrix(expand.grid(rep(list(A),lenS))[,lenS:1],ncol=lenS)
      nrow_subx=nrow(subx)
      Sc=(1:d)[-S]
    }else{
      subx=matrix(0,ncol = 1) # needs to be in this format
      nrow_subx=1
      Sc=1:d
    }
    lenSc <- length(Sc)
    Sc <- sort(Sc,decreasing = TRUE)
    dec_S <- sort(S,decreasing = TRUE) # functions dTV and PI need a decreasing S

    nuj <- numeric(lenSc)
    for (z in 1:lenSc) {
      j <- Sc[z]

      b_Sja <- base_Sja(S=S,j=j,A=A,base=base)
      b_Sj <- base_Sj(S=S,j=j,b_Sja,lenX=lenX,d=d)
      b_S <- base_Sj(S=S,j=NULL,b_Sja,lenX=lenX,d=d)#if S=NULL b_S<-matrix(c(0,lenX-d),ncol=2)
      ncolb_S <- ncol(b_S)

      if( lenS > 0) { PositNx_S <- which(b_S$Nx_Sj>0) }else{ PositNx_S <- 1 }
#PositNx_S is de index of all sequences x_S : N(x_S)>0

      lenPositNx_S <- length(PositNx_S)
      for (k in 1:lenPositNx_S) { #runs in all sequences x_S: N(x_S)>0
        t <- PositNx_S[k]
        cont <- 0

        PIs <- PI(S=dec_S,base=b_Sj,x_S=subx[t,],lenX=lenX,d=d) #receives b_Sj
#(pi(xa_Sj),pi(xb_Sj),pi(xc_Sj),...)
        dTVs <- dTV(S=dec_S,j=j,lenA=lenA,base=b_Sja,A_pairs=A_pairs,x_S=subx[t,]) #receives b_Sja
#(dTv_xS[p(.|a_j),p(.|b_j)],dTv_xS[p(.|a_j),p(.|c_j)]...)

        for (y in 1:nrowA_pairs) {#runs in pairs (b,c) such that b\in A, c\in A and b\neq c
          cont <- cont + prod(PIs[A_pairs[y,]])*dTVs[y]
        }

        PI_xS <- as.numeric(b_S[t,ncolb_S]/(lenX-d))#will be 1 when S=NULL.
        nuj[z] <- nuj[z]+cont/PI_xS
      }
    }
    s <- Sc[which(nuj==max(nuj))]
    S <- c(S,s)
    lenS <- length(S)
  }
  S
}

