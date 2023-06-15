#' The Foward Stepwise method.
#'
#' A function for inference in MTD Markov chains with FS method. Applies Forward Stepwise (FS) algorithm to estimate a relevant lag set \eqn{\Lambda} for MTD models.
#'
#' @param X A mixture transition distribution (MTD) chain sample.
#' @param d An upper bound for the chains order.
#' @param l Stop point for FS algorithm.
#' @param A The States space.
#' @param elbowTest If TRUE the function will choose the point where the \eqn{ max(\nu)} stops decreasing (that is, where the graph of the \eqn{ max(\nu)} changes its direction like an elbow).
#' If l is smaller than this point the function will simply return a set with l elements, so setting elbowTest=TRUE only makes a difference if the elbow point happens in the first l-1 elements of \eqn{ max(\nu)}.
#' @param warning If True may return warnings.
#' @param ... Used to accommodate any extra arguments passed by the [hdMTD()] function.
#'
#' @importFrom utils combn
#'
#' @details The "Forward Stepwise" (FS) algorithm is the first step of the "Forward Stepwise and Cut" (FSC) algorithm for inference in
#' Mixture Transition Distribution (MTD) models. Which consists in the application of the "Forward Stepwise" (FS) step followed by the CUT algorithm.
#' This method and its steps where developed by [Ost and Takahashi](https://arxiv.org/abs/2202.08007) and are specially useful for inference in High order MTD Markov chains.
#' This specific function will apply only the FS step of the algorithm and return an estimated relevant lag set of size l.
#'
#' @return Returns a estimated S set of relevant lags with size l .
#' @export
#' @examples
#' X <- perfectSample(MTDmodel(Lambda=c(2,4),A=c(0,1),w0=0.05),2000)
#'hdMTD_FS(X,A=c(0,1),d=5,l=2)
#'hdMTD_FS(X,d=5,l=2)
#'hdMTD_FS(X,d=4,l=3,elbowTest = TRUE)
#'hdMTD_FS(X,d=4,l=3,elbowTest = TRUE)
#'
hdMTD_FS <- function(X,d,l,A=NULL,elbowTest=FALSE,warning=FALSE,...){
  # Cheking inputs
  while ( is.na(l) || !is.numeric(l) || l%%1 != 0 || l>d ) {
    cat("l value is not valid for FS step. l should be a positive integer less than or equal to d.")
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
  #test if l is valid but creates too big a table
  xa=try(expand.grid(rep(list(A),l)),silent = TRUE)
  if(class(xa)=="try-error"){stop(paste0("The data set with all sequences of size l is too large, try a lower l."))}
  #\.
  #Gathering inputs
  A <- sort(A)
  lenA <- length(A)
  lenX <- length(X)
  base <- countsTab(X=X,d=d)
  A_pairs <- t(utils::combn(A,2))
  A_pairsPos <- t(utils::combn(1:lenA,2))
  nrowA_pairs <- nrow(A_pairs)
  #\.

  S <- NULL
  lenS <- 0
  maxnu <- numeric(l)
  while ( lenS < l ) {

    if( is.numeric(S) ){
      Sc=(1:d)[-S]
    }else{
      Sc <- 1:d
    }

    lenSc <- length(Sc)
    Sc <- sort(Sc,decreasing = TRUE)
    dec_S <- sort(S,decreasing = TRUE)
    nuj <- numeric(lenSc)

    for (z in 1:lenSc) {
      j <- Sc[z]
      b_Sja <- freqTab(S=S,j=j,A=A,countsTab=base)
      b_Sj <- freqTabSj(S=S,j=j,b_Sja,lenX=lenX,d=d)
      b_S <- freqTabSj(S=S,j=NULL,b_Sja,lenX=lenX,d=d)#if S=NULL: b_S<-matrix(c(0,lenX-d),ncol=2)
      ncolb_S <- ncol(b_S)

      if( lenS > 0) {
        PositNx_S <- which(b_S$Nx_Sj>0)
        #PositNx_S = {index of all sequences x_S : N(x_S)>0}
        subx <- b_S[,-ncolb_S]
      }else{
        PositNx_S <- 1
        subx <- matrix(0,ncol = 1) # needs to be in this format
      }

      lenPositNx_S <- length(PositNx_S)
      for (k in 1:lenPositNx_S) { #runs in all sequences x_S: N(x_S)>0
        t <- PositNx_S[k]
        cont <- 0

        PIs <- PI(S=dec_S,freqTabSj=b_Sj,x_S=subx[t,],lenX=lenX,
                  d=d)
        #(pi(xa_Sj),pi(xb_Sj),pi(xc_Sj),...)
        dTVs <- dTV_sample(S=dec_S,j=j,lenA=lenA,base=b_Sja,
                           A_pairs=A_pairs,x_S=subx[t,])
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
    if(elbowTest==TRUE & lenS>1){ ########################
      if(maxnu[lenS]>maxnu[lenS-1]){
        S <- S[-lenS]
        maxnu <- maxnu[-length(maxnu)]
        break
      }
    }
  }
  S
}

