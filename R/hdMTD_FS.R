#' The Forward Stepwise (FS) method for inference in MTD models
#'
#'  A function that estimates the set of relevant lags of an MTD model using the FS method.
#'
#' @param X A vector or single-column data frame containing a chain sample.
#' @param d A positive integer representing an upper bound for the chain order.
#' @param l A positive integer that sets the number of elements in the output vector.
#' @param A A vector with positive integers representing the state space. If not informed,
#' this function will set \code{A=unique(X)}.
#' @param elbowTest Logical. If TRUE, the function will use a special criterion to determine the length
#' of the estimated set of relevant lags. See *Details* below for more information.
#' @param warning Logical. If \code{TRUE}, informs the user if the state space was set as \code{A=unique(X)}.
#' @param ... Other parameters. This is used to accommodate any additional argument passed
#' to [hdMTD_FS()] through the [hdMTD()] function.
#'
#' @importFrom utils combn
#'
#' @details The "Forward Stepwise" (FS) algorithm is the first step of the "Forward Stepwise and Cut"
#' (FSC) algorithm for inference in Mixture Transition Distribution (MTD) models.
#' This method was developed by [Ost and Takahashi](https://arxiv.org/abs/2202.08007)
#' This specific function will only apply the FS step of the algorithm and return an estimated
#' relevant lag set of length \code{l}.
#'
#' @details If the algorithm determines that there are multiple lags equally important and more
#' important than all others, it will sample one of them uniformly.
#'
#' @details If \code{elbowTest=TRUE}, the function will use a new criterion to determine the length of
#' the estimated relevant lag set. In the FS algorithm, a certain quantity \eqn{\nu} is calculated
#' for each lag in \code{1:d}, and the lag with the greatest \eqn{\nu} is deemed important. This lag
#' is included in the output, and using this knowledge, the function proceeds to seek the next
#' important lag (the one with the highest \eqn{\nu} among the remaining ones). The process stops
#' when the output vector reaches length \code{l} if \code{elbowTest=FALSE}. If \code{elbowTest=TRUE},
#' the function will store a vector with the value of these greatest \eqn{\nu} at each step. It will
#' then look at this vector and identify the position of the smallest \eqn{\nu} among them. The
#' output will only keep the lags that came before the one responsible for this \eqn{\nu} value.
#'
#' @references
#' Ost, G. and Takahashi, D. (2022), ‘Sparse markov models for high-dimensional inference’.
#' [ arXiv:2202.08007](https://arxiv.org/abs/2202.08007)
#'
#' @return Returns a vector with the estimated relevant lag set using FS algorithm.
#' @export
#' @examples
#' X <- testChains[,1]
#'hdMTD_FS(X,d=5,l=2)
#'hdMTD_FS(X,d=4,l=3,elbowTest = TRUE)
#'
hdMTD_FS <- function(X,d,l,A=NULL,elbowTest=FALSE,warning=FALSE,...){
  # Cheking inputs
  while ( is.na(l) || !is.numeric(l) || l%%1 != 0 || l>d ) {
    cat("The l value is not valid for FS step. l should be a positive integer less than or equal to d.")
    l <- readline(prompt = "Please enter a valid l : ")
    l <- suppressWarnings(as.numeric(l))
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
  if( !is.numeric(d) || d<2 || (d %% 1)!=0 ){
    stop("The order d must be an integer number greater than 2.")
  }
  #test if l is valid but creates too big a table
  xa=try(expand.grid(rep(list(A),l)),silent = TRUE)
  if(class(xa)=="try-error"){stop(paste0("The dataset with all sequences of length l is too large. Please try a lower value for l."))}


  A <- sort(A)
  lenA <- length(A)
  lenX <- length(X)
  base <- countsTab(X=X,d=d)
  A_pairs <- t(utils::combn(A,2))
  A_pairsPos <- t(utils::combn(1:lenA,2))
  nrowA_pairs <- nrow(A_pairs)


  S <- NULL
  lenS <- 0
  maxnu <- numeric(l) # used if elbow=TRUE
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

    for (z in 1:lenSc) { # runs in all elements in the complement o S i.e. (1:d)\S
      j <- Sc[z]
      b_Sja <- freqTab(S=S,j=j,A=A,countsTab=base)
      b_Sj <- freqTabSj(S=S,j=j,b_Sja,lenX=lenX,d=d)
      b_S <- freqTabSj(S=S,j=NULL,b_Sja,lenX=lenX,d=d) #if S=NULL: b_S<-matrix(c(0,lenX-d),ncol=2)
      ncolb_S <- ncol(b_S)

      if( lenS > 0) {
        PositNx_S <- which(b_S$Nx_Sj>0) #{index of all sequences x_S : N(x_S)>0}
        subx <- b_S[,-ncolb_S]
      }else{
        PositNx_S <- 1
        subx <- matrix(0,ncol = 1) # needs to be in this format
      }

      lenPositNx_S <- length(PositNx_S)
      for (k in 1:lenPositNx_S) { #runs in all sequences x_S such that N(x_S)>0
        t <- PositNx_S[k]
        cont <- 0

        PIs <- PI(S=dec_S,freqTabSj=b_Sj,x_S=subx[t,],lenX=lenX,
                  d=d)

        dTVs <- dTV_sample(S=dec_S,j=j,lenA=lenA,base=b_Sja,
                           A_pairs=A_pairs,x_S=subx[t,])

        for (y in 1:nrowA_pairs) {# runs in pairs (b,c) such that b\in A, c\in A and b!=c
          cont <- cont + prod(PIs[A_pairsPos[y,]])*dTVs[y]
        }

        PI_xS <- as.numeric(b_S[t,ncolb_S]/(lenX-d))#will be 1 when S=NULL.
        nuj[z] <- nuj[z]+cont/PI_xS
      }
    }
    maxnu[lenS+1] <- max(nuj) # used if elbow=TRUE
    posMaxnu <- which( nuj==max(nuj) ) # the position of nu max (can be more than one)
    if( length(posMaxnu) > 1 ){ posMaxnu <- sample(posMaxnu,1) }# if FS chooses more than 1 lag, samples one uniformly.
    s <- Sc[posMaxnu]
    S <- c(S,s)
    lenS <- length(S)
  }
  if(elbowTest==TRUE){
    stp <- which(maxnu==min(maxnu))
    if(length(stp)>1){stp <- stp[length(stp)]}
    if(stp>1){stp <- stp-1}
    S <- S[1:stp]
  }
  S
}

