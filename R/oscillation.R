#' Oscillations of a MTD Markov chain
#'
#' Calculates the oscillations of a MTD model object or estimates the oscillations of a chain sample.
#'
#' @param x Must be a MTD object or a chain sample.
#' @param ... Additional parameters that might be required. Such as:
#'
#' S: If x is a chain sample,
#' the user should provide a set of lags, that must be labeled as S, for which he wishes to estimate the oscillations.
#'
#' A: If x is a chain sample, and there may be elements is A that did not appear in x, the state space A should be
#' specified, and it must be labeled as A.
#'
#'
#' @details The oscillations of a MTD model
#' (\eqn{\delta_k} for \eqn{k in \Lambda}), are the product of the weight \eqn{\lambda_k}
#'   multiplied by the maximum of the total variation distances between the
#'  distributions in the matrix p_k. These values are important because they
#'   measure the influence of a relevant lag k on the model.
#' @return If the x parameter is an MTD object, it will provide the oscillations for
#' each relevant element. In case x is a MTD chain sample, it estimates the oscillations
#'  for a user-inputted set S of lags.
#' @export oscillation
#' @examples oscillation( MTDmodel(Lambda=c(1,4),A=c(2,3) ) )
#' oscillation( MTDmodel(Lambda=c(1,4),A=c(2,3),w0=0.01 ) )
#' oscillation( MTDmodel(Lambda=c(1,4),A=c(2,3),w0=0.01,w_j=c(0.49,0.5)) )
#' oscillation( MTDmodel(Lambda=c(1,4),A=c(2,3),w0=0.01,w_j=c(0.49,0.5),
#' p_j=list(matrix(c(0.1,0.9,0.9,0.1),ncol=2)), single_matrix=TRUE))
#'
oscillation <- function(x,...){ UseMethod("oscillation") }

#' @export
oscillation.MTD <- function(x,...){
  checkMTD(x)
  lenA <- length(x$A) #number of rows/cols in each p_k
  rows <- t(combn(lenA,2))
  y <- x$lambdas[-1]*sapply(x$p_j,dTV_pj,rows)
  names(y) <- paste0("-",x$Lambda)
  class(y) <- "MTDoscillation"
  y
}

#' @export
oscillation.default <- function(x,...){
  x <- checkSample(x)
  params <- list(...)
  S <- params$S
  if(length(S) < 1  ||
     !is.numeric(S) ||
     any( S%%1 != 0) ){
    stop("A set of lags, for which the user intends to estimate the oscillations, must be informed and labeled as S.
         This set should be a number or a vector of positive integer numbers. Try S=c().")
  }
  A <- params$A
  if(length(A)==0){
    A <- unique(x)
  }
  if( length(A)<=1   ||
      any(A%%1 !=0)   )stop("States space A must be a numeric vector with at least two integers.")
  if ( !all( unique(x) %in% A ) ) {
    stop("Check the states space, it must include all states that occur in the sample.")
  }
  A <- sort(A)
  lenX <- length(x)
  S <- sort(S,decreasing = TRUE)
  lenS <- length(S)
  lenA <- length(A)
  A_pairs <- t(utils::combn(A,2))
  lenPairs <- nrow(A_pairs)
  base <- countsTab(X=x,d=S[1])
  y <- numeric(lenS)
  Z <- NULL


  for(i in 1:lenS){
    j <- S[i]
    if(lenS>1){
      Z <- S[-i] #Z are the lags in S\j
    }
    b_Sja <- freqTab(S=Z,j=j,A=A,countsTab=base)
    if(lenS>1){
      b_S <- freqTabSj(S=Z,j=NULL,b_Sja,lenX=lenX,d=S[1])
      PositiveNx_S <- which(b_S$Nx_Sj>0) #position in b_S of the x_Z that appear in sample
      subx <- b_S[PositiveNx_S,-lenS] # list of all x_Z
    }else{
      PositiveNx_S <- 1
      subx <- matrix(0,ncol = 1)
    }

    lenPositiveNx_S <- length(PositiveNx_S)
    dtv_xS <- matrix(0,ncol=lenPairs,nrow = lenPositiveNx_S)
    colnames(dtv_xS) <- apply(A_pairs,1,paste0,collapse="x")
    rownames(dtv_xS) <- apply(subx,1,paste0,collapse="")
    for (t in 1:lenPositiveNx_S) {
        dtv_xS[t,] <- dTV_sample(S=Z,j=j,lenA=lenA,base=b_Sja,
                                    A_pairs=A_pairs,x_S=subx[t,])
    }
    if(lenS>1){
      NxS_dtvxS <- sweep(dtv_xS,MARGIN=1,t(b_S[PositiveNx_S,lenS]),'*')
    }else{
      NxS_dtvxS <- dtv_xS*(lenX-S[1])#since S\j=NULL, N(x_S)=n-d, so N(x_S)/(n-d)=1 in the weighted mean
    }

    h_deltajbc <- apply(NxS_dtvxS, 2, sum)/(lenX-S[1])
    y[i] <- max(h_deltajbc)
  }
  names(y) <- paste0("-",S)
  y <- rev(y)
  y
}

#' @export
print.MTDoscillation <- function(x, ...){
 class(x) <- NULL
 cat("Calculating \U03B4\U2096 = \U03BB\U2096*max{b,c in A: dTV[p\U2096(.|b),p\U2096(.|c)]}, \n")
 cat("for each k in Lambda: \n")
 cat("\n")
 print(x)
 return(invisible(x))
}

