#' Oscillations of a MTD Markov chain
#'
#' Calculates the oscillations of a MTD model object or estimates the oscillations of a chain sample.
#'
#' @param x Must be a MTD object or a chain sample.
#' @param ... bla
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
    stop("S must be informed. S should be a number or a vector of positive integer numbers representing past times for which you want to calculate the oscillations.")
  }
  A <- params$A
  if(length(A)==0){
    A <- unique(x)
  }

  lenX <- length(x)
  S <- sort(S,decreasing = TRUE)
  lenS <- length(S)
  lenA <- length(A)
  A_pairs <- t(utils::combn(A,2))
  base <- countsTab(X=x,d=S[1])
  y <- numeric(lenS)
  Z <- NULL

  for(i in 1:lenS){
    j <- S[i]
    if(lenS>1){
      Z <- S[-i]
    }
    b_Sja <- freqTab(S=Z,j=j,A=A,countsTab=base)
    if(lenS>1){
      b_S <- freqTabSj(S=Z,j=NULL,b_Sja,lenX=lenX,d=S[1])
      PositNx_S <- which(b_S$Nx_Sj>0)
      subx <- b_S[,-lenS]
    }else{
      PositNx_S <- 1
      subx <- matrix(0,ncol = 1)
    }

    lenPositNx_S <- length(PositNx_S)
    dtv_xS <- numeric(lenPositNx_S)
    for (t in 1:lenPositNx_S) {
      dtv_xS[t] <- max(dTV_sample(S=Z,j=j,lenA=lenA,base=b_Sja,
                                  A_pairs=A_pairs,x_S=subx[t,]))
    }
    y[i] <- max(dtv_xS)
  }
  names(y) <- paste0("-",S)
  y <- rev(y)
  y
  #print("The implemented method can currently only calculate oscillations for MTD objects.")
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

