#' Oscillations of an MTD Markov chain
#'
#' Calculates the oscillations of an MTD model object or estimates the oscillations of a chain sample.
#'
#' @param x Must be an MTD object or a chain sample.
#' @param ... Additional parameters that might be required. Such as:
#'
#' \code{S}: If \code{x} is a chain sample the user should provide a set of lags
#' for which he wishes to estimate the oscillations. It must be labeled as \code{S}, an in
#' this scenario the function takes an upper bound for the order as \code{d=max{S}}.
#'
#' \code{A}: If \code{x} is a chain sample, and there may be elements in \code{A} that did not
#' appear in \code{x}, the state space should be specified, and it must be labeled as \code{A}.
#'
#' @details The oscillation for a certain lag \eqn{j} of an MTD model
#' ( {\eqn{\delta_j}: \eqn{j \in \Lambda} }), is the product of the weight \eqn{\lambda_j}
#' multiplied by the maximum of the total variation distance between the distributions in a
#' stochastic matrix \eqn{p_j}.
#' \deqn{\delta_j = \lambda_j\max_{b,c \in \mathcal{A}} d_{TV}(p_j(\cdot | b), p_j(\cdot | c)).}
#' So, if \code{x} is an MTD object, the parameters \eqn{\Lambda}, \eqn{\mathcal{A}}, \eqn{\lambda_j},
#' and \eqn{p_j} are inputted through, respectively, the entries \code{Lambda}, \code{A},
#' \code{lambdas} and the list \code{pj} of stochastic matrices. Hence, an oscillation \eqn{\delta_j}
#' may be calculated for all \eqn{j \in \Lambda}.
#'
#' @details If we wish to estimate the oscillations from a sample, then \code{x} must be a chain,
#' and \code{S}, a vector representing a set of lags, must be informed. This way the transition
#' probabilities can be estimated. Let \eqn{\hat{p}(\cdot| x_S)} symbolize an estimated distribution
#' in \eqn{\mathcal{A}} given a certain past \eqn{x_S} ( which is a sequence of elements of
#'  \eqn{\mathcal{A}} where each element occurred at a lag in \code{S}), and
#'  \eqn{\hat{p}(\cdot|b_jx_S)} an estimated distribution given past \eqn{x_S} and that the symbol
#'  \eqn{b\in\mathcal{A}} occurred at lag \eqn{j}.
#'  If \eqn{N} is the sample size, \eqn{d=}\code{max{S}} and \eqn{N(x_S)} is the
#'  number of times the sequence \eqn{x} appeared in the at the lags in \code{S}, then
#' \deqn{\delta_j = \max_{c_j,b_j \in \mathcal{A}} \frac{1}{N-d}\sum_{x_{S} \in \mathcal{A}^{S}} N(x_S)d_{TV}(\hat{p}(. | b_jx_S), \hat{p}(. | c_jx_S) )}
#' is the estimated oscillation for a lag \eqn{j \in \{1,\dots,d\}\setminus}\code{S}. Note that \eqn{\mathcal{A}^S} is the space of
#' sequences of \eqn{\mathcal{A}} indexed by \code{S}.

#' @return If the \code{x} parameter is an MTD object, it will provide the oscillations for
#' each element in \code{Lambda}. In case \code{x} is a chain sample, it estimates the oscillations
#'  for a user-inputted set \code{S} of lags.
#'
#' @export oscillation
#' @examples
#' oscillation( MTDmodel(Lambda=c(1,4),A=c(2,3) ) )
#' oscillation(MTDmodel(Lambda=c(1,4),A=c(2,3),lam0=0.01,lamj=c(0.49,0.5),
#' pj=list(matrix(c(0.1,0.9,0.9,0.1),ncol=2)), single_matrix=TRUE))
#'
oscillation <- function(x,...){ UseMethod("oscillation") }

#' @export
oscillation.MTD <- function(x,...){
## check Inputs
  checkMTD(x)
  lenA <- length(x$A) #number of rows/cols in each p_j
  rows <- t(utils::combn(lenA,2))
  y <- x$lambdas[-1]*sapply(x$pj,dTV_pj,rows)
  names(y) <- paste0("-",x$Lambda)
  y
}

#' @export
oscillation.default <- function(x,...){
## Check inputs
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
      any(A%%1 !=0)   )stop("States space A must be a vector of integers with at least two numbers.")
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


  for(i in 1:lenS){ #for each element in S
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
    for (t in 1:lenPositiveNx_S) { #calculates total var. dist. of distributions conditioned
# in each past that appears in the sample (but variating a single symbol of said past).
        dtv_xS[t,] <- dTV_sample(S=Z,j=j,lenA=lenA,base=b_Sja,
                                    A_pairs=A_pairs,x_S=subx[t,])
    }
    if(lenS>1){
      NxS_dtvxS <- sweep(dtv_xS,MARGIN=1,t(b_S[PositiveNx_S,lenS]),'*')#dtv|p(.|xz aj)-p(.|xz bj)|*Nxz
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


