#' The total variation distance between distributions
#'
#' Calculates the total variation distance between distributions conditioned
#' in a given past.
#'
#' @details About the \code{base} argument format: this function was created to calculate
#' the total variation distance between distributions estimated through the [freqTab()] function.
#' Hence, typically, \code{base} is an output from [freqTab()], and its format is imperative
#' for this function to work. If you wish to use [dTV_sample()] to estimate total variation
#' distances for distributions acquired otherwise, please note that they must be inputted through
#' the \code{base} argument, following the same pattern as an output from [freqTab()] (i.e.,
#' column names must match, the actual probabilities should be named 'qax_Sj', and the past
#' and present lags should be separated into different columns with appropriate column names.
#' However, the columns 'Nx_Sj' and 'Nxa_Sj' do not need to be present in \code{base}).
#'
#' When you provide the state space \code{A} to this function, it will use it to calculate
#' the arguments \code{lenA <- length(A)} and \code{A_pairs <- t(utils::combn(A, 2))}. While
#' \code{A} is typically used, this function also allows for the input of \code{lenA} and
#' \code{A_pairs} instead of \code{A}. This flexibility is particularly useful within the
#' \code{hdMTD} package functions where [dTV_sample()] may be repeatedly used in loops with
#'  the same \code{A}, optimizing computational efficiency.
#'
#' @param S A numeric vector of positive integers or \code{NULL}. Represents
#' a set of past lags. The distributions from which this function will calculate the total
#' variation distance are conditioned on a fixed sequence indexed by \code{S} ( the user must
#'  also input the sequence through the argument \code{x_S}).
#' @param j An integer number. Represents a lag \code{j} in the \eqn{complement} of \code{S}.
#'  The symbols indexed by \code{j} vary along the state space \code{A}, altering the distribution
#'  through this single lag, and the size of this change is what this function seeks to measure.
#' @param A A vector with positive integers. Must have at least two different entries. \code{A}
#' represents the state space. You may leave \code{A=NULL} (default) if you provide the function
#' with the arguments \code{lenA} and \code{A_pairs} (see description below).
#' @param base  A data frame that must have \code{length(S)+1} columns representing past symbols,
#' one column for the present symbols, and a column with the transition probabilities. Each
#' row depicts a past sequence, a present symbol, and the probability of seeing this symbol after
#' that sequence occurred at times \code{S} and \code{j}. Given a \code{x_S} past sequence,
#' \code{base} must have at least every probability conditioned on \code{x_S} with the element
#' in lag \code{j} and the present varying across \code{A}, so at least \code{lenA^2} rows.
#' Since this function was created to estimate total variation distances between distributions
#' calculated through the [freqTab()] function, \code{base} is typically an output from [freqTab()]
#' (the \code{base} format must match that of an output from [freqTab()] with the argument
#' \code{complete=TRUE}); please check the [freqTab()] documentation and the *Details* section for
#' further information).
#' @param lenA An integer greater than 1 representing the lengt of the state space A, i.e.,
#'  \code{lenA <- length(A)}. If \code{A} is provided, there is no need to input \code{lenA}.
#' @param A_pairs A matrix with all possible pairs of elements from \code{A}. If \code{A} is provided,
#' there is no need to input \code{A_pairs}.
#' @param x_S A vector of length \code{length(S)} or \code{NULL}. If \code{S=NULL}, \code{x_S} will
#'  be set to \code{NULL}. \code{x_S} represents a sequence of symbols from \code{A} indexed by
#'  \code{S}. This sequence remains constant across the conditional distributions to be compared,
#'  representing the fixed configuration of the past.
#'
#' @return Returns a vector of total variation distances, where each entry corresponds to the
#' distance between a pair of distributions conditioned on the same fixed past \code{x_S},
#'  differing only in the symbol indexed by \code{j}, which varies across all distinct pairs
#'  of elements in \code{A}.Therefore, the vector contains as many entries as there are
#'  distinct pairs with elements of \code{A}.
#'
#' @importFrom dplyr %>%
#' @export
#' @examples
#' #creating base argument through freqTab function.
#' pbase <- freqTab(S=c(1,4),j=2,A=c(1,2,3),countsTab = countsTab(testChains[,2],d=5))
#' dTV_sample(S=c(1,2),j=4,A=c(1,2,3),base=pbase,x_S=c(2,3))
#' pbase <- freqTab(S=NULL,j=1,A=c(1,2,3),countsTab = countsTab(testChains[,2],d=5))
#' dTV_sample(S=NULL,j=1,A=c(1,2,3),base=pbase)
dTV_sample <- function(S,j,A=NULL,base,lenA=NULL,A_pairs=NULL,x_S){
  if(length(S)!=0){
      if(!is.vector(S) || any(S%%1 != 0) || any(S<1) || length(S)<0 )stop("S must be a positive integer vector, a number or NULL.")
  }
  if(length(j)!=0){
    if(length(j)!=1 || j%%1!=0 || j %in% S ){stop("j must be a integer number in the complement of S.")}
  }
  if(length(A)==0){
    if(length(lenA)==0 || length(A_pairs)==0){stop("Either the state space A must be provided, or both lenA (the number of elements in A) and A_pairs (all possible pairs of elements from A) must be provided.")
    }else{
      if( length(lenA)!=1 || lenA%%1!=0 ){stop("lenA must be an integer number.")}
      if( !is.matrix(A_pairs) || ncol(A_pairs)!=2 || any(A_pairs%%1 != 0) ){stop("A_pairs must be a two column matrix of integer numbers (even if there is only a single row), in each line should be a pair of different elements of the state space.")}
        }
    }else{
    if( length(A)<=1 ||
        !is.vector(A) ||
        any( A%%1 !=0 ) ||
        all(A==A[1]) ){stop("A must be a vector of integers with at least two distinct values.")}
      if( length(lenA)!=0 || length(A_pairs)!=0 ){warning("Since the state space A was provided, this function will set lenA = length(A) and A_pairs = t(utils::combn(A, 2))}, even if you have provided at least one of them.")}
    lenA <- length(A)
    A_pairs <- t(utils::combn(A,2))
    }
  lenS <- length(S)
  if(lenS==0){x_S <- NULL}else{
    if(length(x_S)!=lenS || !all(x_S %in% A_pairs)){stop("The number of elements in x_S vector must match the number of elements in S, and x_S must be a sequence with elements of A.")}
  }

  nrowA_pairs <- nrow(A_pairs)
  if ( is.numeric(S) ) {
      filtr_S <- paste0("x",S)
      B <- base
      B$test <- apply(B %>% dplyr::select_at(filtr_S),1,is_xS,x_S)
      C <- dplyr::filter(B,test==TRUE)
    }else{
      C <- base
    }
  if(nrow(C)!=lenA^2)stop("Not all distributions conditioned on all possible past symbols at lag j and past sequence x_S are fully present in base.")


  filtr_j <- paste0("x",j)
  disTV <- matrix(0,ncol=nrowA_pairs)
  for (i in 1:nrowA_pairs) {
    D <- C %>% dplyr::filter_at(filtr_j, ~. %in% A_pairs[i,])
    disTV[i] <- sum(abs(D$qax_Sj[1:lenA]-D$qax_Sj[(lenA+1):(2*lenA)]))/2
  }
  colnames(disTV) <- apply(A_pairs, 1, paste0, collapse="x")
  rownames(disTV) <- paste0(x_S,collapse = "")
  disTV
}


is_xS <- function(x,y) {
  return( all( x == y ) )
}



