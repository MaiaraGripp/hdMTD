#' A tibble with sample sequences counts
#'
#' Creates a tibble containing all size \code{d+1} sequences present in the
#' sample along with their absolute frequencies.
#'
#' @details The function will generate a tibble with \code{d+2} columns. For
#' each row, the first \code{d+1} columns will display a sequence
#' (of \code{d+1} elements) that appeared in the sample. The last column,
#' called \code{Nxa}, will contain the number of times each of these
#' sequences appeared in the sample. The number of rows in the outputted
#' tibble depends on \code{X}, since each row represents an unique size
#' \eqn{d+1} sequence that actually appeared in the chain. Let \eqn{|A|} be
#' the number of elements in the state space, and \eqn{d+1} the size of the
#' sequences, the number of rows varies between \eqn{1} and \eqn{|A|^{d+1}}.
#' Note that the function won't accept \code{X} or \code{d} with
#' \code{length(X)<=d+1}.
#'
#' @param X A vector or single-column data frame with a sample from a chain.
#' The first element must be the most recent.
#' @param d Numeric, must be a positive integer. Specifies the number of
#' elements in the sequences, which is \eqn{d+1}. Here, \eqn{d} is often taken as the
#' chain order or serves as an upper limit for it.
#'
#' @return A tibble with every size \code{d+1} sequence present in the sample and its
#'  absolute frequency.
#' @export
#' @importFrom dplyr %>%
#'
#'
#' @examples
#' countsTab(c(1,2,2,1,2,1,1,2,1,2),3)
#'
#' #Using test data.
#' countsTab(testChains[,1],2)
#'
countsTab <-function(X,d){
  ## Checks parameters
  X <- unlist(X)
  X <- checkSample(X)
  if(!is.numeric(d) || d%%1!=0 || d<=0 )stop("d must be a positive integer number")
  if (length(X)<=d+1) { stop("The sample size must be greater than d+1.") }

  n <- length(X)
  d1 <- d+1
  XTab <- NULL
  X <- rev(X)

  if( n-d >=d1 ){
      for (i in 1:d1) { # makes d1 matrices sample sequences
        aux <- (n-(i-1))%%d1
        XTab <- rbind( XTab, matrix( X[i:(n-aux)] ,ncol = d1,byrow = TRUE) )
      }
  }else{
    for (i in 1:(n-d)) { # if n+1<2*d1
      XTab <- rbind( XTab, matrix( X[i:(i+d)] ,ncol = d1,byrow = TRUE) )
    }
  }
  colnames(XTab) <- c(paste("x",seq(d,1),sep="" ),"a")
  count <- rep(1,n-d)
  XTab <- dplyr::as_tibble(cbind(XTab,count))
  XTab <- XTab %>%
    dplyr::group_by(XTab[,-(d1+1)]) %>%
    dplyr::summarise(Nxa=sum(count),.groups="drop")
  XTab
}


