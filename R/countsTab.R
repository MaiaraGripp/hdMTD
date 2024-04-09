#' Table of sample sequence counts.
#'
#' Creates a table with sample sequences and their absolute frequencies.
#'
#' @param X A sample from an MTD Markov Chain.
#' @param d An upper bound for the chain's order.
#'
#' @return A table with every size d+1 sequence in the sample and its absolute frequency.
#' @export
#' @importFrom dplyr %>%
#'
#' @details The function will make a tibble with \eqn{d+2} columns. For each row, the first \eqn{d+1} columns will
#'   have a sequence of size \eqn{d+1} that appeared in the sample. The last column, called \code{Nxa},
#'   will contain the number of times each of these sequences appeared in the sample.
#'
#' @examples
#' countsTab(c(1,2,2,1,2,1,1,2,1,2),3)
#' countsTab(c(0,2,0,2,0,2,1,1,0,0,1,2,1,2,1),4)
#'
countsTab <-function(X,d){
  X <- unlist(X)
  # Checking restrictions
  X <- checkSample(X)
  if (length(X)<d+1) { stop("The sample size must be greater than d+1.") }
  #\.

  n <- length(X)
  d1 <- d+1
  XTab <- NULL
  for (i in 1:d1) {
    aux <- (n-(i-1))%%d1
    XTab <- rbind( XTab, matrix( X[i:(n-aux)] ,ncol = d1,byrow = TRUE) )
  }
  colnames(XTab) <- c(paste("x",seq(d,1),sep="" ),"a")

  count <- NULL
  XTab <- dplyr::as_tibble(cbind(XTab,count=1))
  XTab <- XTab %>%
            dplyr::mutate(count = as.numeric(count))
  XTab <- XTab %>%
            dplyr::group_by(XTab[,1:d1]) %>%
            dplyr::summarise(Nxa=sum(count),.groups="drop")
  XTab
}


