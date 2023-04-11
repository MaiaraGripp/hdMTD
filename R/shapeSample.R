#' Creates a table with sample information
#'
#' @param X A Markov chain sample.
#' @param d The order of the Markov chain.
#'
#' @return A table with every sequence in the sample and its frequency.
#' @export
#' @importFrom dplyr %>%
#'
#' @details The function will make a tibble with \eqn{d+2} columns. For each row, the first \eqn{d+1} columns will
#'   have a sequency of size \eqn{d+1} that appeared in the sample. The last column, called \code{Nxa},
#'   will contain the number of times each of these sequences appeared in the sample.
#'
#'
#' @examples
#' shapeSample(c(1,2,2,1,2,1,1,2,1,2),3)
#' shapeSample(c(0,2,0,2,0,2,1,1,0,0,1,2,1,2,1),4)
shapeSample <-function(X,d){
  # Checking restrictions
  X <- unlist(X)
  checkSample(X)
  if (length(X)<d+1) { stop("The sample size must be greater than d+1.") }
  #\.

  #

  n <- length(X)
  d1 <- d+1
  XTab <- NULL
  for (i in 1:d1) {
    aux <- (n-(i-1))%%d1
    XTab <- rbind( XTab, matrix( X[i:(n-aux)] ,ncol = d1,byrow = TRUE) )
  }
  colnames(XTab) <- c(paste("x",seq(d,1),sep="" ),"a")

  # Adding counts do XTab:
  ##  Column Nxa: how many times each sequence appeared in X.
  count <- NULL #solve problem with global binding in check()
  XTab <- dplyr::as_tibble(cbind(XTab,count=1))
  XTab <- XTab %>%
            dplyr::mutate(count = as.numeric(count))
  XTab <- XTab %>%
            dplyr::group_by(XTab[,1:d1]) %>%
            dplyr::summarise(Nxa=sum(count),.groups="drop")

  XTab
}
