
#' Creates a table with sample informartion.
#'
#' @param X A Markov Chain sample.
#' @param d The order of the Markov chain.
#'
#' @return A table with every sequence in the sample and its frequency.
#' @export
#' @importFrom dplyr %>%
#'
#' @details The function will make a tibble with \eqn{d+2} columns. For each row, the first \eqn{d+1} columns will
#' have a sequency of size \eqn{d+1} that appeared in the sample. The last column, called \code{Nxa},
#' will contain the number of times each of these sequences appeared in the sample.
#'
#'The algorithm to obtain a tibble with all sequences that appear in the sample consists in generating \code{d1} smaller matrices and grouping then by rows.
#' For exemple, let \code{d=3} and the sample \code{X} be a sequence of size \code{n=13}, where
#'\eqn{X=(x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_9,x_{10},x_{11},x_{12},x_{13})}.
#'
#'For the first matrix beggin from de element \eqn{x_1}.
#'
#'\eqn{\begin{bmatrix}x_1 & x_2 & x_3 & x_4\\ x_5 & x_6 & x_7 & x_8\\ x_9 & x_{10} & x_{11} & x_{12}\end{bmatrix}}
#'
#'For the second matrix discard the element \eqn{x_1} and beggin from element \eqn{x_2}. Than attach it to the last matrix.
#'
#' \eqn{\begin{bmatrix} x_1 & x_2 & x_3 & x_4\\ x_5 & x_6 & x_7 & x_8\\ x_9 & x_{10} & x_{11} & x_{12}\end{bmatrix}\\\begin{bmatrix}x_2 & x_3 & x_4 & x_5\\ x_6 & x_7 & x_8 & x_9\\ x_{10} & x_{11} & x_{12} & x_{13}\end{bmatrix}}
#'
#'Repeat \code{d1} times, for the last matrix discard \eqn{x_{d}} and beggin from \eqn{x_{d1}}.
#'
#'\eqn{\begin{bmatrix}x_1 & x_2 & x_3 & x_4\\ x_5 & x_6 & x_7 & x_8\\ x_9 & x_{10} & x_{11} & x_{12}\\ x_2 & x_3 & x_4 & x_5\\ x_6 & x_7 & x_8 & x_9\\ x_{10} & x_{11} & x_{12} & x_{13}\\ \vdots & & & \\ x_4 & x_5 & x_6 & x_7\\ x_8 & x_9 & x_{10} & x_{11} \end{bmatrix}}
#'
#'Note that, for every loop, the last elements of the sample are discarted, if needed, to fit the matrices.
#'
#' @examples
#' shapeSample(c(1,2,2,1,2,1,1,2,1,2),3)
#' shapeSample(c(0,2,0,2,0,2,1,1,0,0,1,2,1,2,1),4)
shapeSample <-function(X,d){
  #checking restrictions
  checkSample(X)
  if (length(X)<d+1){stop("The sample size must be greater than d+1.")}
  if( !is.numeric(d) | d<2 | (d %% 1)!=0 ){
    stop("d must be an integer number greater than 2.")
  }
  #\

  ##generating TableX:
  ##  a table whose rows are the sequences in X with size d+1.
  n <- length(X)
  d1 <- d+1
  TableX <- NULL
  for (i in 1:d1) {
    aux <- (n-(i-1))%%d1 #aux: the number of elements to burn at the end of X
    TableX <- rbind( TableX, matrix( X[i:(n-aux)] ,ncol = d1,byrow = TRUE) )
  }
  colnames(TableX) <- c(paste("x",seq(d,1),sep="" ),"a")

  ##Adding counts do TableX:
  ##  Column Nxa: how many times each sequence appeared in X.
  count=NULL #this is to not generate a problem with global binding in check()
  TableX <- dplyr::as_tibble(cbind(TableX,count=1))
  TableX <- TableX %>% dplyr::mutate(count = as.numeric(count))
  TableX <- TableX %>% dplyr::group_by(TableX[,1:d1]) %>% dplyr::summarise(Nxa=sum(count),.groups="drop")

  return(TableX)
}
