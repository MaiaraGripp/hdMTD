#' MTD samples for tests.
#'
#' A data frame with samples perfectly sampled from a Mixture Transition
#' Distribution model and their parameters. Each chain was perfectly sampled from the same MTD model.
#'
#' Each chain was perfectly sampled from the same MTD model created as follows:
#'
#' set.seed(1)
#'
#' pj <- list("p-1"=matrix(c(0.1,0.1,0.8,0.4,0.4,0.2,0.5,0.3,0.2),
#' byrow = T,ncol = 3),"p-30"=matrix(c(0.05,0.2,0.75,0.4,0.4,0.2,0.3,0.3,0.4),
#' byrow = T,ncol = 3))
#'
#' MTDseed1 <- MTDmodel(Lambda=c(1,30),A=c(1,2,3),lam0=0.05,
#' lamj = c(0.35,0.6),pj=pj)
#'
#'testChainCol1 <- perfectSample(MTDseed1,3000)
#'
#'testChainCol2 <- perfectSample(MTDseed1,3000)
#'
#'testChainCol3 <- perfectSample(MTDseed1,3000)
#'
#' @format A data frame with 3000 rows and 3 perfectly sampled chains.
#'
#' \describe{
#'   \item{testChain1}{perfectSample1,N=3000}
#'   \item{testChain2}{perfectSample2,N=3000}
#'   \item{testChain3}{perfectSample3,N=3000}
#' }
#'
#' @source {Created in-house to serve as example.}
#'
#' @examples
#' data(testChain2)
"testChain2"
