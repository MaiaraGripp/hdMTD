#' MTD samples for tests
#'
#' A tibble with chains perfectly sampled from a MTD model. Each chain
#' was sampled from the same MTD model. Hence the differences between the
#' samples are due to randomness within the perfect sample algorithm.
#'
#' The MTD model from which the chains were sampled was created as follows:
#'
#' set.seed(1)
#'
#' pj <- list("p-1"=matrix(c(0.1,0.1,0.8,0.4,0.4,0.2,0.5,0.3,0.2), byrow = T,ncol = 3),
#'
#' "p-30"=matrix(c(0.05,0.2,0.75,0.4,0.4,0.2,0.3,0.3,0.4), byrow = T,ncol = 3))
#'
#' MTDseed1 <- MTDmodel(Lambda=c(1,30),A=c(1,2,3),lam0=0.05,
#' lamj = c(0.35,0.6),pj=pj)
#'
#'testChain1 <- perfectSample(MTDseed1,3000)
#'
#'testChain2 <- perfectSample(MTDseed1,3000)
#'
#'testChain3 <- perfectSample(MTDseed1,3000)
#'
#'testChains <- dplyr::as_tibble(cbind(testChain1,testChain2,testChain3))
#'
#' @format A tibble with 3000 rows and 3 perfectly sampled chains.
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
#' data(testChains)
"testChains"
