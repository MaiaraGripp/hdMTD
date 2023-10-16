#' MTD samples for tests.
#'
#' A data frame with samples perfectly sampled from a Mixture Transition
#' Distribution model and their parameters. Each chain was perfectly sampled from the same MTD model.
#'
#' The samples were created as follows:
#'
#' set.seed(1)
#'
#' p_j = list("p_-2"=matrix(c(0.73449134, 0.2655087,0.04906214, 0.9509379),
#'  byrow=T,ncol=2),"p_-3"=matrix(c(0.6278761,0.3721239,0.3415812,0.6584188),
#'   byrow=T,ncol=2),"p_-9"=matrix(c(0.4271466,0.5728534,0.6650324, 0.3349676),
#'    byrow=T,ncol=2))
#'
#' MTDseed1 <- MTDmodel(Lambda=c(2,3,9),A=c(0,1),w0=0.05,
#'  w_j = c(0.25,0.3,0.4),p_j=p_j).
#'
#'testChainCol1 <- perfectSample(MTDseed1,3000)
#'
#'testChainCol2 <- perfectSample(MTDseed1,3000)
#'
#'testChainCol3 <- perfectSample(MTDseed1,3000)
#'
#' @format A data frame with 3000 rows and 3 perfectly sampled chains from the same MTD model.
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
#' data(testChain1)
"testChain1"
