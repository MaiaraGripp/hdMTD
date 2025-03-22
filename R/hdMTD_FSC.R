#' Forward Stepwise and Cut method for inference in MTD models
#'
#' A function for inference in MTD Markov chains with FSC method. This function estimates the relevant
#' lag set \eqn{\Lambda} of an MTD model through the FSC algorithm.
#'
#' @param X A vector or single-column data frame containing a chain sample.
#' @param d A positive integer representing an upper bound for the chain order.
#' @param l A positive integer that sets the number of elements in the output vector.
#' @param A A vector with positive integers representing the state space. If not informed,
#' this function will set \code{A <- sort(unique(X))}.
#' @param alpha A positive real number used in the CUT threshold (which determines if two
#' distributions can be considered different). The larger the \code{alpha}, the greater
#' the distance required to consider that there is a difference between a set of distributions.
#' Defaulted to 0.05.
#' @param mu A positive real number such that \eqn{\code{mu}>(e^{\code{mu}}-1)/2}. \code{mu}
#' is also a component of the same threshold as \code{alpha}.
#' @param xi A positive real number, \code{xi} is also a component of the same threshold as
#'  \code{alpha}.
#' @param ... Additional arguments (not used in this function, but maintained for compatibility with [hdMTD()].
#'
#' @details The "Forward Stepwise and Cut" (FSC) is an algorithm for inference in
#' Mixture Transition Distribution (MTD) models. It consists
#' in the application of the "Forward Stepwise" (FS) step followed by the CUT algorithm.
#' This method and its steps where developed by [Ost and Takahashi](http://jmlr.org/papers/v24/22-0266.html)
#' and are specially useful for inference in high-order MTD Markov chains.
#'
#' @references
#' Ost, G. and Takahashi, D. Y. (2023), "Sparse markov models for high-dimensional inference",
#' Journal of Machine Learning Research 24(279), 1-54. [URL: http://jmlr.org/papers/v24/22-0266.html](http://jmlr.org/papers/v24/22-0266.html)
#'
#' @return Returns a vector with the estimated relevant lag set using FSC algorithm.
#' @export
#' @examples
#' X <- testChains[,1]
#' hdMTD_FSC(X,4,3,alpha=0.02)
#' hdMTD_FSC(X,4,2,alpha=0.001)
#'
hdMTD_FSC <- function(X, d, l, alpha = 0.05, mu = 1, xi = 0.5, A = NULL, ...){
    # Validate inputs
    X <- checkSample(X)

    if( length(d) != 1 || !is.numeric(d) || d < 2 || (d%%1) != 0 ){
      stop("The order d must be an integer number greater than 1.")
    }
    # Any other input will be validated within hdMTD_FS() or hdMTD_CUT() functions

    lenX <- length(X)
    if (lenX <= 2 * (d + 1)) {
      stop("The FSC method splits data in two, therefore the sample size must be greater than 2*(d+1).")
    }

    # Set the state space if not provided
    if(is.null(A)) { A <- sort(unique(X)) } else { A <- sort(A) }

    # Split the data into two parts
    m <- lenX%/%2
    Xm <- X[seq_len(m)]
    Xn <- X[(m + 1):lenX]

    # Apply the two methods sequentially
    S <- hdMTD_FS(Xm, d = d, l = l, A = A)
    S <- hdMTD_CUT(Xn, d = d, S = S, A = A, alpha = alpha, mu = mu, xi = xi)

    S
}

