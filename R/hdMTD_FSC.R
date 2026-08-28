#' Forward Stepwise and Cut method for inference in MTD models
#'
#' A function for inference in MTD Markov chains using the FSC method. This function
#' estimates the relevant lag set \eqn{\Lambda} through the FSC algorithm.
#'
#' @param X A vector or single-column data frame containing a chain sample
#' (`X[1]` is the most recent). Missing values (`NA`) are not allowed.
#' @param d A positive integer representing an upper bound for the chain order.
#' @param l An integer between 2 and \code{d} specifying the number of candidate lags
#' selected by the FS step before applying CUT.
#' @param A A vector with distinct nonnegative integers representing the state space.
#' If not informed, this function will set \code{A <- sort(unique(X))}.
#' @param alpha A positive real number used in the CUT threshold (which determines if two
#' distributions can be considered different). The larger the \code{alpha}, the greater
#' the distance required to consider that there is a difference between a set of distributions.
#' Defaulted to 0.05.
#' @param mu A positive real number such that \eqn{\code{mu}>(e^{\code{mu}}-1)/2}. \code{mu}
#' is also a component of the same threshold as \code{alpha}.
#' @param xi A positive real number, \code{xi} is also a component of the same threshold as
#'  \code{alpha}.
#' @param cut_fraction A numeric value strictly between 0 and 1 specifying the proportion
#' of the sample assigned to the CUT step. Since \code{X[1]} is the most recent observation,
#' CUT uses \code{X[1:ceiling(cut_fraction * length(X))]} observations, whereas FS uses the remaining,
#' chronologically older observations. Both subsamples must contain more than \code{d + 1}
#' observations. The default is \code{0.5}.
#' @param ... Additional arguments (not used in this function, but maintained for compatibility with [hdMTD()]).
#'
#' @details The "Forward Stepwise and Cut" (FSC) algorithm combines the
#' Forward Stepwise (FS) and CUT methods. The FS step uses the chronologically
#' older part of the sample to select candidate lags, which are then refined by
#' applying CUT to the most recent part. The proportion of the sample assigned to CUT is
#' controlled by \code{cut_fraction}. The method was developed by [Ost and Takahashi](http://jmlr.org/papers/v24/22-0266.html)
#' and is especially useful for inference in high-order MTD Markov chains.
#'
#' @references
#' Ost, G. & Takahashi, D. Y. (2023).
#' Sparse Markov models for high-dimensional inference.
#' *Journal of Machine Learning Research*, *24*(279), 1-54.
#' \url{http://jmlr.org/papers/v24/22-0266.html}
#'
#' @return Returns a vector with the estimated relevant lag set using FSC algorithm.
#' @export
#' @examples
#' # Simulate a chain from an MTD model
#' set.seed(1)
#' M <- MTDmodel(Lambda = c(1, 3), A = c(1, 2), lam0 = 0.05)
#' X <- perfectSample(M, N = 400)
#'
#' # Forward Stepwise and Cut with different parameters
#' hdMTD_FSC(X, d = 4, l = 2)
#' hdMTD_FSC(X, d = 4, l = 3, alpha = 0.1)
#'
hdMTD_FSC <- function(X, d, l, alpha = 0.05, mu = 1, xi = 0.5, A = NULL, cut_fraction = 0.5, ...){
    # Validate inputs
    X <- checkSample(X)

    if( length(d) != 1 || !is.numeric(d) || d < 2 || (d%%1) != 0 ){
      stop("The order d must be an integer number greater than 1.")
    }

    if ( length(l) != 1 || !is.numeric(l) || is.na(l) ||
         l %% 1 != 0 || l < 2 || l > d) {
      stop("l must be an integer between 2 and d for the FSC method.")
    }

    if (length(cut_fraction) != 1 || !is.numeric(cut_fraction) ||
        is.na(cut_fraction) || !is.finite(cut_fraction) ||
        cut_fraction <= 0 || cut_fraction >= 1) {
      stop("cut_fraction must be a single numeric value strictly between 0 and 1.")
    }
    # Any other input will be validated within hdMTD_FS() or hdMTD_CUT() functions

    lenX <- length(X)

    nCUT <- ceiling(cut_fraction * lenX)

    if (nCUT <= (d + 1) || (lenX - nCUT) <= (d + 1)) {
      stop(
        paste0(
          "The sample split is not valid. Both the FS and CUT subsamples must contain more than d + 1 observations."
        )
      )
    }

    # Set the state space if not provided
    if(is.null(A)) { A <- sort(unique(X)) } else { A <- sort(A) }

    # Split the data into two parts
    X_CUT <- X[seq_len(nCUT)]
    X_FS <- X[(nCUT + 1):lenX]

    # Apply the two methods sequentially
    S <- hdMTD_FS(X_FS, d = d, l = l, A = A)
    S <- hdMTD_CUT(X_CUT, d = d, S = S, A = A, alpha = alpha, mu = mu, xi = xi)

    S
}

