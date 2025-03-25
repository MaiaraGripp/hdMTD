#' The Cut method for inference in MTD models
#'
#' A function that estimates the set of relevant lags of an MTD model using the CUT method.
#'
#' @param X A vector or single-column data frame containing a chain sample (`X[1]` is the most recent).
#' @param d A positive integer representing an upper bound for the chain order.
#' @param S A numeric vector of distinct positive integers from which this function will select
#' a set of relevant lags. Should be a subset of \code{1:d}. Default is \code{1:d}.
#' @param A A vector with positive integers representing the state space. If not informed,
#' this function will set \code{A <- sort(unique(X))}.
#' @param alpha A positive real number used in the CUT threshold (which determines if two
#' distributions can be considered different). The larger the \code{alpha}, the greater
#' the distance required to consider that there is a difference between a set of distributions.
#' @param mu A positive real number such that \eqn{\code{mu}>(e^{\code{mu}}-1)/2}. \code{mu}
#' is also a component of the same threshold as \code{alpha}.
#' @param xi A positive real number, \code{xi} is also a component of the same threshold as
#'  \code{alpha}.
#' @param warning Logical. If \code{TRUE}, the function warns the user when \code{A} is set automatically.
#' @param ... Additional arguments (not used in this function, but maintained for compatibility with [hdMTD()].
#'
#' @details The "Forward Stepwise and Cut" (FSC) is an algorithm for inference in
#' Mixture Transition Distribution (MTD) models. It consists
#' in the application of the "Forward Stepwise" (FS) step followed by the CUT algorithm.
#' This method and its steps where developed by [Ost and Takahashi](http://jmlr.org/papers/v24/22-0266.html)
#' and are specially useful for inference in high-order MTD Markov chains. This specific function
#' will only apply the CUT step of the algorithm and return an estimated relevant lag set.
#'
#' @references
#' Ost, G. and Takahashi, D. Y. (2023), "Sparse markov models for high-dimensional inference",
#' Journal of Machine Learning Research 24(279), 1-54. [URL: http://jmlr.org/papers/v24/22-0266.html](http://jmlr.org/papers/v24/22-0266.html)
#'
#' @return Returns a set of relevant lags estimated using the CUT algorithm.
#' @export
#' @examples
#' X <- testChains[,3]
#' hdMTD_CUT(X,4,alpha=0.02,mu=1,xi=0.4)
#' hdMTD_CUT(X,d=6,S=c(1,4,6),alpha=0.0065)
#'
hdMTD_CUT <- function(X, d, S=1:d, alpha=0.05, mu=1, xi=0.5, A=NULL, warning=FALSE,...){

    # Validate and preprocess the input sample
    X <- checkSample(X)
    check_hdMTD_CUT_inputs(X, d, S, alpha, mu, xi, A, warning)

    # Set the state space if not provided
    if(length(A) == 0) { A <- sort(unique(X)) } else { A <- sort(A) }

    lenA <- length(A)
    lenS <- length(S)
    dec_S <- sort(S, decreasing = TRUE)

    # Generate all possible past sequences of length |S|-1
    subx <- as.matrix(expand.grid(rep(list(A), lenS - 1))[, (lenS - 1):1],
                      ncol = lenS - 1)
    nrow_subx <- nrow(subx)

    # Generate pairs of states for comparison
    A_pairs <- t(utils::combn(A, 2))
    A_pairsPos <- t(utils::combn(seq_len(lenA), 2))
    nrowA_pairs <- nrow(A_pairs)

    # Compute frequency tables
    base <- countsTab(X = X, d = d)
    b_Sja <- freqTab(S = S, A = A, countsTab = base)

    # Compute total variation distances (dTVs) and thresholds
    dTV_txy <- numeric(lenS)
    for (z in seq_len(lenS)) {
      j <- dec_S[z] # selects j (a lag from S)
      Sminusj <- dec_S[dec_S != j]

      Q <- matrix(0,ncol=nrowA_pairs,nrow = nrow_subx)
      R <- matrix(0,ncol = lenA, nrow = nrow_subx)

      for (k in seq_len(nrow_subx)) {
        # Each k refers to a  different sequence of elements from A indexed
        # by S\{j} (x_S\j)
        Q[k, ] <- dTV_sample(S = Sminusj, j = j, lenA = lenA, base = b_Sja,
                            A_pairs = A_pairs, x_S = subx[k, ])
        # Computes dTVs between distributions ( given sequences x_S where the
        # symbol at lag j varies)

        R[k, ] <- sx(S = Sminusj, freqTab = b_Sja, lenA = lenA, x_S = subx[k,],
                     mu = mu, alpha = alpha, xi = xi)
        # sx is a quantity used to calculate thresholds. See function sx() in utils.R.
      }

      txy <- matrix(0, nrow = nrow(R), ncol = nrow(A_pairs))
      for (s in seq_len(nrowA_pairs)) {
          txy[, s] <- rowSums(R[, A_pairsPos[s, ]])
          # txy = sx + sy is the threshold for comparing distributions
          # conditioned in x_S and in y_S ( x_S and y_S are equal sequences
          # except for the symbol in lag j)
      }

      dTV_txy[z] <- max(Q - txy) # The largest dTV minus its threshold referent to lag j
    }

    dec_S[dTV_txy > 0] # Only the lags where the dTV surpasses the threshold remain
}



