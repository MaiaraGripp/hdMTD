#' The Forward Stepwise (FS) method for inference in MTD models
#'
#'  A function that estimates the set of relevant lags of an MTD model using the FS method.
#'
#' @param X A vector or single-column data frame containing a chain sample (`X[1]` is the most recent).
#' @param d A positive integer representing an upper bound for the chain order.
#' @param l A positive integer specifying the number of lags to be selected as relevant.
#' @param A A vector with positive integers representing the state space. If not informed,
#' this function will set \code{A <- sort(unique(X))}.
#' @param elbowTest Logical. If TRUE, the function applies an alternative stopping criterion to
#' determine the length of the set of relevant lags. See *Details* for more information.
#' @param warn Logical. If \code{TRUE}, the function warns the user when \code{A} is set automatically.
#' @param ... Additional arguments (not used in this function, but maintained for compatibility with [hdMTD()].
#'
#'
#' @details The "Forward Stepwise" (FS) algorithm is the first step of the "Forward Stepwise and Cut"
#' (FSC) algorithm for inference in Mixture Transition Distribution (MTD) models.
#' This method was developed by [Ost and Takahashi](http://jmlr.org/papers/v24/22-0266.html)
#' This specific function will only apply the FS step of the algorithm and return an estimated
#' relevant lag set of length \code{l}.
#'
#' This method iteratively selects the most relevant lags based on a certain quantity \eqn{\nu}.
#' In the first step, the lag in \code{1:d} with the greatest \eqn{\nu} is deemed important.
#' This lag is included in the output, and using this knowledge, the function proceeds to seek
#' the next important lag (the one with the highest \eqn{\nu} among the remaining ones).
#' The process stops when the output vector reaches length \code{l} if \code{elbowTest=FALSE}.
#'
#' If \code{elbowTest = TRUE}, the function will store these maximum \eqn{\nu} values
#' at each iteration, and output only the lags that appear before the one with smallest
#' \eqn{\nu} among them.
#'
#' @note Tie-breaking: if multiple lags share the maximum ν, FS picks the most recent
#' lag (smallest j). Up to version 0.1.2 ties were broken at random, which could cause
#' run-to-run differences.
#'
#' @references
#' Ost, G. & Takahashi, D. Y. (2023).
#' Sparse Markov models for high-dimensional inference.
#' *Journal of Machine Learning Research*, *24*(279), 1-54.
#' \url{http://jmlr.org/papers/v24/22-0266.html}
#'
#' @return A numeric vector containing the estimated relevant lag set using FS algorithm.
#'
#' @importFrom utils combn
#' @importFrom utils tail
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' # Simulate a chain from an MTD model
#' set.seed(1)
#' M <- MTDmodel(Lambda = c(1, 3), A = c(1, 2), lam0 = 0.05)
#' X <- perfectSample(M, N = 400)
#'
#' # Forward Stepwise with l = 2
#' hdMTD_FS(X, d = 5, l = 2)
#'
#' # Forward Stepwise with l = 3
#' hdMTD_FS(X, d = 4, l = 3)
#'
hdMTD_FS <- function(X, d, l, A = NULL, elbowTest = FALSE, warn = FALSE,...){

    # Validate and preprocess the input sample
    X <- checkSample(X)
    check_hdMTD_FS_inputs(X, d, l, A, elbowTest, warn) # See validation.R.

    # Set the state space if not provided
    if(length(A) == 0) { A <- sort(unique(X)) } else { A <- sort(A) }

    # Try to generate all possible sequences of length l with elements of A
    if(inherits(try(expand.grid(rep(list(A), l)),silent = TRUE),"try-error")) {
      stop(paste0("The dataset with all sequences of length l is too large. Try reducing l."))
    }

    lenA <- length(A)
    lenX <- length(X)
    base <- countsTab(X = X, d = d)

    A_pairs <- t(utils::combn(A, 2)) # All unique state pairs
    A_pairsPos <- t(utils::combn(seq_len(lenA), 2)) # Corresponding index pairs
    nrowA_pairs <- nrow(A_pairs) # number of pairs

    S <- NULL # Set of selected lags
    lenS <- 0
    maxnu <- numeric(l) # Stores max ν values (used if elbowTest is TRUE)

    while (lenS < l) {

      Sc <- setdiff(seq_len(d), S) # Available lags to select
      Sc <- sort(Sc, decreasing = TRUE)

      dec_S <- sort(S, decreasing = TRUE)
      nuj <- numeric(length(Sc)) # stores nu values for each lag j

      for (z in seq_along(Sc)) {
          j <- Sc[z] # choose one of the available lags

          # Tables with frequencies given pasts in the selected lag set S
          b_Sja <- freqTab(S = S, j = j, A = A, countsTab = base) # given lag j and present state
          b_Sj <- groupTab(S = S, j = j, b_Sja, lenX = lenX, d = d) # given j but without present state
          b_S <- groupTab(S = S, j = NULL, b_Sja, lenX = lenX, d = d) # if S is NULL bs <- [0,lenX-d]
          # groupTab() is defined at utils.R
          ncolb_S <- ncol(b_S)

          # Identify sequences x_S with N(x_S) > 0 (i.e. that appeared in the sample)
          if(lenS > 0) {
            PositNx_S <- which(b_S$Nx_Sj > 0) # Position of x_S with N(x_S) > 0
            subx <- b_S[, -ncolb_S] # all possible x_S
          }else{
            PositNx_S <- 1
            subx <- matrix(0, ncol = 1) # Required format
          }

          # Compute νj
          for (t in PositNx_S) { # runs in all sequences x_S with N(x_S) > 0

            # Compute frequencies given x_S (if lenS>0) and all possible
            # symbols of A at lag j
            PIs <- PI(S = dec_S, groupTab = b_Sj, x_S = subx[t, ],
                      lenX = lenX, d = d) # Must use S = dec_S to match order of symbols in x_S
            # PI() is defined in utils.R

            # Compute total variation distances given x_S
            dTVs <- dTV_sample(S = dec_S, j = j, lenA = lenA, base = b_Sja,
                               A_pairs = A_pairs, x_S = subx[t, ])

            cont <- 0
            for (y in seq_len(nrowA_pairs)) {
              cont <- cont + prod(PIs[A_pairsPos[y, ]]) * dTVs[y]
            } # cont = \sum_{b\in A}\sum_{c\in A} P(x_Sb_j)*P(x_Sc_j)*dTV[q(|xSbj)||q(|xScj)]

            PI_xS <- as.numeric(b_S[t, ncolb_S]/(lenX - d)) # If S is NULL, this
            # evaluates to 1
            nuj[z] <- nuj[z] + cont/PI_xS
          }
      }

      maxnu_val <- max(nuj)
      maxnu[lenS + 1] <- maxnu_val # Store maximum ν, used if elbowTest is TRUE.
      posMaxnu <- utils::tail(which(nuj == maxnu_val), 1)  # Position of most recent max ν lag.

      s <- Sc[posMaxnu]
      S <- c(S, s)
      lenS <- length(S)
    }
    if(elbowTest) {
       stp <- max(which(maxnu == min(maxnu)) - 1, 1)
       S <- S[seq_len(stp)]
    }

    S
}

