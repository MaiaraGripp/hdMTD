#' Oscillations of an MTD Markov chain
#'
#' Calculates the oscillations of an MTD model object or estimates the oscillations of a chain sample.
#'
#' @param x Must be an MTD object or a chain sample.
#' @param ... Additional parameters that might be required. Such as:
#'
#' \code{S}: If \code{x} is a chain sample the user should provide a set of lags
#' for which he wishes to estimate the oscillations. It must be labeled as \code{S}, an in
#' this scenario the function takes an upper bound for the order as \code{d=max{S}}.
#'
#' \code{A}: If \code{x} is a chain sample, and there may be elements in \code{A} that did not
#' appear in \code{x}, the state space should be specified, and it must be labeled as \code{A}.
#'
#' @details The oscillation for a certain lag \eqn{j} of an MTD model
#' ( \eqn{ \{ \delta_j:\ j \in \Lambda \} }), is the product of the weight \eqn{\lambda_j}
#' multiplied by the maximum of the total variation distance between the distributions in a
#' stochastic matrix \eqn{p_j}.
#' \deqn{\delta_j = \lambda_j\max_{b,c \in \mathcal{A}} d_{TV}(p_j(\cdot | b), p_j(\cdot | c)).}
#' So, if \code{x} is an MTD object, the parameters \eqn{\Lambda}, \eqn{\mathcal{A}}, \eqn{\lambda_j},
#' and \eqn{p_j} are inputted through, respectively, the entries \code{Lambda}, \code{A},
#' \code{lambdas} and the list \code{pj} of stochastic matrices. Hence, an oscillation \eqn{\delta_j}
#' may be calculated for all \eqn{j \in \Lambda}.
#'
#' @details If we wish to estimate the oscillations from a sample, then \code{x} must be a chain,
#' and \code{S}, a vector representing a set of lags, must be informed. This way the transition
#' probabilities can be estimated. Let \eqn{\hat{p}(\cdot| x_S)} symbolize an estimated distribution
#' in \eqn{\mathcal{A}} given a certain past \eqn{x_S} ( which is a sequence of elements of
#'  \eqn{\mathcal{A}} where each element occurred at a lag in \code{S}), and
#'  \eqn{\hat{p}(\cdot|b_jx_S)} an estimated distribution given past \eqn{x_S} and that the symbol
#'  \eqn{b\in\mathcal{A}} occurred at lag \eqn{j}.
#'  If \eqn{N} is the sample size, \eqn{d=}\code{max(S)} and \eqn{N(x_S)} is the
#'  number of times the sequence \eqn{x_S} appeared in the sample, then
#' \deqn{\delta_j = \max_{c_j,b_j \in \mathcal{A}} \frac{1}{N-d}\sum_{x_{S} \in \mathcal{A}^{S}} N(x_S)d_{TV}(\hat{p}(. | b_jx_S), \hat{p}(. | c_jx_S) )}
#' is the estimated oscillation for a lag \eqn{j \in \{1,\dots,d\}\setminus}\code{S}. Note that \eqn{\mathcal{A}^S} is the space of
#' sequences of \eqn{\mathcal{A}} indexed by \code{S}.

#' @return If the \code{x} parameter is an MTD object, it will provide the oscillations for
#' each element in \code{Lambda}. In case \code{x} is a chain sample, it estimates the oscillations
#'  for a user-inputted set \code{S} of lags.
#'
#' @export oscillation
#' @examples
#' oscillation( MTDmodel(Lambda=c(1,4),A=c(2,3) ) )
#' oscillation(MTDmodel(Lambda=c(1,4),A=c(2,3),lam0=0.01,lamj=c(0.49,0.5),
#' pj=list(matrix(c(0.1,0.9,0.9,0.1),ncol=2)), single_matrix=TRUE))
#'
oscillation <- function(x,...) {UseMethod("oscillation")}

#' @export
oscillation.MTD <- function(x,...){

  # Validate input
  checkMTD(x)

  lenA <- length(x$A) # The number of rows in each pj matrix
  rows <- t(utils::combn(lenA, 2)) # All possible pairs of pj rows

  dTV_pj <- function(pj, rows) {
      # In a matrix pj, each row is a distribution. dTV_pj calculates the total
      # variation distance between each pair of distributions and returns the maximum.
      apply(rows, 1, function(r) sum(abs(pj[r[1], ] - pj[r[2], ]))/2) |> max()
  }

  # For each j in Lambda, oscillation_j = \lambda_j * dTV_pj
  y <- x$lambdas[-1] * sapply(x$pj, dTV_pj, rows) # oscillations for each j
  names(y) <- paste0("-", x$Lambda)
  y
}

#' @export
oscillation.default <- function(x,...){

    params <- list(...)

    # Validate inputs
    x <- checkSample(x)
    check_oscillation_inputs(x, params)

    S <- params$S
    A <- params$A

    if (length(A) == 0) {A <- sort(unique(x))} else {A <- sort(A)}

    lenX <- length(x)
    S <- sort(S,decreasing = TRUE)
    lenS <- length(S)
    lenA <- length(A)

    A_pairs <- t(utils::combn(A, 2)) # all possible pairs with elements from A
    lenPairs <- nrow(A_pairs) # number of pairs

    base <- countsTab(X = x, d = S[1])
    b_Sja <- freqTab(S = S, j = NULL, A = A,countsTab = base) # Frequency of sequences in sample

    y <- numeric(lenS) # Initiate vector for oscillations

    for(i in seq_len(lenS)){

        j <- S[i] # j is a lag from S

        if (lenS > 1) {
            Z <- setdiff(S, j) # Z <-  S\{j}
            b_S <- groupTab(S = Z, j = NULL, b_Sja, lenX = lenX, d = S[1])
            PositiveNx_S <- which(b_S$Nx_Sj > 0) # Positions of the x_Z that appeared in the sample
            subx <- b_S[PositiveNx_S, -lenS] # List of these x_Z
        } else {
            Z <- NULL
            PositiveNx_S <- 1
            subx <- matrix(0, ncol = 1)
        }
        lenPositiveNx_S <- length(PositiveNx_S)

        # Initiate matrix to store total variation distances (dTV)
        dtv_xS <- matrix(0, ncol = lenPairs, nrow = lenPositiveNx_S)
        colnames(dtv_xS) <- apply(A_pairs, 1, paste0, collapse = "x")
        rownames(dtv_xS) <- apply(subx, 1, paste0, collapse = "")

        # Compute total variation distances
        for (t in seq_len(lenPositiveNx_S)) {
            dtv_xS[t, ] <- dTV_sample(S = Z, j = j, lenA = lenA, base = b_Sja,
                                        A_pairs = A_pairs, x_S = subx[t, ])
        }

        if (lenS > 1) { # Compute weighted total variation distances (dTV)
            NxS_dtvxS <- sweep(dtv_xS, MARGIN = 1, t(b_S[PositiveNx_S, lenS]), '*')
            h_deltajbc <- colSums(NxS_dtvxS)/(lenX - S[1]) # weighted dTV
            # h_deltajbc = \sum_{x_Z:N(x_Z)>0} dTV[\hat{p}(.|x_Zb_j),\hat{p}(.|x_Zc_j)] * N(x_Z)/(n-max(S)) for all b!=c\in A
        } else {
            h_deltajbc <- dtv_xS # Z=NULL, so no weighting by N(x_Z) is needed
        }

        y[i] <- max(h_deltajbc)
    }
    names(y) <- paste0("-", S)
    y <- rev(y)
    y
}


