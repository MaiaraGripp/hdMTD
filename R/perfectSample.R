#' Perfectly samples an MTD Markov chain
#'
#' Samples an MTD Markov Chain from the stationary distribution.
#'
#' @param MTD An MTD object, see [MTDmodel()] for properly generating a MTD object.
#' @param N The sample size. If \code{NULL} sample size will be set to 1000.
#'
#' @return Returns a sample from an MTD model (the first element is the most recent).
#'
#' @details This perfect sample algorithm requires that the MTD model has
#' an independent distribution (p0) with a positive weight (i.e., \code{MTD$lambdas["lam0"]>0} which
#' means \eqn{\lambda_0>0}).
#'
#' @export perfectSample
#'
#' @examples
#' perfectSample(MTDmodel(Lambda=c(1,4), A = c(0,2)), N = 200 )
#' perfectSample(MTDmodel(Lambda=c(2,5), A = c(1,2,3)), N = 1000 )
#'
perfectSample <- function(MTD, N = NULL){
  UseMethod("perfectSample")
}

#' @export
perfectSample.MTD <- function(MTD, N = NULL) {

    # Validate inputs
    checkMTD(MTD)
    if (!is.null(N)) {
        if (!is.numeric(N) || length(N) != 1 || N %% 1 != 0 || N <= max(MTD$Lambda)) {
            stop("N must be NULL or a positive integer greater than max(MTD$Lambda)=",
                 max(MTD$Lambda),".")
        }
    } else {
        warning("Sample size N not provided. N will be set to 1000.")
        N <- 1000
    }

    # Ensure the MTD model is suitable for perfect sampling
    if (sum(MTD$p0) == 0 || MTD$lambdas[1] == 0) {
        stop("This perfect sampling algorithm requires the MTD to have an independent distribution with a positive weight.")
    }

    # Compute cumulative distributions
    cum_lam <- cumsum(MTD$lambdas)
    cum_p0 <- cumsum(MTD$p0)

    chain <- NULL

    # Generate the initial state sequence (length max(MTD$Lambda))
    repeat {
        chain <- c(NA, chain)  # Start the chain with an undefined symbol
        Y <- integer() # Stores positions in the chain that need to be sampled
        t <- 1
        Y[t] <- 1 # First position for which to sample a symbol
        Kt <- 1 # Distance between Y[t] and Y[t-1] (initialized as > 0)

        # Sample a symbol from the independent distribution (i.e. x ~ p0)
        x <- MTD$A[ which(cum_p0 > runif(1))[1] ]

        while (Kt > 0) { # Kt > 0 => Y[t] > Y[t-1], this means there is a position further in
            # the past which is relevant for the present and the algorithm must continue to find
            # a symbol for it ( Kt = 0 iff Y[t] = Y[t-1])

            # Sample the lag that is relevant for sampling chain[Y[t]]
            u <- runif(1)
            Kt <- c(0, MTD$Lambda)[ which(cum_lam > u)[1] ] # K_t ~ lambdas

            t <- t + 1
            Y[t] <- Y[t - 1] + Kt # Y[t] is relevant for sampling chain[Y[t-1]].
            # If Kt = 0, Y[t] = Y[t-1] and the symbol in Y[t-1] does not depend on the past

            # If a symbol already exists at Y[t], use it and stop
            if (!is.na(chain[ Y[t] ])) {
                x <- chain[ Y[t] ]
                Y[t + 1] <- Y[t]
                break
            }
            # If there was nothing previously sampled this loop only stops when Kt = 0

        } # At the end of this loop we have a vector (Y) of positions of the chain for
        # which we must sample symbols.

        # Assign sampled values
        lenY <- length(Y) - 1 # The number of positions, or the number of symbols to be sampled
        chain[Y[lenY]] <- x # Either sampled from p0 or repeated from the chain

        if (lenY > 1)
            {
                for (t in lenY:2) {
                  Kt <- Y[t] - Y[t - 1]  # Kt > 0
                  # Kt = j => chain[Y[t-1]] ~ pj(.|chain[Y[t]]), j\in Lambda
                  p_Kt <- MTD$pj[[ which(MTD$Lambda == Kt) ]]
                  p_KtRow <- which(MTD$A == chain[Y[t]])
                  cum_p_KtRow <- cumsum(p_Kt[p_KtRow, ])
                  chain[Y[t - 1]] <- MTD$A[ which(cum_p_KtRow > runif(1))[1] ]
                }
            }

        if (all(!is.na(chain[ seq_len(max(MTD$Lambda)) ] ))) {break}
        # Stop once the initial sequence is fully sampled
    }

    chain <- chain[seq_len(max(MTD$Lambda))] # Set the initial state

    # Generate the rest of the sample
    for (i in seq_len(N-max(MTD$Lambda))) {
        chain <- c(NA, chain)

        # Samples which lag Kt is relevant for chain[1], Kt \in {0}\cup Lambda
        u <- runif(1)
        Kt <- c(0, MTD$Lambda)[ which(cum_lam > u)[1] ]

        chain[1] <- if (Kt == 0) {
            # Sample from p0, i.e. chain[1]~p0
            MTD$A[which(cum_p0 > runif(1))[1]]
        } else {
            # Kt=j => chain[1]~pj(.|chain[1+Kt]) j \in Lambda
            p_Kt <- MTD$pj[[which(MTD$Lambda == Kt)]]
            p_KtRow <- which(MTD$A == chain[1 + Kt])
            cum_p_KtRow <- cumsum(p_Kt[p_KtRow, ])
            chain[1] <- MTD$A[ which(cum_p_KtRow > runif(1))[1] ]
        }
    }
    chain # Returns X1...XN
}

#' @export
perfectSample.default <- function(MTD,N){
    print("The implemented perfect sample algorithm works only for MTD class objetcs.")
}
