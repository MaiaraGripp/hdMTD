#' Perfectly samples an MTD Markov chain
#'
#' Samples an MTD Markov Chain from the stationary distribution.
#'
#' @name perfectSample
#' @param object An object of class "MTD", see [MTDmodel()] for properly generating
#' @param ... Additional arguments passed to methods.
#' an MTD object.
#'
#' @return Returns a size N sample from an MTD model (the first element is the most recent).
#'
#' @details This perfect sample algorithm requires that the MTD model has
#' an independent distribution (p0) with a positive weight (i.e., \code{lambdas(object)["lam0"]>0}
#' which means \eqn{\lambda_0>0}).
#'
#' @examples
#' M <- MTDmodel(Lambda = c(1, 3, 4), A = c(0, 2))
#' perfectSample(M, N = 200)
#'
#' M <- MTDmodel(Lambda = c(2, 5), A = c(1, 2, 3))
#' perfectSample(M, N = 300)
#'
#' @export
perfectSample <- function(object, ...){
  UseMethod("perfectSample")
}

#' @rdname perfectSample
#' @param N Positive integer. Sample size to generate. Must be > max(Lambda(object)).
#' @export
perfectSample.MTD <- function(object, N, ...) {

  # Validate inputs
  checkMTD(object)
  d <- max(Lambda(object))
  if (missing(N) || length(N) != 1L || !is.numeric(N) || !is.finite(N) || N <= d || N %% 1 != 0) {
    stop(
      paste0(
        "Argument N is missing or invalid. ",
        "The intended sample size N is required and must be a positive integer > max(Lambda(object)) = ",
        d, ". ",
        "Example: perfectSample(object, N = 1000) returns a sample of size 1000."
      )
    )
  }
  N <- as.integer(N)

  # Ensure the MTD model is suitable for perfect sampling
  if (sum(p0(object)) == 0 || lambdas(object)[1] == 0) {
    stop("This perfect sampling algorithm requires the MTD object to have an independent distribution p0 with a positive weight lam0 > 0.")
  }

  # Compute cumulative distributions
  cum_lam <- cumsum(lambdas(object))
  cum_p0 <- cumsum(p0(object))

  chain <- NULL

  # Generate the initial state sequence (length d)
  repeat {
    chain <- c(NA, chain)  # Start the chain with an undefined symbol
    Y <- integer() # Stores positions in the chain that need to be sampled
    t <- 1
    Y[t] <- 1 # First position for which to sample a symbol
    Kt <- 1 # Distance between Y[t] and Y[t-1] (initialized as > 0)

    # Sample a symbol from the independent distribution (i.e. x ~ p0)
    x <- states(object)[ which(cum_p0 > runif(1))[1] ]

    while (Kt > 0) { # Kt > 0 => Y[t] > Y[t-1], this means there is a position further in
      # the past which is relevant for the present and the algorithm must continue to find
      # a symbol for it ( Kt = 0 iff Y[t] = Y[t-1])

      # Sample the lag that is relevant for sampling chain[Y[t]]
      u <- runif(1)
      Kt <- c(0, Lambda(object))[ which(cum_lam > u)[1] ] # K_t ~ lambdas

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
        p_Kt <- pj(object)[[ which(Lambda(object) == Kt) ]]
        p_KtRow <- which(states(object) == chain[Y[t]])
        cum_p_KtRow <- cumsum(p_Kt[p_KtRow, ])
        chain[Y[t - 1]] <- states(object)[ which(cum_p_KtRow > runif(1))[1] ]
      }
    }

    if ( !anyNA( chain[seq_len(d)] ) ) {break}
    # Stop once the initial sequence is fully sampled
  }

  chain <- chain[seq_len(d)] # Set the initial state

  # Generate the rest of the sample
  for (i in seq_len(N - d)) {
    chain <- c(NA, chain)

    # Samples which lag Kt is relevant for chain[1], Kt \in {0}\cup Lambda
    u <- runif(1)
    Kt <- c(0, Lambda(object))[ which(cum_lam > u)[1] ]

    chain[1] <- if (Kt == 0) {
      # Sample from p0, i.e. chain[1]~p0
      states(object)[which(cum_p0 > runif(1))[1]]
    } else {
      # Kt=j => chain[1]~pj(.|chain[1+Kt]) j \in Lambda
      p_Kt <- pj(object)[[which(Lambda(object) == Kt)]]
      p_KtRow <- which(states(object) == chain[1 + Kt])
      cum_p_KtRow <- cumsum(p_Kt[p_KtRow, ])
      chain[1] <- states(object)[ which(cum_p_KtRow > runif(1))[1] ]
    }
  }
  chain
}

