#' Methods for objects of class \code{"MTDest"}
#'
#' @description
#' Printing, summarizing, and extracting information from EM fits of Mixture
#'  Transition Distribution (MTD) models.
#'
#' @details
#' These methods handle objects returned by \code{\link{MTDest}} (class \code{"MTDest"}):
#' \itemize{
#'   \item \code{print.MTDest()} displays a compact summary of the fitted model:
#'    the lag set (\code{S}), the state space (\code{A}), the final
#'     log-likelihood, and, if available, the number of EM updates performed.
#'
#'   \item \code{summary.MTDest()} collects key components of the object into a
#'   readable summary list (class \code{"summary.MTDest"}), including lambdas,
#'   transition matrices, independent distribution (if present), log-likelihood,
#'    oscillations (if available), and iteration diagnostics.
#'
#'   \item \code{print.summary.MTDest()} prints the summary in a readable format,
#'    including lambdas, transition matrices, independent distribution,
#'     log-likelihood, oscillations (if available), and iteration diagnostics
#'      (if available).
#'
#'   \item \code{coef.MTDest()} extracts the estimated mixture weights
#'    (\code{lambdas}), the list of transition matrices (\code{pj}), and the
#'     independent distribution (\code{p0}).
#'
#'   \item \code{logLik.MTDest()} returns the log-likelihood as an object of class
#'    \code{"logLik"}, with attributes \code{df} (number of free parameters
#'     under the multimatrix model) and \code{nobs} (effective sample size).
#' }
#'
#' @param x An object of class \code{"MTDest"} or \code{"summary.MTDest"},
#'  depending on the method.
#' @param object An object of class \code{"MTDest"} (used by summary, coef, and
#'  logLik methods).
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return
#' \describe{
#'  \item{\code{print.MTDest}}{Invisibly returns the \code{"MTDest"} object, after
#'     displaying its lag set, state space, final log-likelihood, and iteration
#'     count (if available).}
#'   \item{\code{print.summary.MTDest}}{Invisibly returns the \code{"summary.MTDest"} object,
#'     after displaying its contents: lambdas; transition matrices; independent
#'     distribution; log-likelihood; oscillations (if available); and
#'     iteration diagnostics (if available).}
#'   \item{\code{summary.MTDest}}{An object of class \code{"summary.MTDest"}.}
#'   \item{\code{coef.MTDest}}{A list with estimated \code{lambdas}, \code{pj}, and \code{p0}.}
#'   \item{\code{logLik.MTDest}}{ An object of class \code{"logLik"} with attributes
#'     \code{df} (number of free parameters) and \code{nobs} (effective sample size).}
#' }
#'
#' @seealso \code{\link{MTDest}}, \code{\link{as.MTD}}
#'
#' @examples
#' set.seed(1)
#' MTD <- MTDmodel(Lambda = c(1, 3), A = c(0, 1), lam0 = 0.01)
#' X <- perfectSample(MTD, N = 200)  # small N to keep examples fast
#' init <- list(
#'   p0 = c(0.4, 0.6),
#'   lambdas = c(0.05, 0.45, 0.5),
#'   pj = list(
#'     matrix(c(0.2, 0.8, 0.45, 0.55), byrow = TRUE, ncol = 2),
#'     matrix(c(0.25, 0.75, 0.3, 0.7),  byrow = TRUE, ncol = 2)
#'   )
#' )
#' fit <- MTDest(X, S = c(1, 3), init = init, iter = TRUE)
#' print(fit)
#' summary(fit)
#' coef(fit)
#' logLik(fit)
#' BIC(fit)
#'
#' @name MTDest-methods
NULL

#' @exportS3Method print MTDest
print.MTDest <- function(x, ...) {
  cat("An object of class 'MTDest' (EM estimation of MTD model)\n")
  cat("  Lags (-S):", paste(-x$S, collapse = ", "), "\n")
  cat("  State space (A):", paste(x$A, collapse = ", "), "\n")
  cat("  Log-likelihood:", format(x$logLik, digits = 6), "\n")
  if (!is.null(x$iterations)) {
    cat("  Number of updates:", x$iterations, "\n")
  }
  cat("  Use summary() for full description.\n")
  cat("  Accessors: lambdas(), pj(), p0(), lags(), S(), states().\n")
  cat("  Methods:   coef(), probs(), as.MTD(), logLik().\n")
  invisible(x)
}

#' @exportS3Method summary MTDest
summary.MTDest <- function(object, ...) {
  stopifnot(inherits(object, "MTDest"))

  lenA <- length(object$A)
  lenS <- length(object$S)
  indep <- any(object$p0 != 0)

  out <- list(
    call  = object$call,
    S     = object$S,
    A     = object$A,
    lambdas = object$lambdas,
    pj      = object$pj,
    p0      = if (indep) object$p0 else NULL,
    logLik  = object$logLik,
    oscillations = object$oscillations,
    iterations  = if (!is.null(object$iterations)) object$iterations else NA_integer_,
    lastComputedDelta = if (!is.null(object$lastComputedDelta)) object$lastComputedDelta else NA_real_,
    deltaLogLik = if (!is.null(object$deltaLogLik)) object$deltaLogLik else NA_real_
  )
  class(out) <- "summary.MTDest"
  out
}

#' @exportS3Method print summary.MTDest
print.summary.MTDest <- function(x, ...) {
  cat("Summary of EM estimation for MTD model:\n")

  if (!is.null(x$call)) cat("\nCall:\n"); if (!is.null(x$call)) print(x$call)

  cat("\nLags (-S):", paste(-x$S, collapse = ", "),
      "\nState space (A):", paste(x$A, collapse = ", "), "\n")

  cat("\nlambdas (weights):\n")
  print(x$lambdas)

  if (!is.null(x$p0)) {
    cat("\nIndependent distribution p0:\n")
    print(x$p0)
  } else {
    cat("\nIndependent distribution p0: (none)\n")
  }

  cat("\nTransition matrices pj (one per lag):\n")
  for (i in seq_along(x$pj)) {
    cat(sprintf(" \n pj for lag j = %s:\n", -x$S[i]))
    print(x$pj[[i]])
  }

  cat("\nLog-likelihood:", format(x$logLik, digits = 6), "\n")

  if (!is.null(x$oscillations)) {
    cat("\nOscillations:\n")
    print(x$oscillations)
  }

  if (!is.na(x$iterations) || !is.na(x$lastComputedDelta)) {
    cat("\nIterations Report:\n")
    cat("Number of updates:", x$iterations, "\n")
    cat("Last compared difference of logLik:", format(x$lastComputedDelta, digits = 6), "\n")
  }

  invisible(x)
}

#' @exportS3Method coef MTDest
coef.MTDest <- function(object, ...) {
  stopifnot(inherits(object, "MTDest"))
  out <- list(
    lambdas = object$lambdas,
    pj      = object$pj,
    p0      = object$p0
  )
  out
}

#' @exportS3Method logLik MTDest
logLik.MTDest <- function(object, ...) {
  stopifnot(inherits(object, "MTDest"))

  indep <- any(object$p0 != 0)

  df <- n_parameters(
    Lambda        = object$S,
    A             = object$A,
    single_matrix = FALSE,
    indep_part    = indep,
    zeta          = length(object$S)
  )

  n_eff <- if (!is.null(object$n_eff)) object$n_eff else NA_integer_

  structure(
    as.numeric(object$logLik),
    df   = as.integer(df),
    nobs = as.integer(n_eff),
    class = "logLik"
  )
}

