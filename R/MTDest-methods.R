#' Methods for objects of class \code{"MTDest"}
#'
#' @description
#' Printing method for EM fits of Mixture Transition Distribution (MTD) models.
#'
#' @details
#' The \code{print.MTDest()} method displays a compact summary of the fitted model:
#' the lag set (\code{S}), the state space (\code{A}), the final log-likelihood,
#' and, if available, the number of EM updates performed.
#'
#' @param x An object of class \code{"MTDest"} or \code{"summary.MTDest"},
#'  depending on the method.
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
#' @seealso \code{\link{MTDest}}, \code{\link{summary.MTDest}},
#'   \code{\link{print.summary.MTDest}}, \code{\link{coef.MTDest}},
#'    \code{\link{logLik}}, \code{\link{AIC}}, \code{\link{BIC}}
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
#' @rdname MTDest-methods
#' @exportS3Method print MTDest
print.MTDest <- function(x, ...) {
  cat("An object of class 'MTDest' (EM estimation of MTD model)\n")
  cat("  Lags (-S):", paste(-x$S, collapse = ", "), "\n")
  cat("  State space (A):", paste(x$A, collapse = ", "), "\n")
  cat("  Log-likelihood:", format(x$logLik, digits = 6), "\n")
  if (!is.null(x$iterations)) {
    cat("  Number of updates:", x$iterations, "\n")
  }

  invisible(x)
}
#' @description Summary method for EM fits of MTD models.
#'
#' @details The \code{summary.MTDest()} method collects key fields from an
#' \code{"MTDest"} object into a compact list (class \code{"summary.MTDest"})
#' suitable for printing.
#'
#' @param object An object of class \code{"MTDest"} (used by methods that operate
#' on the fitted model).
#'
#' @rdname MTDest-methods
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

#' @description Printing method for \code{"summary.MTDest"} objects.
#'
#' @details
#' The \code{print.summary.MTDest()} method prints that summary in a readable format,
#' including lambdas, transition matrices, independent distribution,
#' log-likelihood, oscillations (if available), and iteration diagnostics (if available).
#'
#' @rdname MTDest-methods
#' @exportS3Method print summary.MTDest
print.summary.MTDest <- function(x, ...) {
  cat("Summary of EM estimation for MTD model:\n")
  if (!is.null(x$call)) cat("Call: "); if (!is.null(x$call)) print(x$call)

  cat("\nlambdas (weights):\n")
  print(x$lambdas)

  cat("\nTransition matrices pj (one per lag):\n")
  for (i in seq_along(x$pj)) {
    cat(sprintf(" \n pj for lag j = %s:\n", -x$S[i]))
    print(x$pj[[i]])
  }

  if (!is.null(x$p0)) {
    cat("\nIndependent distribution p0:\n")
    print(x$p0)
  } else {
    cat("\nIndependent distribution p0: (none)\n")
  }

  cat("\nLog-likelihood:", format(x$logLik, digits = 6), "\n")

  if (!is.null(x$oscillations)) {
    cat("\nOscillations:\n")
    print(x$oscillations)
  }

  cat("\nLags (-S):", paste(-x$S, collapse = ", "),
      "\nState space (A):", paste(x$A, collapse = ", "), "\n")

  if (!is.na(x$iterations) || !is.na(x$lastComputedDelta)) {
    cat("\nIterations Report:\n")
    cat("Number of updates:", x$iterations, "\n")
    cat("Last compared difference of logLik:", format(x$lastComputedDelta, digits = 6), "\n")
  }

  invisible(x)
}

#' @description Extract coefficients from an \code{"MTDest"} fit.
#'
#' @details The \code{coef.MTDest()} method returns the fitted mixture weights
#' (\code{lambdas}), the list of transition matrices (\code{pj}), and (if present)
#' the independent distribution \code{p0}.
#'
#' @rdname MTDest-methods
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

#' @description Extract log-likelihood from an \code{"MTDest"} fit.
#'
#' @details
#' The \code{logLik.MTDest()} computes the log-likelihood and returns an object
#' of class \code{"logLik"} with attributes: \code{nobs}, the effective sample
#' size used for estimation, and \code{df}, number of free parameters estimated
#' supposing all transition matrices pj to be distinct (multimatrix model).
#' Note: in the returned "logLik" object, \code{df} denotes the number of free
#' parameters (not residual degrees of freedom).
#'
#' @rdname MTDest-methods
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

