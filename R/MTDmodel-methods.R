#' Methods for objects of class \code{"MTD"}
#'
#' @description
#' Printing, summarizing, and coefficient-extraction methods for Mixture
#' Transition Distribution (MTD) model objects.
#'
#' @details
#' \code{print.MTD()} displays the relevant lag set (shown as negative integers)
#' and the state space. For a detailed overview including mixture weights and a
#' compact preview of the global transition matrix \eqn{P}, use \code{summary()}.
#'
#' @param x An object of class \code{"MTD"} or \code{"summary.MTD"},
#'  depending on the method.
#' @param object An object of class \code{"MTD"}.
#' @param X A vector or single-column data frame containing an MTD chain sample
#' (values must be in the model's state space).
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return
#' \describe{
#'   \item{\code{print.MTD}}{Invisibly returns the \code{"MTD"} object, after
#'         displaying its relevant lag set and state space.}
#'   \item{\code{summary.MTD}}{An object of class \code{"summary.MTD"} with fields:
#'         \code{order}, \code{states}, \code{lags}, \code{indep},
#'         \code{lambdas}, \code{p0} (or \code{NULL}),
#'         \code{P_dim}, and \code{P_head}.}
#'   \item{\code{print.summary.MTD}}{Invisibly returns the
#'         \code{"summary.MTD"} object after printing its contents.}
#'   \item{\code{coef.MTD}}{A list with model parameters:
#'         \code{lambdas}, \code{pj}, and \code{p0}.}
#'   \item{\code{logLik.MTD}}{An object of class \code{"logLik"} with attributes
#'     \code{nobs} (number of transitions) and \code{df} (free parameters),
#'     honoring model constraints such as \code{single_matrix} and the independent
#'     component (\code{indep_part}).}
#' }
#'
#' @seealso
#' \code{\link{MTDmodel}}, \code{\link{MTDest}},
#' \code{\link{transitP}}, \code{\link{lambdas}}, \code{\link{pj}},
#' \code{\link{p0}}, \code{\link{lags}}, \code{\link{Lambda}}, \code{\link{states}},
#' \code{\link{summary.MTDest}}, \code{\link{coef.MTDest}},
#' \code{\link{oscillation}}, \code{\link{perfectSample}},
#' \code{\link[stats]{logLik}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' m <- MTDmodel(Lambda = c(1, 3), A = c(0, 1), lam0 = 0.05)
#'
#' print(m)       # compact display: lags (Z^-) and state space
#' s <- summary(m)
#' print(s)
#'
#' coef(m)        # list(lambdas = ..., pj = ..., p0 = ...)
#' transitP(m)    # global transition matrix P
#' pj(m); p0(m); lambdas(m); lags(m); Lambda(m); states(m)
#' }
#'
#' @name MTD-methods
NULL

# --------------------------- print.MTD ---------------------------------

#' @rdname MTD-methods
#' @exportS3Method print MTD
print.MTD <- function(x, ...) {
  lg <- lags(x)
  A  <- states(x)

  cat("An object of class 'MTD'\n")
  cat("  Relevant lags: ", fmt_vec(lg), "\n", sep = "")
  cat("  State space (A): ", fmt_vec(A),  "\n", sep = "")
  cat("  Use summary(x) for full description.\n")
  cat("  Accessors: transitP(), lambdas(), pj(), p0(), lags(), Lambda(), states().\n")
  cat("  Methods:   coef(), probs(), oscillation(), perfectSample(), logLik().\n")
  invisible(x)
}

# ------------------------- summary.MTD ---------------------------------

#' @rdname MTD-methods
#' @exportS3Method summary MTD
summary.MTD <- function(object, ...) {
  checkMTD(object)  # robust validation

  w   <- lambdas(object)
  p0 <- p0(object)
  indep_flag <- (sum(p0) > 0) && (w[1] > 0)
  P <- transitP(object)

  out <- list(
    call    = if (!is.null(object$call)) object$call else NULL,
    order   = max(Lambda(object)),
    Lambda  = Lambda(object),
    states  = states(object),
    lags    = lags(object),
    indep   = indep_flag,
    lambdas = w,
    pj      = pj(object),
    p0      = if (indep_flag) p0 else NULL,
    P_dim   = dim(P),
    P       = P
  )
  class(out) <- "summary.MTD"
  out
}

#' @rdname MTD-methods
#' @exportS3Method print summary.MTD
print.summary.MTD <- function(x, ...) {
  cat("Mixture Transition Distribution (MTD) model \n")
  if (!is.null(x$call)) { cat("\nCall:\n"); print(x$call) }

  cat("\nRelevant lags: ", fmt_vec(x$lags), "\n", sep = "")
  cat("State space: ", fmt_vec(x$states), "\n", sep = "")

  cat("\nlambdas (weights):\n"); print(x$lambdas)

  if (!is.null(x$p0)) {
    cat("\nIndependent distribution p0:\n"); print(x$p0)
  }

  cat("\nTransition matrices pj (one per lag):\n")
  for (i in seq_along(x$pj)) {
    cat(sprintf(" \n pj for lag j = %s:\n", -x$Lambda[i]))
    print(x$pj[[i]])
  }

  cat(sprintf("\nTransition matrix P: %d x %d\n", x$P_dim[1], x$P_dim[2]))
  cat("- Preview of first rows of P:\n"); print(utils::head(x$P, n = min(6L, nrow(x$P))))

  ## ---- Reading guide for P (right-to-left interpretation) ----

  cat("\nReading guide for P:\n")
  if (!is.null(x$Lambda)) {
    cat("Rows list past contexts from oldest to newest, matching lags ",
        paste0("(", paste(-sort(x$Lambda, decreasing = TRUE), collapse = ", "), ")"), ".\n", sep = "")
  } else {
    cat("Rows list past contexts from oldest to newest.\n")
  }
  invisible(x)
}

# --------------------------- coef.MTD ----------------------------------

#' @rdname MTD-methods
#' @exportS3Method coef MTD
coef.MTD <- function(object, ...) {
  checkMTD(object)
  list(
    lambdas = lambdas(object),
    pj      = pj(object),
    p0      = p0(object)
  )
}

# --------------------------- logLik.MTD ----------------------------------

#' @rdname MTD-methods
#' @exportS3Method logLik MTD
logLik.MTD <- function(object, X,...) {
  checkMTD(object)
  X <- checkSample(X)

  L <- Lambda(object)
  d <- max(L)
  if (length(X) <= d) stop("Insufficient sample size: length(X) must be > max(Lambda(object)).")

  A <- states(object)
  if (!all(X %in% A)) {
    stop("Sample contains values outside the model state space.")
  }

  # Sufficient statistics from data
  ct <- countsTab(X, d)
  ft <- freqTab(S = L, j = NULL, A =  A, countsTab = ct)
  pos <- which(ft$Nxa_Sj > 0)

  P <- transitP(object)
  ll <- sum(log(as.vector(t(P))[pos]) * ft$Nxa_Sj[pos])

  nobs <- length(X) - d
  indep_part <- lambdas(object)[1] > 0
  single_matrix <- isTRUE(object$single_matrix)

  df <- n_parameters(Lambda = L, A = A,
                     single_matrix = single_matrix,
                     indep_part = indep_part)

  structure(ll, nobs = nobs, df = df, class = "logLik")
}
