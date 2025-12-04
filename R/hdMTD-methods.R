#' Methods for objects of class \code{"hdMTD"}
#'
#' @description
#' Printing and summarizing methods for lag-selection results returned by
#' \code{\link{hdMTD}}.
#'
#' @details
#' An object of class \code{"hdMTD"} is an integer vector \eqn{S} (selected lags,
#' elements of \eqn{\mathbb{N}^+}) with attributes:
#' \itemize{
#'   \item \code{method}: one of \code{"FS"}, \code{"FSC"}, \code{"CUT"}, \code{"BIC"}.
#'   \item \code{d}: upper bound for the order used in the call.
#'   \item \code{call}: the matched call that produced the object.
#'   \item \code{settings}: a (method-specific) list of the arguments actually used.
#'   \item \code{A}: the state space actually used (provided or inferred).
#'   \item (optional) \code{BIC_selected_value} for \code{method="BIC"} with \code{BICvalue=TRUE}.
#'   \item (optional) \code{extras$BIC_out} for \code{method="BIC"} (exactly the
#'  output of \code{hdMTD_BIC()}).
#' }
#' \code{print()} shows the method, \code{d}, and the selected set of lags in
#' \eqn{\mathbb{N}^+}. \code{summary()} prints the call, the estimated lag
#' set, optional BIC diagnostics and (optionally) the method-specific settings
#' when \code{settings = TRUE}.
#'
#' @param x An object of class \code{"hdMTD"} used in \code{print.hdMTD()}.
#' @param object An object of class \code{"hdMTD"} used in \code{summary.hdMTD()}.
#' @param settings Logical (\code{summary.hdMTD()} only). If \code{TRUE}, the
#'  printed summary includes the method-specific \code{settings} list.
#'  Default \code{FALSE}.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return
#' \describe{
#'   \item{\code{print.hdMTD}}{Invisibly returns the \code{"hdMTD"} object.}
#'   \item{\code{summary.hdMTD}}{Invisibly returns a named list with fields:
#'         \code{call}, \code{S}, \code{lags}, \code{A},
#'         \code{method}, \code{d}, \code{settings}, \code{BIC_selected} and
#'         \code{BIC_out}. Relevant information is printed to the console in
#'          a readable format.}
#' }
#'
#' @seealso \code{\link{hdMTD}}, \code{\link{hdMTD_FS}}, \code{\link{hdMTD_FSC}},
#'   \code{\link{hdMTD_CUT}}, \code{\link{hdMTD_BIC}}, \code{\link{S}}, \code{\link{lags}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' M <- MTDmodel(Lambda = c(1, 4), A = c(1, 3), lam0 = 0.05)
#' X <- perfectSample(M, N = 400)
#' S_hat <- hdMTD(X, d = 5, method = "FS", l = 2)
#' print(S_hat)
#' summary(S_hat)
#' S(S_hat); lags(S_hat)
#' }
#'
#' @name hdMTD-methods
NULL

#' @exportS3Method print hdMTD
print.hdMTD <- function(x, ...) {
  Spos   <- S(x)
  method <- attr(x, "method", exact = TRUE)
  d      <- attr(x, "d",      exact = TRUE)

  cat("An object of class 'hdMTD' (lag selection)\n")
  cat("  Method: ", method, "\n", sep = "")
  cat("  d:      ", d,      "\n", sep = "")
  cat("  Selected S set: ", fmt_vec(Spos),   "\n", sep = "")
  cat("  Accessors: S() and/or lags().\n")
  cat("  Use summary(x, settings = TRUE) for full description.\n\n")
  invisible(x)
}

# ------------------------- summary.hdMTD ---------------------------------

#' @exportS3Method summary hdMTD
summary.hdMTD <- function(object, settings = FALSE,...) {

  # Extract optional BIC info if present
  extras <- attr(object, "extras", exact = TRUE)
  BIC_out <- if (!is.null(extras) && !is.null(extras$BIC_out)) extras$BIC_out else NULL
  BIC_sel    <- attr(object, "BIC_selected_value", exact = TRUE)

  out <- list(
    S        = S(object),
    lags     = lags(object),
    A        = attr(object, "A", exact = TRUE),
    method   = attr(object, "method",   exact = TRUE),
    d        = attr(object, "d",        exact = TRUE),
    settings = attr(object, "settings", exact = TRUE),
    call     = attr(object, "call",     exact = TRUE),
    BIC_selected = BIC_sel,
    BIC_out      = BIC_out
  )
  attr(out, "include_settings") <- isTRUE(settings)  # control printing
  print_hdMTD_summary(out)
  invisible(out)
}

#' @keywords internal
#' @noRd
print_hdMTD_summary <- function(object) {
  cat("hdMTD lag selection\n")
  if (!is.null(object$call)) {
    cat("\nCall:\n"); print(object$call)
  }
  cat("\nMethod: ", object$method, "\n", sep = "")
  cat("Order upper bound (d): ", object$d, "\n", sep = "")
  cat("Selected S set: ", fmt_vec(object$S), "\n", sep = "")

  if (!is.null(object$BIC_selected)) {
    cat("\nBIC at selected set: ", object$BIC_selected, "\n", sep = "")
  }

  lag_str <- if (length(object$lags)) fmt_vec(object$lags) else "empty"
  cat("\nRelevant lag set estimated by ", object$method, " method : ",
      lag_str, "\n", sep = "")

  if (!is.null(object$BIC_out)) {
    cat("\nBIC output (as returned by hdMTD_BIC):\n")
    print(object$BIC_out)
  }

  if (isTRUE(attr(object, "include_settings")) && !is.null(object$settings)) {
    cat("\nSettings used:\n")
    utils::str(object$settings, give.attr = FALSE, no.list = TRUE)
  }
  invisible(object)
}
