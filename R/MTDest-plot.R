#' Plot method for MTDest objects
#'
#' Produces bar plots for an \code{MTDest} object. By default, it shows
#' two plots in sequence: (i) oscillations by relevant lag, (ii) mixture weights
#' \eqn{\lambda_j} (including \code{lam0} if \code{> 0}), and - if available -
#' (iii) log-likelihood variation per update. When \code{type} is
#' specified, only the requested plot is drawn.
#'
#' @details
#' Produces the same bar plots as \code{\link{plot.MTD}} for a fitted object of
#' class \code{"MTDest"}. Internally, the object is converted to an \code{"MTD"}
#' via \code{as.MTD(x)} and then passed to \code{plot.MTD()}.
#' If EM iteration diagnostics are available (i.e., the object was fitted with
#' \code{iter = TRUE} and \code{length(deltaLogLik) > 0}), a third plot showing
#' the log-likelihood variation per update is displayed automatically when
#' \code{type} is missing. You can also request it explicitly
#' with \code{type = "convergence"}.
#'
#' @param x An object of class \code{"MTDest"}.
#' @param type Character string indicating what to plot: \code{"oscillation"},
#'   \code{"lambdas"}, or \code{"convergence"}. If \code{type} is missing,
#'   oscillations and lambdas are shown sequentially (press Enter to proceed),
#'   and—if available—the convergence plot is shown last.
#' @param main Optional main title. If missing, an informative default is used.
#'   (For \code{type = "convergence"}, the default title is
#'   \code{"EM convergence: logLik variation per update"}).
#' @param ylim Optional y-axis limits. If missing, determined from the data.
#'   (For \code{type = "convergence"}, limits are chosen from the diagnostics.)
#' @param col Color for bars (for bar plots). Ignored for \code{type = "convergence"}.
#'   Defaults to \code{"gray70"}.
#' @param border Bar border color (for bar plots). Defaults to \code{NA}.
#' @param ... Further graphical parameters passed to the underlying plotting functions.
#'
#' @return
#' If \code{type} is provided, invisibly returns the numeric vector that was plotted
#' (\code{oscillation}, \code{lambdas} or \code{deltaLogLik}). When \code{type}
#' is missing invisibly returns a list with components \code{oscillation},
#' \code{lambdas}, and - if available - \code{deltaLogLik}.
#'
#' @seealso \code{\link{plot.MTD}}, \code{\link{as.MTD}}, \code{\link{MTDest}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' M <- MTDmodel(Lambda = c(1, 3), A = c(0, 1), lam0 = 0.05)
#' X <- perfectSample(M, N = 300)
#' fit <- MTDest(X, S = c(1, 3), init = coef(M), iter = TRUE)
#'
#' plot(fit)                        # oscillations, lambdas, then (if available) convergence
#' plot(fit, type = "convergence")  # convergence panel only
#' }
#'
#' @importFrom graphics par barplot plot lines
#' @exportS3Method plot MTDest
plot.MTDest <- function(x, type, main, ylim, col = "gray70", border = NA, ...) {
  # Base plots (oscillation / lambdas) come from plot.MTD on as.MTD(x)
  m <- as.MTD(x)

  # Convergence data (may be absent or length 0)
  conv <- x$deltaLogLik
  has_conv <- !is.null(conv) && length(conv) > 0L

  # If no explicit type: show seq oscillation -> lambdas -> (optionally) convergence
  if (missing(type)) {
    old_ask <- par("ask")
    par(ask = TRUE)
    on.exit(par(ask = old_ask))

    out1 <- plot(m, type = "oscillation",
                 main = if (missing(main)) "Oscillations by lag" else main,
                 ylim = ylim, col = col, border = border, ...)

    out2 <- plot(m, type = "lambdas",
                 main = if (missing(main)) "MTD weights by lag" else main,
                 ylim = ylim, col = col, border = border, ...)

    if (has_conv) {
      # Convergence panel
      y <- conv
      # Choose a safe ylim
      ymax <- max(y, 0)
      pad  <- max(0.05, 0.08 * ymax)
      ylim3 <- if (missing(ylim)) c(0, ymax + pad) else ylim

      plot(seq_along(y), y, type = "b",
           xlab = "Iteration", ylab = "logLik variation",
           main = if (missing(main)) "EM convergence: logLik variation per update" else main,
           ylim = ylim3, ...)
      invisible(list(oscillation = out1, lambdas = out2, deltaLogLik = y))
    } else {
      invisible(list(oscillation = out1, lambdas = out2))
    }

  } else {
    type <- match.arg(type, c("oscillation", "lambdas", "convergence"))

    if (type %in% c("oscillation", "lambdas")) {
      # Delegate to plot.MTD
      invisible(plot(m, type = type, main = main, ylim = ylim, col = col, border = border, ...))

    } else { # type == "convergence"
      if (!has_conv) {
        stop("No EM diagnostics available: 'deltaLogLik' is missing or empty. This may happen if the algorithm converged in the first iteration or made no updates. If you did not set iter = TRUE, please refit with iter = TRUE to record diagnostics.")
      }
      y <- conv
      ymax <- max(y, 0)
      pad  <- max(0.05, 0.08 * ymax)
      if (missing(ylim)) ylim <- c(0, ymax + pad)

      plot(seq_along(y), y, type = "b",
           xlab = "Iteration", ylab = "logLik variation",
           main = if (missing(main)) "EM convergence: logLik variation per update" else main,
           ylim = ylim, ...)
      invisible(y)
    }
  }
}
