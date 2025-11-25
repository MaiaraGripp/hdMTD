#' Plot method for MTDest objects
#'
#' Produces plots for an \code{MTDest} object. By default, it shows in sequence:
#' (i) barplot of oscillations by relevant lag, (ii) barplot of mixture weights
#' \eqn{\lambda_j} (including \code{lam0} if \code{> 0}), and (iii) graphs of
#' \code{pj} (one graph for each lag in \code{Lambda}). If EM diagnostics are
#' available, a convergence panel (log-likelihood variation per update) is shown
#' last. When \code{type} is specified, only the requested plot is drawn.
#'
#' @details
#' For \code{type = "oscillation"}, the function calls \code{oscillation(as.MTD(x))}
#' to obtain \eqn{\delta_j = \lambda_j \max_{b,c} d_{TV}(p_j(\cdot|b), p_j(\cdot|c))}
#' for each lag in \code{Lambda(as.MTD(x))}, and draws a bar plot named by the lags.
#'
#' For \code{type = "lambdas"}, it plots the mixture weights \eqn{\lambda_j} by lag.
#' If \code{lam0 > 0}, the weight for the independent component is included and
#' labeled \code{"0"}.
#'
#' For \code{type = "pj"}, the function draws the directed, weighted graph of the
#' transition matrices \code{pj} taken from \code{pj(x)}. Vertices correspond to the
#' states \code{states(x)}. A directed edge \code{a -> b} carries weight
#' \code{p_j(b | a)}. Edge widths and edge labels are proportional to the transition
#' probabilities (labels shown in the plot are rounded to two decimals).
#' By default, self-loops (\code{a -> a}) are not drawn. The self-loop probability
#' at a state \code{a} can be inferred as \deqn{1 - \sum_{b: b \ne a} p_j(b \mid a).}
#' For \code{type = "pj"}, a specific matrix can be selected via \code{pj_index}
#' (e.g., \code{pj_index = 2} plots \code{pj(x)[[2]]}). In automatic mode (when
#' \code{type} is missing), graphs for all \code{pj} matrices are shown
#' sequentially.
#'
#' If EM iteration diagnostics are available (i.e., the object was fitted with
#' \code{iter = TRUE} and \code{length(deltaLogLik) > 0}), a convergence panel
#' showing the log-likelihood variation per update is displayed automatically
#' when \code{type} is missing. You can also request it explicitly with
#' \code{type = "convergence"}.
#'
#' @param x An object of class \code{"MTDest"}.
#' @param type If \code{type} is missing, \code{oscillation}, \code{lambdas},
#' \code{pj} graphs and—if available—the \code{convergence} plots are shown
#' sequentially (press Enter to proceed). Else, \code{type} is a character string
#' indicating what to plot: \code{"oscillation"}, \code{"lambdas"}, \code{"pj"} or
#' \code{"convergence"}.
#' @param main Optional main title. When \code{type} is missing, panel-specific
#' defaults are used and \code{main} is ignored.
#' @param ylim Optional y-axis limits. When \code{type} is missing, panel-specific
#' defaults are used and \code{ylim} is ignored.
#' @param col Bar fill color (passed to \code{barplot}). Defaults to \code{"gray70"}.
#' @param border Bar border color (for bar plots). Defaults to \code{NA}.
#' @param pj_index Integer index of \code{pj(x)} to plot when
#'   \code{type = "pj"} (default \code{1}).
#' @param ... Further graphical parameters: passed to \code{barplot()} when
#'   \code{type \%in\% c("oscillation","lambdas")}, to \code{plot.igraph()} when
#'   \code{type = "pj"}, and to \code{plot()} when \code{type = "convergence"}.
#'   Ignored when \code{type} is missing.
#'
#' @return
#' If \code{type} is \code{"oscillation"} or \code{"lambdas"}, invisibly returns
#' the numeric vector that was plotted. If \code{type = "pj"}, invisibly returns
#' the selected transition matrix \code{pj(x)[[pj_index]]}. If
#' \code{type = "convergence"}, invisibly returns the numeric vector
#' \code{deltaLogLik}. If \code{type} is missing, invisibly returns the list
#' produced by \code{plot.MTD(m)} (typically with components \code{oscillation}
#' and \code{lambdas}); if \code{deltaLogLik} is available, the returned list also
#' includes \code{deltaLogLik}.
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
#' plot(fit)
#' plot(fit, type = "pj", pj_index = 2)
#' plot(fit, type = "convergence")
#' }
#'
#' @importFrom graphics par barplot plot lines
#' @exportS3Method plot MTDest
plot.MTDest <- function(x, type, main, ylim, col = "gray70", border = NA,
                        pj_index = 1, ...) {
  # Convert to MTD and delegate base plots to plot.MTD
  m <- as.MTD(x)

  # Convergence data (may be absent or length 0)
  conv <- x$deltaLogLik
  has_conv <- !is.null(conv) && length(conv) > 0L

  if (missing(type)) {
    old_ask <- par("ask")
    par(ask = TRUE)
    on.exit(par(ask = old_ask))

    out_m <- plot(m, main = main, ylim = ylim, col = col, border = border)

    if (has_conv) {
      y <- conv
      ymax <- max(y, 0)
      pad  <- max(0.05, 0.08 * ymax)
      ylim3 <- c(0, ymax + pad)

      plot(seq_along(y), y, type = "b",
           xlab = "Iteration", ylab = "logLik variation",
           main = "EM convergence: logLik variation per update",
           ylim = ylim3)
      invisible(c(out_m, list(deltaLogLik = y)))
    } else {
      invisible(out_m)
    }

  } else {
    type <- match.arg(type, c("oscillation", "lambdas", "pj", "convergence"))

    if (type %in% c("oscillation", "lambdas", "pj")) {
      # Delegate directly to plot.MTD
      res <- plot(m, type = type, main = main, ylim = ylim,
                  col = col, border = border,
                  pj_index = pj_index, ...)
      invisible(res)

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
