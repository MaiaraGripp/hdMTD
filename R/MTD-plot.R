#' Plot method for MTD objects
#'
#' Produces bar plots for an \code{MTD} object. By default, it shows
#' two plots in sequence: (i) oscillations by relevant lag, and (ii) mixture weights
#' \eqn{\lambda_j} (including \code{lam0} if \code{> 0}). When \code{type} is
#' specified, only the requested plot is drawn.
#'
#' @param x An object of class \code{"MTD"}.
#' @param type Character string indicating what to plot: \code{"oscillation"}
#'   for the oscillations by lag, or \code{"lambdas"} for the weights by lag.
#'   If \code{type} is missing, both plots are shown sequentially (press Enter to proceed).
#' @param main Optional main title. If missing, an informative default is used.
#' @param ylim Optional y-axis limits. If missing, determined from the data.
#' @param col Bar fill color (passed to \code{barplot}). Defaults to \code{"gray70"}.
#' @param border Bar border color (passed to \code{barplot}). Defaults to \code{NA}.
#' @param ... Further graphical parameters passed to \code{barplot()}.
#'
#' @details
#' For \code{type = "oscillation"}, the function calls \code{oscillation(x)}
#' to obtain \eqn{\delta_j = \lambda_j \max_{b,c} d_{TV}(p_j(\cdot|b), p_j(\cdot|c))}
#' for each lag in \code{Lambda(x)}, and draws a bar plot named by the lags.
#'
#' For \code{type = "lambdas"}, it plots the mixture weights \eqn{\lambda_j} by lag.
#' If \code{lam0 > 0}, the weight for the independent component is included and
#' labeled \code{"0"}.
#'
#'
#' @return
#' If \code{type} is provided, invisibly returns the numeric vector that was plotted
#' (oscillation or lambdas). If \code{type} is missing, invisibly returns a list with
#' components \code{oscillation} and \code{lambdas}.
#'
#' @seealso \code{\link{oscillation}}, \code{\link{lambdas}}, \code{\link{Lambda}}
#' @importFrom graphics par barplot
#' @examples
#' \dontrun{
#' m <- MTDmodel(Lambda = c(1, 3), A = c(0, 1))
#'
#' ## Automatic mode: shows oscillations then lambdas (press Enter between plots)
#' plot(m)
#'
#' ## Single plot:
#' plot(m, type = "oscillation")
#' plot(m, type = "lambdas")
#' }
#'
#' @exportS3Method plot MTD
plot.MTD <- function(x, type = c("oscillation", "lambdas"),
                     main, ylim, col = "gray70", border = NA, ...) {
  checkMTD(x)

  if (missing(type)) {
    oldpar <- par(ask = TRUE)
    on.exit(par(oldpar))

    ## 1) Oscillations
    y1 <- oscillation(x)
    main1 <- if (missing(main)) "Oscillations by lag" else main
    if (missing(ylim)) {
      ymax1 <- max(y1, 0)
      pad1 <- max(0.05, 0.08 * ymax1)
      ylim1 <- c(0, ymax1 + pad1)
    } else ylim1 <- ylim
    barplot(y1, names.arg = names(y1), ylim = ylim1,
            main = main1, xlab = "Relevant lags",
            col = col, border = border, ...)

    ## 2) Lambdas
    lj <- lambdas(x)
    lam0 <- as.numeric(lj[1])
    lag_names <- paste0(lags(x))
    if (lam0 > 0) {lag_names <- c("0", lag_names)}
    y2 <- if (lam0 > 0) lj else lj[-1]; names(y2) <- lag_names
    main2 <- if (missing(main)) "MTD weights by lag" else main
    ymax2 <- max(y2, 0)
    pad2 <- max(0.05, 0.08 * ymax2)
    ylim2 <- c(0, ymax2 + pad2)
    barplot(y2, names.arg = names(y2), ylim = ylim2,
            main = main2, xlab = "Relevant lags",
            col = col, border = border, ...)

    invisible(list(oscillation = y1, lambdas = y2))

  } else {
    type <- match.arg(type)

    if (type == "oscillation") {
      y <- oscillation(x)
      if (missing(main)) main <- "Oscillations by lag"
      if (missing(ylim)) {
        ymax <- max(y, 0)
        pad <- max(0.05, 0.08 * ymax)
        ylim <- c(0, ymax + pad)
      }
      barplot(y, names.arg = names(y), ylim = ylim,
              main = main, xlab = "Relevant lags",
              col = col, border = border, ...)
      invisible(y)

    } else { # type == "lambdas"
      lj <- lambdas(x)
      lam0 <- as.numeric(lj[1])
      lag_names <- paste0(lags(x))
      if (lam0 > 0) {lag_names <- c("0", lag_names)}
      y <- if (lam0 > 0) lj else lj[-1]
      names(y) <- lag_names
      if (missing(main)) main <- "MTD weights by lag"
      if (missing(ylim)) {
        ymax <- max(y, 0)
        pad <- max(0.05, 0.08 * ymax)
        ylim <- c(0, ymax + pad)
      }
      barplot(y, names.arg = names(y), ylim = ylim,
              main = main, xlab = "Relevant lags",
              col = col, border = border, ...)
      invisible(y)
    }
  }
}
