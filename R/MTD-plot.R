#' Plot method for MTD objects
#'
#' Produces plots for an \code{MTD} object. By default, it shows the following
#' sequence of plots: (i) barplot of oscillations by relevant lag,
#' (ii) barplot of mixture weights \eqn{\lambda_j} (including \code{lam0} if \code{> 0})
#' and (iii) graphs of \code{pj}, one graph for each lag in \code{Lambda}. When
#' \code{type} is specified, only the requested plot is drawn.
#'
#' @param x An object of class \code{"MTD"}.
#' @param type  If \code{type} is missing, all plots are shown sequentially
#'  (press Enter to proceed). Else, \code{type} is a character string indicating
#'  what to plot: \code{"oscillation"}, \code{"lambdas"} or \code{"pj"}.
#' @param main Optional main title. When \code{type} is missing, panel-specific
#' defaults are used and \code{main} is ignored.
#' @param ylim Optional y-axis limits. When \code{type} is missing, panel-specific
#'   defaults are used and \code{ylim} is ignored.
#' @param col Bar fill color (passed to \code{barplot}). Defaults to \code{"gray70"}.
#' @param border Bar border color (passed to \code{barplot}). Defaults to \code{NA}.
#' @param pj_index Integer index of \code{pj(x)} to plot when \code{type = "pj"}
#'  (default \code{1}).
#' @param ... Further graphical parameters: passed to \code{barplot()} when
#'   \code{type \%in\% c("oscillation","lambdas")} and to \code{plot.igraph()}
#'   when \code{type = "pj"}. Ignored when \code{type} is missing.
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
#' For \code{type = "pj"}, the function draws the directed, weighted graph of a
#' transition matrix \code{pj} taken from \code{pj(x)}. Vertices correspond to the
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
#' @return
#' If \code{type} is \code{"oscillation"} or \code{"lambdas"}, invisibly returns the numeric
#' vector that was plotted. If \code{type = "pj"}, invisibly returns the selected transition
#' matrix. If \code{type} is missing, invisibly returns a list with components
#' \code{oscillation} and \code{lambdas}.
#'
#' @seealso \code{\link{oscillation}}, \code{\link{lambdas}}, \code{\link{Lambda}}
#' @importFrom graphics par barplot
#' @importFrom igraph graph_from_adjacency_matrix layout_nicely
#' @examples
#' \dontrun{
#' m <- MTDmodel(Lambda = c(1, 3), A = c(0, 1))
#'
#' ## Automatic mode (press Enter between plots)
#' plot(m)
#'
#' ## Single plot:
#' plot(m, type = "oscillation")
#' plot(m, type = "lambdas")
#' plot(m, type = "pj", pj_index = 2)
#' }
#'
#' @exportS3Method plot MTD
plot.MTD <- function(x, type, main, ylim, col = "gray70", border = NA,
                     pj_index = 1,...) {
  checkMTD(x)

  if (missing(type)) {
    old_ask <- par("ask")
    par(ask = TRUE)
    on.exit(par(ask = old_ask))

    ## 1) Oscillations
    y1 <- oscillation(x)
    main1 <- "Oscillations by lag"
    ymax1 <- max(y1, 0)
    pad1 <- max(0.05, 0.08 * ymax1)
    ylim1 <- c(0, ymax1 + pad1)

    barplot(y1, names.arg = names(y1), ylim = ylim1,
            main = main1, xlab = "Relevant lags",
            col = col, border = border)

    ## 2) Lambdas
    lj <- lambdas(x)
    lam0 <- as.numeric(lj[1])
    lag_names <- paste0(lags(x))
    if (lam0 > 0) {lag_names <- c("0", lag_names)}
    y2 <- if (lam0 > 0) lj else lj[-1]
    names(y2) <- lag_names
    main2 <- "MTD weights by lag"
    ymax2 <- max(y2, 0)
    pad2 <- max(0.05, 0.08 * ymax2)
    ylim2 <- c(0, ymax2 + pad2)

    barplot(y2, names.arg = names(y2), ylim = ylim2,
            main = main2, xlab = "Relevant lags",
            col = col, border = border)

    ## 3) Graphs of pj
    P <- pj(x) # list
    pjNames <- names(P)
    for (k in seq_along(P)) {
      mat <- P[[k]]
      vNames <- colnames(mat) #states(x)
      if (is.null(vNames) || !identical(vNames, as.character(states(x)))) {
        vNames <- as.character(states(x))
      }
      g <- igraph::graph_from_adjacency_matrix(
        mat, mode = "directed",
        weighted = TRUE, diag = FALSE
      )
      if (!is.null(vNames)) igraph::V(g)$name <- vNames

      main3 <- if (!is.null(pjNames) && nzchar(pjNames[k])) {
        paste0("Transition graph (pj): ", pjNames[k])
      } else {
        paste0("Transition graph (pj): p-", Lambda(x)[k])
      }

      plot(g,
           main = main3,
           layout = igraph::layout_nicely(g),
           edge.width = 5 * igraph::E(g)$weight,
           edge.label = round(igraph::E(g)$weight, 2),
           edge.color = "gray",
           edge.curved = 0.2,
           edge.arrow.size = 0.45,
           vertex.size = 30,
           vertex.label.color = "white",
           vertex.label.font = 2,
           vertex.color = "darkblue"
           )

    }

    invisible(list(oscillation = y1, lambdas = y2))

  } else {

    type <- match.arg(type, c("oscillation", "lambdas", "pj"))

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

    } else if (type == "lambdas") {
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

    } else { #type = pj
      P <- pj(x)
      if (length(P) == 0) stop("pj(x) returned an empty list.")
      if (pj_index < 1 || pj_index > length(P) || pj_index%%1 != 0) {
        stop(sprintf("pj_index must be an integer in 1..%d", length(P)))
      }
      mat <- P[[pj_index]]
      vNames <- colnames(mat)
      if (is.null(vNames) || !identical(vNames, as.character(states(x)))) {
        vNames <- as.character(states(x))
      }

      g <- igraph::graph_from_adjacency_matrix(
        mat, mode = "directed", weighted = TRUE, diag = FALSE
      )
      if (!is.null(vNames)) igraph::V(g)$name <- vNames

      if (missing(main)) {
        pjNames <- names(P)
        main <- if (!is.null(pjNames) && nzchar(pjNames[pj_index])) {
          paste0("Transition graph: ", pjNames[pj_index])
        } else {
          paste0("Transition graph: p-", Lambda(x)[pj_index])
        }
      }

      ig_defaults <- list(
        main = main,
        layout = igraph::layout_nicely(g),
        edge.width = 5 * igraph::E(g)$weight,
        edge.label = round(igraph::E(g)$weight, 2),
        edge.color = "gray",
        edge.curved = 0.2,
        edge.arrow.size = 0.45,
        vertex.size = 30,
        vertex.label.color = "white",
        vertex.label.font = 2,
        vertex.color = "darkblue"
      )

      dots <- list(...)

      if (length(dots)) {
        ig_defaults[names(dots)] <- NULL
      }

      do.call(plot, c(list(g), ig_defaults, dots))

      invisible(mat)
    }
  }
}
