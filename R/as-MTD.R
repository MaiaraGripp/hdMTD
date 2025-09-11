#' Coerce an EM fit to an MTD model
#'
#' @description
#' Convenience coercion to rebuild an object of class \code{"MTD"} from an
#' \code{"MTDest"} fit (or its summary). This simply feeds the estimated
#' parameters back into \code{\link{MTDmodel}}.
#'
#' @param x An object of class \code{"MTDest"} or \code{"summary.MTDest"}.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return An object of class \code{"MTD"} as returned by \code{\link{MTDmodel}}.
#'
#' @seealso \code{\link{MTDest}}, \code{\link{MTDmodel}}
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   MTD <- MTDmodel(Lambda = c(1, 3), A = c(0, 1), lam0 = 0.05) # generates MTD model
#'   X <- perfectSample(MTD, N = 400) # generates MTD sample
#'   init <- list(
#'     p0 = c(0.4, 0.6),
#'     lambdas = c(0.05, 0.45, 0.5),
#'     pj = list(
#'       matrix(c(0.2, 0.8, 0.45, 0.55), byrow = TRUE, ncol = 2),
#'       matrix(c(0.25, 0.75, 0.3, 0.7),  byrow = TRUE, ncol = 2)
#'     )
#'   )
#'   fit <- MTDest(X, S = c(1, 3), init = init) # estimates parameters from sample
#'   m <- as.MTD(fit) # generates an MTD model from estimated parameters
#'   str(m, max.level = 1)
#' }
#'
#' @export
as.MTD <- function(x, ...) UseMethod("as.MTD")

#' @exportS3Method as.MTD MTDest
as.MTD.MTDest <- function(x, ...) {
  stopifnot(inherits(x, "MTDest"))
  MTDmodel(
    Lambda = x$S,
    A      = x$A,
    lam0   = x$lambdas[1],
    lamj   = x$lambdas[-1],
    pj     = x$pj,
    p0     = x$p0,
    single_matrix = FALSE,
    indep_part    = any(x$p0 != 0)
  )
}

#' @exportS3Method as.MTD summary.MTDest
as.MTD.summary.MTDest <- function(x, ...) {
  stopifnot(inherits(x, "summary.MTDest"))
  MTDmodel(
    Lambda = x$S,
    A      = x$A,
    lam0   = x$lambdas[1],
    lamj   = x$lambdas[-1],
    pj     = x$pj,
    p0     = if (is.null(x$p0)) rep(0, length(x$A)) else x$p0,
    single_matrix = FALSE,
    indep_part    = !is.null(x$p0) && any(x$p0 != 0)
  )
}
