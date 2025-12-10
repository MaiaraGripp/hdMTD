#' Accessors for objects of classes \code{"MTD"}, \code{"MTDest"}, and \code{"hdMTD"}
#'
#' @description
#' Public accessors that expose model components without relying on the internal
#' list structure. These accessors are available for \code{"MTD"} (model
#' objects), \code{"MTDest"} (EM fits), and \code{"hdMTD"} (lag selection)..
#'
#' @details
#' Returned lag sets follow the package convention and are shown as negative
#' integers via \code{lags()} (elements of \eqn{\mathbb{Z}^-}). For convenience,
#' positive-index accessors are also provided:
#' \code{Lambda()} for \code{"MTD"} objects and \code{"MTDest"} fits, and
#' \code{S()} for \code{"MTDest"} and \code{"hdMTD"} objects
#' (elements of \eqn{\mathbb{N}^+}).
#'
#' The function \code{transitP()} returns the global transition matrix, i.e., a
#' convex combination of the independent distribution \code{p0} and the
#' lag-specific transition matrices \code{pj}, weighted by \code{lambdas}. For
#' \code{"MTDest"} objects, this automatically calls \code{\link{as.MTD}} to
#' first convert the fit into a proper \code{"MTD"} model.
#'
#' @note
#' Naming conventions reflect the conceptual distinction:
#' \itemize{
#'   \item \code{Lambda}: The true (known) relevant lag set in \code{"MTD"} models
#'   \item \code{S}: An estimated or candidate lag set in \code{"MTDest"} and
#'         \code{"hdMTD"} objects
#' }
#' All accessors for \code{"MTD"} work for \code{"MTDest"} also.
#'
#' @param object An object of class \code{"MTD"}, \code{"MTDest"} or \code{"hdMTD"}
#'  (as supported by each accessor).
#'
#' @return
#' \describe{
#'   \item{\code{pj(object)}}{A \code{list} of stochastic matrices (one per lag).}
#'   \item{\code{p0(object)}}{A numeric probability vector for the independent component.}
#'   \item{\code{lambdas(object)}}{A numeric vector of mixture weights that sums to 1.}
#'   \item{\code{lags(object)}}{The lag set (elements of \eqn{\mathbb{Z}^-}).}
#'   \item{\code{Lambda(object)}}{The set of relevant lags as positive integers
#'    (elements of \eqn{\mathbb{N}^+}).}
#'   \item{\code{S(object)}}{For \code{"MTDest"} and \code{"hdMTD"}, the set of
#'   candidate/estimated lags as positive integers (elements of \eqn{\mathbb{N}^+}).}
#'   \item{\code{states(object)}}{The state space.}
#'   \item{\code{transitP(object)}}{The global transition matrix \eqn{P}.}
#' }
#'
#' @seealso \code{\link{MTDmodel}}, \code{\link{MTDest}}, \code{\link{hdMTD}}, \code{\link{as.MTD}}
#'
#' @examples
#' \dontrun{
#' ## For generating an MTD model
#' set.seed(1)
#' m <- MTDmodel(Lambda = c(1, 3), A = c(0, 1))
#' pj(m); p0(m); lambdas(m); lags(m); Lambda(m); states(m)
#' transitP(m)
#' ## For an EM fit (using coef(m) as init for simplicity):
#' X <- perfectSample(m, N = 800)
#' fit <- MTDest(X, S = c(1, 3), init = coef(m))
#' pj(fit); p0(fit); lambdas(fit); lags(fit); S(fit); states(fit)
#' transitP(fit)
#' ## For lag selection:
#' S_hat <- hdMTD(X, d = 5, method = "FS", l = 2)
#' S(S_hat); lags(S_hat)
#' }
#'
#' @name MTD-accessors
NULL

# ===== Generics (exported) =====

#' @rdname MTD-accessors
#' @export
pj <- function(object) {
  UseMethod("pj")
}

#' @rdname MTD-accessors
#' @export
p0 <- function(object) {
  UseMethod("p0")
}

#' @rdname MTD-accessors
#' @export
lambdas <- function(object) {
  UseMethod("lambdas")
}

#' @rdname MTD-accessors
#' @export
lags <- function(object) {
  UseMethod("lags")
}

#' @rdname MTD-accessors
#' @export
Lambda <- function(object) {
  UseMethod("Lambda")
}

#' @rdname MTD-accessors
#' @export
S <- function(object) {
  UseMethod("S")
}

#' @rdname MTD-accessors
#' @export
states <- function(object) {
  UseMethod("states")
}

#' @rdname MTD-accessors
#' @export
transitP <- function(object) {
  UseMethod("transitP")
}

# ===== MTD obj methods =====

#' @exportS3Method pj MTD
pj.MTD <- function(object) {
  object$pj
}

#' @exportS3Method p0 MTD
p0.MTD <- function(object) {
  object$p0
}

#' @exportS3Method lambdas MTD
lambdas.MTD <- function(object) {
  object$lambdas
}

#' @exportS3Method Lambda MTD
Lambda.MTD <- function(object) {
  if(inherits(object, "MTDest")) {
    return(object$S)
  } else {
    return(object$Lambda)
  }
}

#' @exportS3Method lags MTD
lags.MTD <- function(object) {
  -Lambda(object) # lags are in Z^-
}

#' @exportS3Method states MTD
states.MTD <- function(object) {
  object$A
}

#' @exportS3Method transitP MTD
transitP.MTD <- function(object) {
  if(inherits(object, "MTDest")) {object <- as.MTD(object)}
  object$P
}

# ===== MTDest obj methods =====

#' @exportS3Method S MTDest
S.MTDest <- function(object) {
  object$S
}

# ===== hdMTD obj methods =====

#' @exportS3Method S hdMTD
S.hdMTD <- function(object) {
  as.integer(unclass(object))
}

#' @exportS3Method lags hdMTD
lags.hdMTD <- function(object) {
  -S(object)
}
