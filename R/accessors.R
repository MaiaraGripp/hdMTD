#' Accessors for objects of classes \code{"MTD"}, \code{"MTDest"}, and \code{"hdMTD"}
#'
#' @description
#' Public accessors that expose model components without relying on the internal
#' list structure. These accessors are available for \code{"MTD"} (model
#' objects), \code{"MTDest"} (EM fits), and \code{"hdMTD"} (lag selection),
#' except \code{transitP()} which only applies to \code{"MTD"}.
#'
#' @details
#' Returned lag sets follow the package convention and are shown as negative
#' integers via \code{lags()} (elements of \eqn{\mathbb{Z}^-}). For convenience,
#' positive-index accessors are also provided:
#' \code{Lambda()} for \code{"MTD"} objects and \code{S()} for \code{"MTDest"} and
#' \code{"hdMTD"} objects (elements of \eqn{\mathbb{N}^+}). Internally, lags may
#' be stored as positive integers in \code{Lambda} or \code{S}.
#'
#' For computing the global transition matrix of an EM fit the user can first
#' coerce the MTDest object into an MTD using \link{as.MTD} and then access the
#' matrix with \code{transitP(as.MTD(object))}.
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
#'   \item{\code{Lambda(object)}}{For \code{"MTD"}, the lag set as positive integers
#'    (elements of \eqn{\mathbb{N}^+}).}
#'   \item{\code{S(object)}}{For \code{"MTDest"} and \code{"hdMTD"}, the lag set as
#'    positive integers (elements of \eqn{\mathbb{N}^+}).}
#'   \item{\code{states(object)}}{The state space.}
#'   \item{\code{transitP(object)}}{For \code{"MTD"} objects only, the global
#'   transition matrix \eqn{P}. Not available for \code{"MTDest"}.}
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
#' transitP(as.MTD(fit))
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

#' @exportS3Method lags MTD
lags.MTD <- function(object) {
  -object$Lambda  # lags are in Z^-
}

#' @exportS3Method Lambda MTD
Lambda.MTD <- function(object) {
  object$Lambda
}

#' @exportS3Method states MTD
states.MTD <- function(object) {
  object$A
}

#' @exportS3Method transitP MTD
transitP.MTD <- function(object) {
  object$P
}

# ===== MTDest obj methods =====

#' @exportS3Method pj MTDest
pj.MTDest <- function(object) {
  object$pj
}

#' @exportS3Method p0 MTDest
p0.MTDest <- function(object) {
  object$p0
}

#' @exportS3Method lambdas MTDest
lambdas.MTDest <- function(object) {
  object$lambdas
}

#' @exportS3Method lags MTDest
lags.MTDest <- function(object) {
  -object$S  # lags are in Z^-
}

#' @exportS3Method S MTDest
S.MTDest <- function(object) {
  object$S
}

#' @exportS3Method states MTDest
states.MTDest <- function(object) {
  object$A
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
