#' Predictive probabilities for MTD / MTDest
#'
#' @description
#' Compute one-step-ahead predictive probabilities under an MTD model or an MTDest fit.
#'
#' Conventions:
#' - Samples are read most recent first: `x[1] = X_{t-1}`, `x[2] = X_{t-2}`, etc.
#' - The global transition matrix `P` is indexed by row labels that list the past
#'   context from oldest to newest. A cell at row `"s_k...s_1"` and column `"a"`
#'   is read as `p(a | s_1...s_k)`.
#' - If both `newdata` and `context` are missing, `probs()` returns the full
#'   global transition matrix (`transitP(object)` for `MTD`; `transitP(as.MTD(object))`
#'   for `MTDest`).
#'
#' @param object An `MTD` or `MTDest` object.
#' @param context Optional vector or matrix/data.frame of contexts (rows). By default,
#'   each row follows the "most recent first" convention; set `oldLeft = TRUE` if rows
#'   are supplied oldest to newest. Must have exactly `length(Lambda(object))` columns
#'   for `MTD` or `length(S(object))` for `MTDest` (one symbol per lag).
#' @param newdata Optional vector or matrix/data.frame of samples (rows). Columns follow
#'   the "most recent first" convention. Must have at least `max(Lambda(object))` columns
#'   for `MTD` or `max(S(object))` for `MTDest`. Only one of `newdata` or `context`
#'   can be provided at a time. When there are extra columns, only the columns at the
#'   model lags are used.
#' @param oldLeft Logical. If `TRUE`, interpret rows in `newdata`/`context` as
#' oldest to newest (e.g. leftmost = `newdata[ ,1]` = oldest). If `FALSE` (default),
#' rows are most recent first.
#'
#' @return A numeric matrix of predictive probabilities with one row per input context and
#'   columns indexed by `states(object)`. Row names are the context labels (oldest to newest)
#'   formed by concatenating state symbols without a separator. If both `newdata` and
#'   `context` are missing, the full global transition matrix is returned.
#'
#' @details
#' All entries of `newdata`/`context` must belong to the model's state space `states(object)`.
#' For `MTDest`, returning the full matrix materializes `transitP(as.MTD(object))`, which can be
#' large for big state spaces or many lags.
#'
#' @examples
#' set.seed(1)
#' m <- MTDmodel(Lambda = c(1,3), A = c(0,1), lam0 = 0.1)
#'
#' # Full matrix
#' P <- probs(m)
#'
#' # Using a sample row (most recent first): newdata has >= max(Lambda) columns
#' new_ctx <- c(1, 0, 1, 0)      # X_{t-1}=1, X_{t-2}=0, X_{t-3}=1, ...
#' probs(m, newdata = new_ctx)   # one row of probabilities
#'
#' # Explicit contexts (exactly |Lambda| symbols per row)
#' probs(m, context = c(0, 1), oldLeft = FALSE)  # most recent first
#' probs(m, context = c(0, 1), oldLeft = TRUE)   # oldest to newest
#'
#' # Multiple contexts (rows)
#' ctxs <- rbind(c(1,0,1), c(0,1,1), c(1,1,0))
#' probs(m, newdata = ctxs)
#'
#' @seealso \code{\link{transitP}}, \code{\link{states}}, \code{\link{Lambda}},
#'   \code{\link{S}}, \code{\link{as.MTD}}, \code{\link{empirical_probs}}
#'
#' @name probs
#' @export
probs <- function(object, context = NULL, newdata = NULL, oldLeft = FALSE) UseMethod("probs")
#'
#' @rdname probs
#' @export
probs.MTD <- function(object, context = NULL, newdata = NULL, oldLeft = FALSE) {
  stopifnot(inherits(object, "MTD"))

  P <- transitP(object)

  if (is.null(newdata) && is.null(context)) { # Full transition matrix already available
    return(P)
  }

  if (!is.null(newdata) && !is.null(context)) stop("Only one of `newdata` or `context` can be provided at a time.")

  A        <- states(object)
  S        <- sort(as.integer(Lambda(object)))  # positive lags

  if (!is.null(newdata)) {
    # --- Coerce newdata to a matrix ---
    if (is.vector(newdata)) {
      newdata <- matrix(as.vector(newdata), nrow = 1)
    } else if (is.data.frame(newdata)) {
      newdata <- as.matrix(newdata)
    } else if (!is.matrix(newdata)) {
      stop("`newdata` must be either NULL, a vector, a matrix, or a data.frame.")
    }

    if(!all(unique(as.vector(newdata)) %in% A)) stop("Some `newdata` elements are not in the state space.")
    if (ncol(newdata)<max(S)) {
      stop("`newdata` must have at least ", max(S), " symbols from the state space.")
    }

    contextP <- rownames(P)


    if(oldLeft){# User gives sample with OLDEST on the LEFT.

      cols_needed <- rev(ncol(newdata) - S + 1)
      ctx_S_order <- newdata[, cols_needed, drop = FALSE]

    } else {# User gives sample with MOST RECENT on the LEFT.

      ctx_S_order <- newdata[, rev(S), drop = FALSE]   # older to newer to compare with P
    }

    labels <- apply(ctx_S_order, 1, paste0, collapse = "")
    rows   <- match(labels, contextP)

    out <- P[rows, , drop = FALSE]
    rownames(out) <- labels
    return(out)
  }

  if (!is.null(context)){
    # --- Coerce context to a matrix with rows = contexts ---
    if (is.vector(context)) {
      context <- matrix(as.vector(context), nrow = 1)
    } else if (is.data.frame(context)) {
      context <- as.matrix(context)
    } else if (!is.matrix(context)) {
      stop("`context` must be either NULL, a vector, a matrix, or a data.frame.")
    }

    if(!all(unique(as.vector(context)) %in% A)) stop("Some `context` elements are not in the state space.")
    if (ncol(context) != length(S)) {
      stop("`context` must have exactly ", length(S), " symbols (one per lag in Lambda).")
    }
    contextP <- rownames(P)


    if(oldLeft){# User gives context with OLDEST on the LEFT.

      ctx_oldfirst <- context

    } else {# User gives sample with MOST RECENT on the LEFT.

      ctx_oldfirst <- context[, ncol(context):1, drop = FALSE]   # older to newer to compare with P

    }
    labels <- apply(ctx_oldfirst, 1, paste0, collapse = "")
    rows   <- match(labels, contextP)
    out <- P[rows, , drop = FALSE]
    rownames(out) <- labels
    return(out)
  }
}
#'
#' @rdname probs
#' @export
probs.MTDest <- function(object, context = NULL, newdata = NULL, oldLeft = FALSE) {
  stopifnot(inherits(object, "MTDest"))

  P <- transitP(as.MTD(object))

  if (is.null(newdata) && is.null(context)) { # Full transition matrix
    return(P)
  }

  if (!is.null(newdata) && !is.null(context)) stop("Only one of `newdata` or `context` can be provided at a time.")

  A        <- states(object)
  S        <- sort(as.integer(S(object)))  # positive lags

  if (!is.null(newdata)) {
    # --- Coerce newdata to a matrix ---
    if (is.vector(newdata)) {
      newdata <- matrix(as.vector(newdata), nrow = 1)
    } else if (is.data.frame(newdata)) {
      newdata <- as.matrix(newdata)
    } else if (!is.matrix(newdata)) {
      stop("`newdata` must be either NULL, a vector, a matrix, or a data.frame.")
    }

    if(!all(unique(as.vector(newdata)) %in% A)) stop("Some `newdata` elements are not in the state space.")
    if (ncol(newdata)<max(S)) {
      stop("`newdata` must have at least ", max(S), " symbols from the state space.")
    }

    contextP <- rownames(P)


    if(oldLeft){# User gives sample with OLDEST on the LEFT.

      cols_needed <- rev(ncol(newdata) - S + 1)
      ctx_S_order <- newdata[, cols_needed, drop = FALSE]

    } else {# User gives sample with MOST RECENT on the LEFT.

      ctx_S_order <- newdata[, rev(S), drop = FALSE]   # older to newer to compare with P
    }

    labels <- apply(ctx_S_order, 1, paste0, collapse = "")
    rows   <- match(labels, contextP)

    out <- P[rows, , drop = FALSE]
    rownames(out) <- labels
    return(out)
  }

  if (!is.null(context)){
    # --- Coerce context to a matrix with rows = contexts ---
    if (is.vector(context)) {
      context <- matrix(as.vector(context), nrow = 1)
    } else if (is.data.frame(context)) {
      context <- as.matrix(context)
    } else if (!is.matrix(context)) {
      stop("`context` must be either NULL, a vector, a matrix, or a data.frame.")
    }

    if(!all(unique(as.vector(context)) %in% A)) stop("Some `context` elements are not in the state space.")
    if (ncol(context) != length(S)) {
      stop("`context` must have exactly ", length(S), " symbols (one per lag in S).")
    }
    contextP <- rownames(P)


    if(oldLeft){# User gives context with OLDEST on the LEFT.

      ctx_oldfirst <- context

    } else {# User gives sample with MOST RECENT on the LEFT.

      ctx_oldfirst <- context[, ncol(context):1, drop = FALSE]   # older to newer to compare with P

    }
    labels <- apply(ctx_oldfirst, 1, paste0, collapse = "")
    rows   <- match(labels, contextP)
    out <- P[rows, , drop = FALSE]
    rownames(out) <- labels
    return(out)
  }

}
