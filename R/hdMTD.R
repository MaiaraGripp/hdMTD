#' Inference in MTD models
#'
#' This function estimates the relevant lag set in a Mixture Transition Distribution (MTD) model
#' using one of the available methods. By default, it applies the Forward Stepwise ("FS") method,
#' which is particularly useful in high-dimensional settings.
#'  The available methods are:
#' - "FS" (Forward Stepwise): selects the lags by a criterion that depends on their oscillations.
#' - "CUT": a method that selects the relevant lag set based on a predefined threshold.
#' - "FSC" (Forward Stepwise and Cut): applies the "FS" method followed by the "CUT" method.
#' - "BIC": selects the lag set using the Bayesian Information Criterion.
#'
#' The function dynamically calls the corresponding method function (e.g., [hdMTD_FSC()] for
#' \code{method = "FSC"}). Additional parameters specific to each method can be provided via `...`,
#' and default values are used for unspecified parameters.
#'
#' @details
#' This function serves as a wrapper for the method-specific functions:
#' - [hdMTD_FS()], for \code{method = "FS"}
#' - [hdMTD_FSC()], for \code{method = "FSC"}
#' - [hdMTD_CUT()], for \code{method = "CUT"}
#' - [hdMTD_BIC()], for \code{method = "BIC"}
#'
#' Any additional parameters (`...`) must match those accepted by the corresponding method function.
#' If a parameter value is not explicitly provided, a default value is used.
#' The main default parameters are:
#' - \code{S = seq_len(d)}: Used in "BIC" or "CUT" methods.
#' - \code{alpha = 0.05}, \code{mu = 1}. Used in "CUT" or "FSC" methods.
#' - \code{xi = 0.5}.  Used in "CUT", "FSC" or "BIC" methods.
#' - \code{minl = 1}, \code{maxl = length(S)}, \code{byl = FALSE}. Used in "BIC" method.
#' All default values are specified in the documentation of the method-specific functions.
#'
#' @param X A vector or single-column data frame containing a chain sample.
#' @param d A positive integer representing an upper bound for the chain order.
#' @param method  A character string indicating the method for estimating the relevant lag set.
#' The available methods are: "FS" (default), "FSC", "CUT", and "BIC". See the *Details* section
#' and the documentation of the corresponding method functions for more information.
#' @param ... Additional arguments for the selected method. If not specified, default values
#' will be used (see *Details* ).
#'
#' @return A vector containing the estimated relevant lag set.
#' @export
#'
#' @examples
#' # Simulate a chain from an MTD model
#' set.seed(1)
#' M <- MTDmodel(Lambda = c(1, 4), A = c(1, 3), lam0 = 0.05)
#' X <- perfectSample(M, N = 400)

#' # Fit using Forward Stepwise (FS)
#' hdMTD(X = X, d = 5, method = "FS", l = 2)

#' # Fit using Bayesian Information Criterion (BIC)
#' hdMTD(X = X, d = 5, method = "BIC", xi = 1, minl = 3, maxl = 3)
#'
hdMTD <- function(X, d, method = "FS", ...){

  method <- toupper(method)
  if(!method %in% c("FSC", "FS", "CUT", "BIC")){
    stop("The chosen method is unknown")
  }

  fmtd <-  match.fun(paste0("hdMTD_", method))
  fmtd_params <- names(formals(fmtd)) # names of parameters need for method

  params <- list(...)
  nms <- names(params)

  # Capture extra args
  if (length(params)) {
    if (is.null(nms) || any(is.na(nms)) || any(!nzchar(nms))) {
      stop("All additional arguments in '...' must be named (non-empty, non-NA).")
    }
    if (any(duplicated(nms))) {
      dup <- unique(nms[duplicated(nms)])
      stop("Duplicated argument name(s) in '...': ", paste(dup, collapse = ", "))
    }
    unknown <- setdiff(nms, fmtd_params)
    if (length(unknown)) {
      stop("Unknown argument(s) for method ", method, ": ",
           paste(unknown, collapse = ", "),
           ". See ?hdMTD_", method, " for accepted parameters.")
    }
  }

  # List of default parameters
  dparams <- list(S = seq_len(d), l = NULL, alpha = 0.05, mu = 1, xi = 0.5,
                  minl = 1, maxl = d, A = NULL, byl = FALSE, BICvalue = FALSE,
                  single_matrix = FALSE, indep_part = TRUE, zeta = d,
                  elbowTest = FALSE, warning = FALSE)

  # overlay user-provided values
  if(length(params) > 0){
      dparams[nms] <- params # Replace default parameters with informed ones

      if ( method == "BIC"){
          if(!("maxl" %in% nms)){dparams$maxl <- length(dparams$S)} # Default maxl to length(S)
          if(!("zeta" %in% nms)){dparams$zeta <- dparams$maxl} # Default zeta to maxl
      }
  }

  # Check l parameter for FS and FSc
  if (method == "FS") {
    if (isTRUE(dparams$elbowTest) && is.null(dparams$l)) {
      dparams$l <- d
    }
    # Else, l is mandatory for FS
    if (!isTRUE(dparams$elbowTest) && is.null(dparams$l)) {
      stop("For method 'FS', please supply 'l' (number of selected lags).")
    }
  }
  # l is mandatory for FSC
  if (method == "FSC") {
    if (is.null(dparams$l)) {
      stop("For method 'FSC', please supply 'l' (number of selected lags).")
    }
  }

  fmtd(X = X, d = d, S = dparams$S, l = dparams$l, alpha = dparams$alpha,
       mu = dparams$mu, xi = dparams$xi, minl = dparams$minl, maxl = dparams$maxl,
       A = dparams$A, byl = dparams$byl, BICvalue = dparams$BICvalue,
       single_matrix = dparams$single_matrix, indep_part = dparams$indep_part,
       zeta = dparams$zeta, elbowTest = dparams$elbowTest, warning = dparams$warning)
}
