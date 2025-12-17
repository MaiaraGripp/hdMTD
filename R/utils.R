# utils.R - Internal functions for the package
#
# This file contains auxiliary functions used within the package.
# These functions are not exported and are used internally to run some
# simple algorithm or calculation.

# 1 - groupTab()
# 2 - n_parameters()
# 3 - fmt_vec()
# 4 - compute_transitP()

# At the end of each auxiliary function below there is
# a note naming the functions that use it.

###########################################################
###########################################################
###########################################################

# groupTab: groups and summarizes a frequency table by a given set of lags.
#
# This function groups a tibble `freqTab` by a subset of lags `S` and possibly
# an additional lag `j`. It then sums the frequency counts (`Nxa_Sj`) for
# each unique combination of these lags.
#
# Arguments:
# - S: A numeric vector of past lags.
# - j: A single lag (integer).
# - freqTab: A tibble containing frequency counts of sequences in a column named `Nxa_Sj`.
# - lenX: Sample size (integer).
# - d: Maximum lag order (integer).
#
# Returns:
# - A tibble grouped by lags in `Sj` with summed absolute frequencies.
# - If `S` and `j` are NULL or empty, returns a 1-row matrix [0, lenX - d].

groupTab <- function(S, j, freqTab, lenX, d){

    Sj <- sort(c(S, j), decreasing = TRUE) # The set of lags

    if (length(Sj) > 0) {
        # Summarizes freqTab by lags in Sj
        groupTab <- freqTab %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(paste0("x", Sj)))) %>%
            dplyr::summarise(Nx_Sj = sum(Nxa_Sj), .groups="drop")

        return(groupTab)

    } else {
      # If no lags are provided, returns a default dataframe.
      return( data.frame(x = 0, Nx_Sj = lenX - d) )
    }
}
# groupTab is used in: oscillation.R, hdMTD_FS.R.


###########################################################
###########################################################
###########################################################

# n_parameters: Computes the number of parameters in an MTD model.
#
# This function calculates the number of parameters in an MTD model based on the
# relevant lag set `Lambda`, the state space `A`, and model constraints.
#
# Arguments:
# - Lambda: A numeric vector of positive integers representing the relevant lag set.
#   The elements are sorted in increasing order, where the smallest number
#   represents the most recent past time, and the largest represents the earliest past time.
# - A: A vector of positive integers representing the state space.
# - single_matrix: Logical. If TRUE, the MTD model assumes a single stochastic
#   matrix `p_j` shared across all lags in `Lambda`, reducing the number of parameters.
# - indep_part: Logical. If FALSE, the independent distribution is not included,
#   meaning `lambda_0 = 0`, which reduces the number of parameters.
# - zeta: A positive integer indicating the number of distinct matrices `p_j`
#   in the MTD model. The default value is `length(Lambda)`, meaning each lag has
#   its own matrix. If `zeta = 1`, all matrices `p_j` are identical; if `zeta = 2`,
#   there are two distinct types, and so on.
#
# Returns:
# - An integer representing the total number of parameters in the MTD model.
#
# Notes:
# - If `single_matrix = TRUE`, `zeta` is automatically set to 1.

n_parameters <- function(Lambda, A, single_matrix = FALSE, indep_part = TRUE, zeta = length(Lambda))
{
    lenA <- length(A)
    lenL <- length(Lambda)
    if (single_matrix) {zeta <- 1}

    n_parameters <- lenL - 1 # Number of free weight parameters if lam0 = 0
    if (indep_part) {
        n_parameters <- n_parameters + lenA  # adds 1 ( for lam0) plus lenA - 1 ( for p0)
    }
    n_parameters <- n_parameters + lenA * (lenA - 1) * zeta
    # lenA*(lenA-1) is the number of free parameters in each matrix pj. zeta is the number of distinct matrices pj
    n_parameters
} # n_parameters is used in hdMTD_BIC.R, MTDest-methods.R and MTDmodel-methods.R.


###########################################################
###########################################################
###########################################################

# fmt_vec: Pretty-prints a vector into a compact one-line string
#
# This function formats an arbitrary vector into a concise, human-readable
# string. It joins elements with commas and, when the vector is longer than
# `max_items`, it truncates the output and appends a suffix indicating the
# total length (e.g., ", ... (57 total)").
#
# Arguments:
# - x: A vector of any atomic mode (numeric, integer, character, logical, etc.).
# - max_items: Integer (default = 10). Maximum number of elements to show
#   before truncating with the "..., (n total)" suffix.
# - digits: Integer or NULL (default = NULL). If non-NULL and `x` is numeric,
#   values are formatted with `format(x, digits = digits)` prior to joining.
# - empty: Character (default = "∅"). String returned when `x` has length 0 or is NULL.
#
# Returns:
# - A single character string with the compact representation of `x`.

fmt_vec <- function(x, max_items = 10, digits = NULL, empty = "empty") {
  if (is.null(x) || length(x) == 0) return(empty)
  x <- as.vector(x)
  if (!is.null(digits) && is.numeric(x)) x <- format(x, digits = digits)
  n <- length(x)
  if (n <= max_items) {
    paste(x, collapse = ", ")
  } else {
    paste0(paste(x[seq_len(max_items)], collapse = ", "),
           ", ... (", n, " total)")
  }
}# Used in: MTD-methods.R, MTDest-methods.R and hdMTD-methods.R


###########################################################
###########################################################
###########################################################

# compute_transitP: Computes the global transition matrix of an MTD model
#
# This function constructs the global transition matrix P of an MTD model from
# the model parameters. The matrix P has |A|^|Lambda| rows (past contexts) and
# |A| columns (next-state values). Each row corresponds to a past context
# (ordered from oldest to newest, matching the ordering of Lambda) and gives
# the conditional distribution of the next state given that context.
#
# Arguments:
# - Lambda: A numeric vector of positive integers representing the relevant lag set.
#   The order of the elements must match the order of the list pj.
# - A: A vector representing the state space. The columns of P are named by A.
# - lambdas: Numeric vector of mixture weights of length length(Lambda) + 1.
#   The first element is the weight of the independent component (p0), and the
#   remaining elements correspond to each lag-specific transition matrix in pj.
# - pj: List of row-stochastic matrices (one per lag in Lambda). Each matrix must be
#   length(A) x length(A), with rows indexing the previous state and columns indexing
#   the next state.
# - p0: Numeric vector of length(A) giving the independent distribution on A. If the
#   independent component is inactive, p0 should be a zero vector (typically with
#   lambdas[1] = 0).
#
# Details:
# The returned transition probabilities are given by
#   P(a | x) = lambda_0 p0(a) + sum_{j in Lambda} lambda_j p_j(a | x_j),
# where x = (x_j)_{j in Lambda} denotes the past context.
#
# Returns:
# - A numeric matrix P with nrow(P) = length(A)^length(Lambda) and ncol(P) = length(A).
# - Columns are named by A.
# - If length(Lambda) > 1, row names encode past contexts by concatenating the values
#   in A from oldest to newest.
#
# Notes:
# - This function can be expensive when length(A)^length(Lambda) is large, since it
#   enumerates all past contexts via expand.grid().
compute_transitP <- function(Lambda, A, lambdas, pj, p0) {
  lenL <- length(Lambda)
  lenA <- length(A)
  lenAL <- lenA^lenL

  # Generate all possible size lenL sequences with digits from 1 to lenA
  subx <- try(expand.grid(rep(list(seq_len(lenA)), lenL)), silent = TRUE)
  if(inherits(subx,"try-error")) {
    stop(paste0("For length(Lambda)=",lenL," the dataset with all pasts sequences (x of length(Lambda)) with elements of A is too large."))
  }
  subx <- subx[, order(lenL:1)]

  P <- matrix(0, ncol = lenA, nrow = lenAL)

  if (lenL == 1) {
    for (i in seq_len(lenAL)) { # runs in all lines of P
      P[i, ] <- lambdas %*% rbind(p0, pj[[1]][i, ])
    }
    rownames(P) <- A
  } else {
    for (i in seq_len(lenAL)) {
      aux <- matrix(0, ncol = lenA, nrow = lenL)
      for (j in seq_len(lenL)) {
        aux[j, ] <- pj[[j]][subx[i, (lenL + 1 - j)], ]
        # The lines in aux are each from a different pj
      }
      P[i, ] <- lambdas %*% rbind(p0, aux)
    }
  }
  colnames(P) <- A
  if(lenL > 1){
    subx <- as.matrix(expand.grid(rep(list(A), lenL))) # Elements from A
    subx <- subx[, order(lenL:1)]
    rownames(P) <- apply(subx, 1, paste0, collapse = "")
  }
  return(P)
}  # Used in: MTDmodel.R and accessors.R - transitP().
