# utils.R - Internal functions for the package
#
# This file contains auxiliary functions used within the package.
# These functions are not exported and are used internally to run some
# simple algorithm or calculation.

# 1 - groupTab()
# 2 - PI()
# 3 - sx()
# 4 - n_parameters()

# At the end of each auxiliary function bellow there is
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

  groupTab <- function(S,j,freqTab,lenX,d){

    Sj <- sort(c(S,j), decreasing = TRUE) # The set of lags

    if ( length(Sj) > 0 ) {
      # Summarizes freqTab by lags in Sj
      groupTab <- freqTab %>%
        dplyr::group_by_at(paste0("x", Sj)) %>%
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

# PI: Estimates the empirical stationary distribution for a given sequence.
#
# This function computes the stationary distribution for sequences stored in a
# frequency table (`groupTab`). It filters sequences that match `x_S` in the
# specified lags `S`, then normalizes their frequencies to estimate stationary
# probabilities.
#
# Arguments:
# - S: A numeric vector of past lags. Determines which columns in `groupTab` should
#   be used for filtering.
# - groupTab: A tibble containing sequence frequencies (`Nx_Sj` column).
# - x_S: A vector representing a specific sequence of states in lags `S`.
# - lenX: Sample size (integer).
# - d: Maximum lag order (integer).
#
# Returns:
# - A numeric matrix (column vector) with estimated stationary probabilities.
#   The column name corresponds to the concatenated elements of `x_S`.

  PI <- function(S, groupTab, x_S, lenX, d) {

    if ( length(S) > 0 ) {
    # Filters groupTab by x_S.
      filtr_S <- paste0("x", S)
      groupTab <- groupTab %>%
        dplyr::mutate(match = purrr::pmap_lgl(dplyr::pick(dplyr::all_of(filtr_S)), ~ all(c(...) == x_S))) %>%
        dplyr::filter(match) %>%
        dplyr::select(-match)
    }
    PI <- matrix(groupTab$Nx_Sj/(lenX-d),ncol = 1)
    colnames(PI) <- "freq"
    PI
  }
  # PI is used in hdMTD_FS.R.

###########################################################
###########################################################
###########################################################

# sx: Computes thresholds for the CUT method in MTD models
#
# This function calculates the threshold used in the CUT step of the MTD inference algorithm.
# It computes a quantity that determines whether the variation in distributions across lagged states
# is significant, based on the provided parameters `mu`, `alpha`, and `xi`.
#
# Arguments:
# - S: A numeric vector of past lags. Determines which columns in `freqTab` should be used for filtering.
# - freqTab: An output of freqTab() function, i.e. a tibble containing frequency counts
#  (`Nx_Sj` column) and conditional probabilities (`qax_Sj`).
# - lenA: The number of distinct states in the state space.
# - x_S: A vector representing a specific sequence of states in lags `S`.
# - mu: A positive real number between 0 and 3, influencing the threshold calculation.
# - alpha: A positive real number controlling the sensitivity of the threshold.
# - xi: A positive real number affecting the threshold scaling.
#
# Returns:
# - A numeric vector of length `lenA`, where each entry represents the computed threshold
# for a specific state.


  sx <- function(S, freqTab, lenA, x_S, mu, alpha, xi){
    filtr_S <- paste0("x",S)
    C <- freqTab %>%
      mutate(match = purrr::pmap_lgl(pick(all_of(filtr_S)), ~ all(c(...) == x_S))) %>%
      filter(match) %>%
      select(-match)

    Nx <- C$Nx_Sj[seq(1, nrow(C), by = lenA)]

    prob_adjusted <- sqrt((C$qax_Sj + alpha / rep(Nx, each = lenA)) * mu / (2 * mu - exp(mu) + 1))
    sum <- colSums(matrix(prob_adjusted, nrow = lenA))

    return(sqrt(0.5 * alpha * (1 + xi) / Nx) * sum + alpha * lenA / (6 * Nx))
  } # sx is used in hdMTD_CUT.R.

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

  n_parameters <- function(Lambda, A,
                           single_matrix = FALSE,
                           indep_part = TRUE,
                           zeta = ifelse(single_matrix, 1, length(Lambda))){
    lenA <- length(A)
    lenL <- length(Lambda)

    n_parameters <- lenL-1 # Number of free weight parameters if lam0 = 0
    if (indep_part) {
      n_parameters <- n_parameters + lenA  # adds 1 ( for lam0) plus lenA - 1 ( for p0)
    }
    n_parameters <- n_parameters + lenA * (lenA - 1) * zeta
    # lenA*(lenA-1) is the number of free parameters in each matrix pj. zeta is the number of distinct matrices pj
    n_parameters
  } # n_parameters is used in hdMTD_BIC.R.


###########################################################
###########################################################
###########################################################

# prodinf: Computes the vector product between to vectors
#
# This function computes the product vector between to vectors x and y.
# However, for each i such that |x[i]|=inf and y[i]=0, x[i] * y[i] = 0.
#
# Arguments:
# - x: A numeric vector, may have inf values
# - y: A numeric vector with the same length as x.
#
# Returns:
# - A size x vector with the product between x and y. Whenever inf*0 occurs the output
# is set to 0.

  prodinf <- function(x,y){
    prinf <- numeric(length(x))
    for ( i in seq_len(length(x)) ) {
      if (is.infinite(x[i]) && y[i] == 0) {
        prinf[i] <- 0
      } else {
        prinf[i] <- x[i] * y[i]
      }
    }
    prinf
  }


