# utils.R - Internal functions for the package
#
# This file contains auxiliary functions used within the package.
# These functions are not exported and are used internally to run some
# simple algorithm or calculation.

# 1 - groupTab()
# 2 - PI()
# 3 - is_xS()
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
  # Note. groupTab is used in: oscillation.R, hdMTD_FS.R.

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
  # Note. PI is used in hdMTD_FS.R.


  is_xS <- function(x,y) {
    return( all( x == y ) )
  }


