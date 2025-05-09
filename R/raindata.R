#' Rain data set for the city of Canberra, Australia
#'
#' A data frame with the rainfall history in the city of Canberra, Australia.
#' The data spans from 01/11/2007 to 25/06/2017.
#'
#'
#' @format A data frame with 3525 rows and 2 columns. Each row corresponds
#' to a day specified in column 1 ("Date"). The value in column 2 ("RainToday")
#'  is 0 if no rain was recorded in the city of Canberra that day, and 1 otherwise.
#'
#' \describe{
#'   \item{Date}{Date in YYYY-MM-DD format.}
#'   \item{RainToday}{Binary indicator (0 = no rain, 1 = rain).}
#' }
#'
#' @source
#' **Original data source**: Australian Bureau of Meteorology (BOM).
#' **Accessed via**: Kaggle (\url{https://www.kaggle.com/datasets/jsphyg/weather-dataset-rattle-package}).
#'
#' @note
#' The original BOM data is subject to their Terms of Use.
#' For direct access, visit the BOM website manually:
#' \url{https://www.bom.gov.au/climate/data/} (may require browser access).
#'
#' @examples
#' data(raindata)
"raindata"
