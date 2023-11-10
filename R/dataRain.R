#' Rain data set for the city of Canberra, Australia.
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
#'   \item{Date}{The day}
#'   \item{RainToday}{1 if it rained that day, 0 otherwise.}
#' }
#'
#' @details Since this dataset was acquired with some missing observations,
#' certain actions were taken. The entire months of April 2011, December 2012,
#' and February 2013 were missing, and they were replaced with a repetition of
#' the same months from the previous year. Apart from these months, there
#' were a total of 18 days with missing data that where considered days without rain.
#'
#' @source {This data set was acquired at https://www.kaggle.com/datasets/jsphyg/weather-dataset-rattle-package
#' but it's original source is Data source: http://www.bom.gov.au/climate/dwo/ and http://www.bom.gov.au/climate/data.
#' Copyright Commonwealth of Australia 2010, Bureau of Meteorology.}
#'
#' @examples
#' data(raindata)
"raindata"
