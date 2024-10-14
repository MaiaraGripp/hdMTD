#' Maximum temperatures in the city of Brasília, Brazil.
#'
#' A data frame with the maximum temperature of the last
#' hour, by each hour, in the city of Brasília, Brazil.
#' The data spans from 01/01/2003 to 31/08/2024.
#'
#'
#' @format A data frame with 189936 rows and 3 columns. Each row corresponds
#' to a time in a day specified in columns 2 ("TIME") and 1 ("DATE")
#' respectively. The value in column 3 ("MAXTEMP")
#' is the maximum temperature measured in the last hour, in
#' Celsius (Cº), in the city of Brasília, the capital of Brazil, located in the
#' central-western part of the country.
#'
#' \describe{
#'   \item{DATE}{The day, from 01/01/2003 to 31/08/2024}
#'   \item{TIME}{The time, form 00:00 to 23:00 each day}
#'   \item{MAXTEMP}{The maximum temperature measured in the last hour in Celsius}
#' }
#'
#'
#' @source {This data set was acquired at INMET, the
#' National Institute of Meteorology in Brazil
#' https://bdmep.inmet.gov.br/. The measurements were made by
#' an automatic station in Brasília (latitude -15.79, longitude -470.93,
#'  altitude 1159.54)}
#'
#' @examples
#' data(tempdata)
"tempdata"
