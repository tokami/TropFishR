#' @title Date - Year conversion
#'
#' @description Convert dates to numeric years
#'
#' @param date a date (class 'Date')
#'
#' @examples
#' date2yeardec(Sys.Date())
#'
#' @return a scalar (class 'numeric')

date2yeardec <- function(date){
  as.POSIXlt(date)$year + 1900 + (as.POSIXlt(date)$yday) / 365
}

