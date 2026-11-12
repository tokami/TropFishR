#' @title Year - Date conversion
#'
#' @description Convert numeric years to dates
#'
#' @param yeardec numeric year
#'
#' @examples
#'  yeardec2date(2014.14)
#'
#' @return date in format "\%Y-\%m-\%d" (class 'Date').
#'
#' @export

yeardec2date <- function(yeardec){
  # adapted from lubridate::date_decimal
  start <- as.POSIXct(paste0(trunc(yeardec),  "/01/01"), tz="UTC")
  end   <- as.POSIXct(paste0(trunc(yeardec)+1,"/01/01"), tz="UTC")
  secs <- as.numeric(difftime(end, start, units="secs")) * (yeardec - trunc(yeardec))
  ## round to milliseconds so that floating point error (e.g. 23:59:59.999998)
  ## does not move the date back by one day
  res <- as.Date(start + round(secs, 3))
  return(res)
}
