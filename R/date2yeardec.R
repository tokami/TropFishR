#' @title Date - Year conversion
#'
#' @description Convert dates to numeric years with decimal as fraction of a year
#'
#' @param date a date (class 'Date')
#'
#' @examples
#' date2yeardec(Sys.Date())
#'
#' @return a scalar (class 'numeric')
#' @export

date2yeardec <- function(date){
  # adapted from lubridate::decimal_date.default
  if (any(!inherits(date, c("POSIXt", "POSIXct", "POSIXlt", "Date")))){
    stop("date(s) not in POSIXt or Date format")
  }
  if(inherits(date, "Date")){
    ## fast path (same result): day of the year / days of the year
    lt <- as.POSIXlt(date)
    Y <- lt$year + 1900
    leap <- (Y %% 4 == 0 & Y %% 100 != 0) | Y %% 400 == 0
    return(Y + (lt$yday * 86400) / ((365 + leap) * 86400))
  }
  Y <- as.numeric(format(date, format="%Y"))
  start <- as.POSIXct(paste0(Y,  "/01/01"), tz="UTC")
  end   <- as.POSIXct(paste0(Y+1,"/01/01"), tz="UTC")
  sofar <- as.numeric(difftime(date, start, units = "secs"))
  total <- as.numeric(difftime(end, start, units = "secs"))
  res <- Y + sofar/total
  return(res)
}
