#' @title Year - Date conversion
#'
#' @description Convert numeric years to dates
#'
#' @param yeardec numeric year
#'
#' @examples
#'  yeardec2date(2014.14)
#'
#' @return date in format "%Y-%m-%d" (class 'Date').


yeardec2date <- function(yeardec){

  as.Date(strptime(paste(yeardec %/% 1,
                         ceiling(yeardec %% 1 * 365 + 1), sep="-"),
                   format = "%Y-%j"))

}

