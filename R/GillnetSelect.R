#' @title Gillnet selctivity
#
#' @description  This function estimates the selecitvity of gillnets nets.
#'
#' @param classes Midpoints of the length class as vector (length frequency data) or ages as vector (age composition data).
#' @param catch Catch as vector, or a matrix with catches of subsequent years if the catch curve with constat time intervals should be applied.
#' @param datatype Type of data which is used for analysis, either 'length' or 'age', for length frequency or age composition data, respectively
#' @param Linf Infinite length for investigated species in cm [cm].
#' @param K Growth coefficent for investigated species per year [1/year].
#' @param t0 Theoretical time zero, at which individuals of this species hatch (default: 0).
#'
#' @examples
#' data("ex.GillnetSelect")
#' with(ex.GillnetSelect,GillnetSelect(midLengths))
#'
#'
#' @details To calculate selection factor (SF), L25, L50 and L75 for trawl nets /fisheries.
#'
#' @references
#'
#'
#' @export


GillnetSelect <- function(classes){

  df.TS <- cbind(classes,numCodend,numCover)
  df.TS <- as.data.frame(df.TS)
  df.TS$classes <- as.character(df.TS$classes)

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(df.TS$classes, split="\\+"))
  df.TS$classes.num <- as.numeric(classes.num[,1])

  #


}
