#' @title Beverton & Holt's Z-Equations
#
#' @description A method to estimate the instantaneous total mortality rate (Z) based
#'    on a method derived by Beverton and Holt (1956).
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item  \code{midLengths} or \code{age}: midpoints of length groups
#'   (length-frequency data) or ages (age composition data),
#'   \item \code{Linf}: infinite length for investigated species in cm [cm],
#'   \item \code{K}: growth coefficent for investigated species per year [1/year],
#'   \item \code{t0}: theoretical time zero, at which individuals of this species hatch,
#'   \item \code{catch}: catch as vector, or a matrix with catches of subsequent years;
#' }
#' @param catch_column optional; in case catch is a matrix or data.frame, a number
#'    indicating which column of the matrix should be analysed (Default: \code{NA}).
#'
#' @keywords function mortality Z
#'
#' @examples
#' # based on length-frequency data
#' data(synLFQ2)
#' Z_BevertonHolt(synLFQ2, catch_column = 2)
#'
#' # based on age composition data
#' data(synCAA1)
#' Z_BevertonHolt(synCAA1, catch_column = 3)
#'
#' @details  The first length group or age class within the list object \code{midLengths} or
#'    \code{age} will be used as the Lprim or tprime (length of recruitment to fishery).
#'
#' @return A list with the input parameters and following objects:
#' \itemize{
#'   \item \strong{tmean} or \strong{Lmean}: mean age or length of fish,
#'   \item \strong{tprime} or \strong{Lprime}: some age or length for which all fish of
#'      that length and longer are under full exploitation,
#'   \item \strong{Z}: total mortality.
#' }
#'
#'
#' @references
#' Beverton R.J.H and S.J. Holt, 1956. A review of methods of estimating mortality rates
#' in exploited fish populations, with special reference to sources of bias in catch
#' sampling. \emph{Rapp.P.-v.Reun.CIEM}, 140:67-83
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

Z_BevertonHolt <- function(param, catch_column = NA){
  res <- param
  catch <- res$catch

  if(class(catch) == "data.frame" | class(catch) == "matrix"){
    if(is.na(catch_column)) stop("Please provide a number indicating which column of the catch matrix should be analysed!")
    catch <- catch[,catch_column]
  }

  #   Length based equation
  if("midLengths" %in% names(res)){

    classes <- as.character(res$midLengths)
    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    Linf <- res$Linf
    K <- res$K

    # Error message if catch and age do not have same length
    if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
      if(length(classes) != length(catch[,1])) stop("Ages and catch do not have the same length!")
    }else if(class(catch) == 'numeric'){
      if(length(classes) != length(catch)) stop("Ages and catch do not have the same length!")
    }

    # calculate  C * (L1 + L2) / 2
    c_midlength <- catch * classes.num

    # calculate L mean
    Lmean <- sum(c_midlength) / sum(catch)

    # calculate L prime
    Lprime <- classes.num[1] -
      ((classes.num[2] - classes.num[1]) / 2 )

    Z = K * (Linf - Lmean) / (Lmean - Lprime)

    #save all in list
    ret <- c(res,list(
      Lmean = Lmean,
      Lprime = Lprime,
      Z = Z
    ))
    return(ret)
  }

  #     Aged based equation
  if("midAge" %in% names(res) | "age" %in% names(res)){

    if("midAge" %in% names(res)) classes <- as.character(res$midAge)
    if("age" %in% names(res)) classes <- as.character(res$age)
    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    # Error message if catch and age do not have same length
    if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
      if(length(classes) != length(catch[,1])) stop("Ages and catch do not have the same length!")
    }else if(class(catch) == 'numeric'){
      if(length(classes) != length(catch)) stop("Ages and catch do not have the same length!")
    }

    sample.size <- sum(catch,na.rm=T)
    sum.age.number <- sum((catch * classes.num), na.rm=T)
    tmean <- sum.age.number/sample.size
    interval <- (classes.num[2] - classes.num[1]) / 2
    tprime <- classes.num[1] - interval

    Z.BH <- 1 / (tmean - tprime)

    #save all in list
    ret <- c(res,list(
      tmean = tmean,
      tprime = tprime,
      Z = Z.BH
    ))
    return(ret)
  }

}
