#' @title Estimate Z from Catch per unit of effort (CPUE)
#'
#' @description Method to estimate the instantaneous total mortality rate (Z) from
#'    CPUE data according to standard, Heincke's, or Robson & Chapman's method.
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{cohort}: a vector with with a cohort label,
#'   \item \strong{age}: a vector with ages,
#'   \item \strong{CPUE}: a vector with Catch per unit of effort (CPUE) values;
#' }
#' @param method a character string indicating which assessment method should be used:
#'    \code{"standard"}, \code{"Heincke"}, or \code{"RobsonChapman"}.
#' @param omit_age1 logical; if \code{TRUE} the first age group is omitted (Default \code{FALSE}).
#'
#' @keywords function mortality CPUE
#'
#' @examples
#' # load data
#' data(synCPUE)
#'
#' # run model with standard method
#' Z_CPUE(synCPUE, method = "standard")
#'
#' # run model with Heincke's method
#' Z_CPUE(synCPUE, method = "Heincke")
#'
#' # run model with Robson and Chapman's method
#' Z_CPUE(synCPUE, method = "RobsonChapman", omit_age1 = TRUE)
#'
#' @details In Heincke's and RobsonChapman's method age groups older than 4 are lumped,
#'   because age groups older than 3 or 4 are said to be hard to seperate (Ricker, 1975).
#'   Sparre and Venema (1998) recommend to omit the first age group in case it is not fully
#'   exploited by the fishery.
#'
#' @return A list with input parameters and a Z value or matrix.
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' Sparre, P., Venema, S.C., 1999. Introduction to tropical fish stock assessment.
#' Part 2. Excercises. FAO Fisheries Technical Paper, (306.2, Rev. 2). 94 p.
#'
#' Ricker, W.E., 1975. Computation and interpretation of biological statistics of fish
#' populations. \emph{Bull.Fish.Res.Board Can.}, (191):382 p.
#'
#' @export

Z_CPUE <- function(param, method = "standard", omit_age1 = FALSE){

  res <- param
  cohort <- res$cohort
  classes <- as.character(res$age)
  CPUE <- res$CPUE

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
  classes.num <- as.numeric(classes.num[,1])

  switch(method,

         "standard" ={
           df.HZ <- data.frame(cohort = cohort[1:(length(cohort)-1)])
           result_Z <- list()
           for(i in 2:length(classes.num)){
             Zi <- rep(NA,(length(classes.num)-1))
             for(k in 1:(i-1)){
               Zi[k] <- round((1 / (classes.num[i] - classes.num[k])) *
                                (log(CPUE[k] /CPUE[i])),digits = 2)
             }
             result_Z[[i-1]] <- Zi
           }
           for(i in 1:length(result_Z)){
             df.HZ[[paste(cohort[i+1])]] <- result_Z[[i]]
           }
           ret <- c(res,list(
             Z_mat = df.HZ
           ))
           return(ret)
         },

         "Heincke" ={
           CPUE[4] <- sum(CPUE[4:length(CPUE)])
           if(omit_age1) CPUE <- CPUE[-1]
           CPUE.H.n <- CPUE[2:length(CPUE)]
           CPUE.H.d <- CPUE[1:length(CPUE)]
           Z.H = - log( (sum(CPUE.H.n)) /
                          (sum(CPUE.H.d)))
           ret <- c(res,list(
             Z = Z.H
           ))
           return(ret)
         },
         "RobsonChapman" ={
           CPUE[4] <- sum(CPUE[4:length(CPUE)])
           if(omit_age1) CPUE <- CPUE[-1]
           sum_CPUE.H.n <- sum(CPUE[2:length(CPUE)])
           sum_CPUE.H.d <- sum(CPUE[1:length(CPUE)])

           Z.H = - log((sum_CPUE.H.n) /
                          (sum_CPUE.H.d + sum_CPUE.H.n - 1))

           ret <- c(res,list(
             Z = Z.H
           ))
           return(ret)
         })
}
