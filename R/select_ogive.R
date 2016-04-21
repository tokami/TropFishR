#' @title Selectivity patterns
#'
#' @description  Based on a few parameters, this function estimates the fraction per
#'    length group retained in the net. Different selection curves can be used for the
#'    estimation.
#'
#' @param s_list a list with selectivity parameters dependent on the type of
#'    selection curve:
#' \itemize{
#'   \item \code{selecType} type of
#'    selection curve used for estimation (options:
#'    "knife_edge",
#'    "trawl_ogive",
#'    "lognormal",
#'    "normal_fixed"),
#'   \item \code{Lc} length-at-first-capture,
#'   \item \code{meshSizes} a vector with mesh sizes in increasing order,
#'   \item \code{select_ps} selectivity parameters,
#'   \item \code{L75} length at which individuals are caught with a
#'    probability of 75%.
#' }
#' @param Lt a vector with lengths corresponding to age classes
#' @param Lc length-at-first-capture (Default: NA)
#'
#' @examples
#' # create list with selectivity information
#' select.list <- list(selecType = 'knife_edge',
#'    Lc = 34, tc = 5, meshSizes = c(8.1,9.1),
#'    select_ps = c(21.1,23.8))
#'
#' # create vector with mid lengths
#' Lt <- c(12,22,29,34,38,41,43,45,46,47,48,49,49)
#'
#' # run model
#' select_ogive(select.list, Lt)
#'
#' @details This function is embedded within \code{\link{predict_mod}}. \code{selecType}
#'    "knife_edge" only requires a Lc value. "trawl_ogive" requires a Lc (L50) and
#'    a L75 value. "lognormal" and "normal_fixed" require two mesh sizes with two
#'    corresponding "select_ps" values.
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

select_ogive <- function(s_list, Lt, Lc = NA){
  type <- s_list$selecType

  switch(type,

         'knife_edge' = {
           if(is.na(Lc)){
             if("L50" %in% names(s_list)){L50 <- s_list$L50
             }else L50 <- s_list$Lc
           }
           sel <- rep(0, length(Lt))
           sel[as.numeric(as.character(Lt)) >= Lc] <- 1
         },

         'trawl_ogive' = {
           # if no new Lc (L50) value is provided:
           if(is.na(Lc)){
             if("L50" %in% names(s_list)){L50 <- s_list$L50
             }else L50 <- s_list$Lc

             L75 <- s_list$L75
           }
           # if new Lc (L50) value is provided:
           if(!is.na(Lc)){
             if("L50" %in% names(s_list)){L50old <- s_list$L50
             }else L50old <- s_list$Lc
             L75old <- s_list$L75
             L50 <- Lc
             # new L75 based on relation of L50 values
             L75 <- L75old * (L50/L50old)
           }

           sel <- 1 / (1 + exp(- (Lt - L50)/
                                 ((2*(L75 - L50))/(log(0.75/(1-0.75))-
                                                     log(0.25/(1-0.25)))))) #or: sel <- 1 / (1 + exp(S1.TS - S2.TS * Lt))

         },

         "lognormal" = {
           select_ps <- s_list$select_ps
           meshSizes <- s_list$meshSizes
           sel <- (1/Lt) * exp(select_ps[1] +
                                 log(meshSizes[1]/meshSizes[2]) -
                                 (select_ps[2]^2/2) -
                                 (((log(Lt) - select_ps[1] -
                                      log(meshSizes[1]/meshSizes[2]))^2)/
                                    (2 *select_ps[2]^2)))
         },

         "normal_fixed" = {
           select_ps <- s_list$select_ps
           meshSizes <- s_list$meshSizes
           sel <- exp(-((Lt - meshSizes[1] *
                           select_ps[1])^2/
                          (2 * select_ps[2]^2)))
         },

         stop(paste("\n",type, "not recognised, possible curve types are \n",
                    "\"knife_edge\", \"trawl_ogive\", \"lognormal\" \n",
                    " and \"normal_fixed\""))
  )
  return(sel)
}

