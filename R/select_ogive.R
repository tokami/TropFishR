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
#'   \item \code{Lc} length-at-first-capture (also called L50),
#'   \item \code{meshSizes} a vector with mesh sizes in increasing order,
#'   \item \code{select_p1} selectivity parameter 1 (see Millar and Holst (1997)),
#'   \item \code{select_p2} selectivity parameter 2 (see Millar and Holst (1997)),
#'   \item \code{L75} length at which individuals are caught with a
#'    probability of 75%.
#' }
#' @param Lt a vector with lengths corresponding to age classes
#' @param Lc length-at-first-capture (Default: NA)
#'
#' @examples
#' # create list with selectivity information
#' select.list <- list(selecType = 'knife_edge',
#'    Lc = 34, L75 = 37, tc = 5, meshSizes = c(60,80),
#'    select_p1 = 2.7977, select_p2 = 0.1175)
#'
#'
#' # create vector with mid lengths
#' Lt <- seq(5, 50, 0.01)
#'
#' # knife edge selectivity
#' sel_ke <- select_ogive(select.list, Lt)
#'
#' # trawl ogive selectivity
#' select.list$selecType = "trawl_ogive"
#' sel_to <- select_ogive(select.list, Lt)
#'
#' plot(Lt, sel_ke, type = 'l')
#' lines(Lt, sel_to, col = 'blue')
#'
#'
#' # Gillnet selectivity ("lognormal" and "normal_fixed")
#' select.list$selecType <- "lognormal"
#' sel_log <- select_ogive(select.list, Lt)
#'
#' select.list$selecType <- "normal_fixed"
#' select.list$select_p1 <- 0.2
#' select.list$select_p2 <- 1.5
#' sel_nf <- select_ogive(select.list, Lt)
#'
#' plot(Lt, sel_log, type = 'l')
#' lines(Lt, sel_nf, col = 'blue')
#'
#' @details This function is embedded within \code{\link{predict_mod}}. \code{selecType}
#'    "knife_edge" only requires a Lc value. "trawl_ogive" requires a Lc (L50) and
#'    a L75 value. "lognormal" requires two mesh sizes, an estimate of mu and of sigma.
#'    "normal_fixed" requires two mesh sizes with an estimate of the selection factor (SF) and an
#'    estimate of sigma.
#'
#' @references
#' Millar, R. B., Holst, R. (1997). Estimation of gillnet and hook selectivity using log-linear
#' models. ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
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
           }else L50 <- Lc
           sel <- rep(0, length(Lt))
           sel[as.numeric(as.character(Lt)) >= L50] <- 1
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
           mu <- s_list$select_p1
           sigma <- s_list$select_p2
           meshSizes <- s_list$meshSizes
           sel <- (1/Lt) * exp(mu + log(meshSizes[2]/meshSizes[1]) -
                                 (sigma^2/2) - (((log(Lt) - mu -
                                      log(meshSizes[2]/meshSizes[1]))^2)/(2 *sigma^2)))
         },

         "normal_fixed" = {
           SF <- s_list$select_p1
           sigma <- s_list$select_p2
           meshSizes <- s_list$meshSizes
           sel <- exp(-((Lt - meshSizes[1] * SF)^2/(2 * sigma^2)))
         },

         stop(paste("\n",type, "not recognised, possible curve types are \n",
                    "\"knife_edge\", \"trawl_ogive\", \"lognormal\" \n",
                    " and \"normal_fixed\""))
  )
  return(sel)
}

