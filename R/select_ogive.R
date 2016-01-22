#' @title Selectivity function for prediction models
#'
#' @description  This function estimates the fraction retained in the net for corresponding
#'    lengths of fish. Different selection curves can be used.
#'
#' @param s_list A list with a range of selectivity parameters dependent of the type of
#'    selection curve:
#' \itemize{
#'   \item \strong{$selecType} the type of
#'    selection curves to use for the calculation (options: "knife_edge",
#'    "trawl_ogive", "lognormal", "normal_fixed"),
#'   \item \strong{Lc} length-at-first-capture,
#'   \item \strong{meshSizes} a vector with mesh sizes in increasing order,
#'   \item \strong{select_ps} selectivity parameters,
#'   \item \strong{L75} length at which individuals are caught with a
#'    probability of 75%;
#' }
#' @param Lt A vector with lengths corresponding to age classes
#' @param Lc length-at-first-capture (Default: NA)
#'
#' @examples
#' # create list with selectivity information
#' select.list <- list(selecType = 'knife_edge',
#'    Lc = 34, tc = 5, mesh_size = 8.1, mesh_size1 = 9.1,
#'    select_p1 = 21.1, select_p2 = 23.8)
#'
#' # create vector with mid lengths
#' Lt <- c(12,22,29,34,38,41,43,45,46,47,48,49,49)
#'
#' # run model
#' select_ogive(select.list, Lt)
#'
#' @details Used within predicition models.
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

select_ogive <- function(s_list, Lt, Lc = NA){
  type <- s_list$selecType
  switch(type,
         'knife_edge' = {
           if(is.na(Lc)) Lc <- as.numeric(as.character(s_list$Lc))
           sel <- rep(0, length(Lt))
           sel[as.numeric(as.character(Lt)) >= Lc] <- 1
         },
         'trawl_ogive' = {
           if(is.na(Lc)) Lc <- s_list$L50
           L50 <-  Lc    #correct???
           #and L75 ???   -> same relationship between L50 and L75 ? then:
           L75 <-  Lc / (s_list$L50 / s_list$L75)

           sel <- 1 / (1 + exp(- (Lt - L50)/
                                 ((2*(L75 - L50))/(log(0.75/(1-0.75))-
                                                     log(0.25/(1-0.25))))))

           #         #alternative:
           #         sel <- 1 / (1 + exp(S1.TS - S2.TS * Lt))

         },
         "lognormal" = {
           sel <- (1/Lt) * exp(s_list$select_p1 +
                                 log(s_list$mesh_size/s_list$mesh_size1) -
                                 (s_list$select_p2^2/2) -
                                 (((log(Lt) - s_list$select_p1 -
                                      log(s_list$mesh_size/s_list$mesh_size1))^2)/
                                    (2 *s_list$select_p2^2)))

#            sel <- (1/classes) * exp(s_list$select_p1 +
#                                       log(s_list$mesh_size/s_list$mesh_size1) -
#                                       (s_list$select_p2^2/2) -
#                                       (((log(classes) -
#                                            s_list$select_p1 -
#                                            log(s_list$mesh_size/s_list$mesh_size1))^2)/
#                                          (2 *s_list$select_p2^2)))
         },
         "normal_fixed" = {
           sel <- exp(-((Lt - s_list$mesh_size *
                           s_list$select_p1)^2/
                          (2 * s_list$select_p2^2)))
         })
  return(sel)
}
