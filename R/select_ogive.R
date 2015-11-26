#' @title Selectivity function
#
#' @description  This function estimates the selecitvity.
#'    First sentence. second sentence.
#'
#' @param s_list list of selectivity parameters
#' @param classes age or size classes
#' @param Lt lengths corresponding to age classes
#' @param Lc length-at-first-capture Default = \code{NA}
#'
#' @examples
#' # create list with selectivity information
#' select.list <- list(selecType = 'knife_edge',  #or 'gillnet' or 'trawl_ogive'
#'    Lc = 34,tc = 5,selecDist = 'lognormal',    #or 'normal_fixed'
#'    mesh_size = 8.1,mesh_size1 = 9.1,select_p1 = 21.1,select_p2 = 23.8)
#'
#' # create vector with mid lengths
#' Lt <- c(12,22,29,34,38,41,43,45,46,47,48,49,49)
#'
#' # run model
#' select_ogive(select.list, seq(1,13,1),Lt)
#'
#' @details To calculate selection factor (SF), L25, L50 and L75 for trawl nets /fisheries.
#'
#' @references xx
#'
#'
#' @export

select_ogive <- function(s_list, classes, Lt, Lc = NA){

  if(s_list$selecType == 'knife_edge'){
    if(is.na(Lc)) Lc <- s_list$Lc
    sel <- rep(0, length(Lt))
    sel[Lt >= Lc] <- 1

  }else if(s_list$selecType == 'trawl_ogive'){
    if(is.na(Lc)) Lc <- s_list$L50
    L50 <-  Lc    #correct???
    #and L75 ???   -> same relationship between L50 and L75 ? then:
    L75 <-  Lc / (s_list$L50 / s_list$L75)

    sel <- 1 / (1 + exp(- (Lt - L50)/
                          ((2*(L75 - L50))/(log(0.75/(1-0.75))-
                                              log(0.25/(1-0.25))))))

    #         #alternative:
    #         sel <- 1 / (1 + exp(S1.TS - S2.TS * Lt))

  }else if (s_list$selecType == "lognormal") {
    sel <- (1/classes) * exp(s_list$select_p1 +
                                             log(s_list$mesh_size/s_list$mesh_size1) -
                                             (s_list$select_p2^2/2) -
                                             (((log(classes) -
                                                  s_list$select_p1 -
                                                  log(s_list$mesh_size/s_list$mesh_size1))^2)/
                                                (2 *s_list$select_p2^2)))
  }
  if (s_list$selecType == "normal_fixed") {
    sel <- exp(-((classes - s_list$mesh_size *
                    s_list$select_p1)^2/
                   (2 * s_list$select_p2^2)))
  }
  return(sel)
}
