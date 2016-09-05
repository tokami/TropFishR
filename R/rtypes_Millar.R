#' @title Millar's selectivity types
#
#' @description  This function returns a function corresponding to the type of curve
#'    which was selected to represent the selectivity of nets or hooks.
#'
#' @param rtype a character string indicating which method for the estimation of selection
#'    curves should be used:
#'    \code{"norm.loc"} for normal with common spread method,
#'    \code{"norm.sca"} for normal with variable spread method,
#'    \code{"lognorm"} for lognormal method,
#'    \code{"binorm.sca"} for bi-normal method,
#'    \code{"bilognorm"} for bi-lognormal method,
#'    \code{"tt.logistic"} for control and logistic method,
#'    \code{"gamma"} for gamma method.
#'
#' @source https://www.stat.auckland.ac.nz/~millar/selectware/
#'
#' @details Function adapted from the selectivity functions provided by
#'   Prof. Dr. Russell Millar (https://www.stat.auckland.ac.nz/~millar/).
#'   Until now following curves are incorporated:
#'   \code{"norm.loc"} for a normal curve with common spread,
#'    \code{"norm.sca"} for a normal curve with variable spread,
#'    \code{"lognorm"} for a lognormal curve,
#'    \code{"binorm.sca"} for a bi-normal curve,
#'    \code{"bilognorm"} for a bi-lognormal curve,
#'    \code{"tt.logistic"} for a control and logistic curve.
#'
#' @references
#'  Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#'  54(3), 471-477.

rtypes_Millar <- function(rtype) {
# Adapted R code from Russell Millar (https://www.stat.auckland.ac.nz/~millar/selectware/)
  switch(rtype,

         "norm.loc" = {
           r <- function(classes, meshSizes, th) {
             relsize <- meshSizes / meshSizes[1]
             seln <- exp(-(classes - th[1] * relsize)^2 / (2 * th[2]^2))
             return(seln)
           }
         },

         "norm.sca" = {
           r <- function(classes, meshSizes, th) {
             relsize <- meshSizes / meshSizes[1]
             seln <- exp(-(classes - th[1]*relsize)^2 / (2 * th[2]^2 * relsize^2))
             return(seln)
             }
           },

         "lognorm" = {
           r <- function(classes, meshSizes, th) {
             relsize <- meshSizes / meshSizes[1]
             seln <- (relsize / classes) * exp(th[1] - th[2]^2 / 2)
             seln <- seln * exp( -(log(classes) - th[1] - log(relsize))^2 / (2*th[2]^2))
             return(seln)
             }
           },

         "binorm.sca" = {
           r <- function(classes, meshSizes, th) {
             relsize <- meshSizes / meshSizes[1]
             seln1 <- exp(-(classes - th[1] * relsize)^2 / (2 * th[2]^2 * relsize^2))
             seln2 <- exp(-(classes - th[3] * relsize)^2 / (2 * th[4]^2 * relsize^2))
             p <- exp(th[5]) / (1 + exp(th[5])) #i.e., th[5]=logit(p)
             seln <- p * seln1 + (1 - p) * seln2
             return(seln)
             }
           },

         "bilognorm" = {
           r <- function(classes, meshSizes, th) {
             relsize <- meshSizes / meshSizes[1]
             seln1 <- (relsize / classes) * exp(th[1] - th[2]^2 / 2)
             seln1 <- seln1*exp( -(log(classes) - th[1] - log(relsize))^2 / (2 * th[2]^2))
             seln2 <- (relsize / classes) * exp(th[3] - th[4]^2 / 2)
             seln2 <- seln2 * exp( -(log(classes) - th[3] - log(relsize))^2 / (2 * th[4]^2))
             p <- exp(th[5]) / (1 + exp(th[5])) #i.e., th[5]=logit(p)
             seln <- p * seln1 + (1 - p) * seln2
             return(seln)
             }
           },

         "tt.logistic" = {
           r <- function(classes, meshSizes, th) {
             control <- (meshSizes == meshSizes[1])
             p <- exp(th[3]) / (1 + exp(th[3])) #i.e., th[3]=logit(p)
             wk <- exp(th[1] + th[2] * classes)
             lselect <- wk / (1 + wk)
             seln <- (1 - p) * control + p * lselect * (1 - control)
             return(seln)
             }
           },
         "gamma" = {
           r <- function(classes, meshSizes, th){
             seln <- (classes / ((th[1] - 1) * th[2] * meshSizes)) ^ (th[1] - 1) *
               exp(th[1] - 1 - classes / (th[2] * meshSizes))
             return(seln)
           }
         },

         stop(paste("\n",rtype, "not recognised, possible curve types are \n",
                    "\"norm.loc\", \"norm.sca\", \"lognorm\" \n",
                    "\"binorm.sca\", \"bilognorm\", \"tt.logistic\" and \"gamma\""))
  )
  return(r)
}


