#' @title Millar's selectivity types
#
#' @description  This function estimates the selecitvity of a gillnet from an experimental
#'    catch with two gillnets with different mesh sizes.
#'
#' @param param A list with following parameters: vector with midlengths of size classes
#'    (\code{$midLengths}), vector with meshSizesss in increasing order (\code{$meshSizesss}),
#'    and a matrix with the number of individuals caught with each sized mesh
#'    (\code{$CatchPerNet_mat}).
#' @param model A character string indicating which type fo model for the estimation
#'   of the selection curves should be used: \code{"fixed_normal"}, \code{"normal"},
#'   \code{"gamma"}, or \code{"lognormal"}.
#'
#' @examples
#'
#'
#' @source https://www.stat.auckland.ac.nz/~millar/selectware/
#'
#' @details Model adapted from the selectivity functions provided by Prof. Dr. Russell Millar
#'   (https://www.stat.auckland.ac.nz/~millar/).
#'
#'   This function estimates the fractions retained by each net (SNet1 and SNet2), the
#'   optimum lengths for each net, the selection factor (SF), and the standard deviation
#'   of the factor (stand.dev).
#'   Assumptions of this method are, that (i) the optimum length Lm is proportional to the mesh
#'   size (Lm equals SF times m), (ii) the selection curves are normally distributed with a common
#'   standard deviation, (iii) the nets have the same fishing power (same dimensions and material).
#'   Requirements for the experimental set-up are: selection curves corresponding to the two
#'   mesh sizes have to overlap, and the nets have to be set in the same area, during the
#'   same time.
#'
#' @references
#'  Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
#'  Holt, S. J. 1963. A method for determining gear selectivity and its application.
#'  ICNAF Special Publication, 5: 106-115.
#'
#'
#' @export



#Curves to be added include:
#tt.richards, for richards fit to trouser trawl data
#gamma, for net selectivity.

rtypes_Millar <- function(rtype) {
  switch(rtype,
         "norm.loc"={
           r=function(classes,meshSizes,th) {
             relsize=meshSizes/meshSizes[1]
             seln=exp(-(classes-th[1]*relsize)^2/(2*th[2]^2))
             return(seln) } },
         "norm.sca"={
           r=function(classes,meshSizes,th) {
             relsize=meshSizes/meshSizes[1]
             seln=exp(-(classes-th[1]*relsize)^2/(2*th[2]^2*relsize^2))
             return(seln) } },
         "lognorm"={
           r=function(classes,meshSizes,th) {
             relsize=meshSizes/meshSizes[1]
             seln=(relsize/classes)*exp(th[1]-th[2]^2/2)
             seln=seln*exp( -(log(classes)-th[1]-log(relsize))^2/(2*th[2]^2) )
             return(seln) } },
         "binorm.sca"={
           r=function(classes,meshSizes,th) {
             relsize=meshSizes/meshSizes[1]
             seln1=exp(-(classes-th[1]*relsize)^2/(2*th[2]^2*relsize^2))
             seln2=exp(-(classes-th[3]*relsize)^2/(2*th[4]^2*relsize^2))
             p=exp(th[5])/(1+exp(th[5])) #i.e., th[5]=logit(p)
             seln=p*seln1+(1-p)*seln2
             return(seln) } },
         "bilognorm"={
           r=function(classes,meshSizes,th) {
             relsize=meshSizes/meshSizes[1]
             seln1=(relsize/classes)*exp(th[1]-th[2]^2/2)
             seln1=seln1*exp( -(log(classes)-th[1]-log(relsize))^2/(2*th[2]^2) )
             seln2=(relsize/classes)*exp(th[3]-th[4]^2/2)
             seln2=seln2*exp( -(log(classes)-th[3]-log(relsize))^2/(2*th[4]^2) )
             p=exp(th[5])/(1+exp(th[5])) #i.e., th[5]=logit(p)
             seln=p*seln1+(1-p)*seln2
             return(seln) } },
         "tt.logistic"={
           r=function(classes,meshSizes,th) {
             control=(meshSizes==meshSizes[1])
             p=exp(th[3])/(1+exp(th[3])) #i.e., th[3]=logit(p)
             wk=exp(th[1]+th[2]*classes)
             lselect=wk/(1+wk)
             seln=(1-p)*control+p*lselect*(1-control)
             return(seln) } },
         stop(paste("\n",rtype, "not recognised, possible curve types are \n",
                    "\"norm.loc\", \"norm.sca\", \"lognorm\" \n",
                    "\"binorm.sca\", \"bilognorm\", and \"tt.logistic\""))
  )#End of switch
  return(r) }

