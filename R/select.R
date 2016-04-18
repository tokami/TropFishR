#' @title Selectivity model
#
#' @description  This function estimates the selecitvity of gillnets and
#'    trawl nets from experimental catches.
#'
#' @param param a list with following parameters:
#' \itemize{
#'   \item \code{type} a string indicating which type of gear was used (options:
#'   "gillnet" or "trawl_net"),
#'   \item \code{midLengths} midlengths of size classes
#'   \item \code{meshSizes} a vector with mesh sizes in increasing order,
#'   \item \code{CatchPerNet_mat} a matrix with the catches per net in corresponding
#'   order of mesh sizes
#' }
#' @param plot logical; should the results be displayed graphically?
#'
#' @examples
#' # Gillnet selectivity
#' data(tilapia)
#' select(param = tilapia)
#'
#' # Trawl selectivity
#' data(bream)
#' select(bream)
#'
#' @return A list with the input parameters and following list objects for type = "gillnet":
#' #' \itemize{
#'   \item \strong{classes.num}: numeric vector with length classes without a plus group,
#'   \item \strong{SNet1}: selection ogive net 1,
#'   \item \strong{SNet2}: selection ogive net 2,
#'   \item \strong{LmNet1}: optimum length net 1,
#'   \item \strong{LmNet2}: optimum length net 2,
#'   \item \strong{SF}: selection factor,
#'   \item \strong{stand.dev}: standard deviation of selection factor.
#' }
#' and following objects for type = "trawl_net":
#' \itemize{
#'   \item \strong{classes.num}: numeric vector with length classes,
#'   \item \strong{SLobs}: observed selection ogive,
#'   \item \strong{SLest}: estimated selection ogive,
#'   \item \strong{S1}: constant of selection curve,
#'   \item \strong{S2}: another constant of selection curve,
#'   \item \strong{L25}: length at which 25% of the fish are retained in the codend,
#'   \item \strong{L50}: length at which 50% of the fish are retained in the codend,
#'   \item \strong{L75}: length at which 75% of the fish are retained in the codend,
#'   \item \strong{SF}: selection factor.
#' }
#'
#' @details This function estimates the fractions retained by each net, the
#'   optimum lengths for each net, the selection factor (SF), and the standard deviation
#'   of the factor (stand.dev). Calculations are based on a normal distribution with common spread.
#'   Assumptions of this method are, that (i) the optimum length Lm is proportional to the mesh
#'   size (Lm = SF * m), (ii) the selection curves are normally distributed with a common
#'   standard deviation, (iii) the nets have the same fishing power (same dimensions and material).
#'   Requirements for the experimental set-up are: selection curves corresponding to the two
#'   mesh sizes have to overlap, and the nets have to be set in the same area, during the
#'   same time.
#'   To calculate selection factor (SF), L25, L50 and L75 for trawl nets /fisheries.
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

select <- function(param, plot = TRUE){
  res <- param
  classes <- as.character(res$midLengths)

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
  classes.num <- as.numeric(classes.num[,1])

  #  Gillnet
  if(is.null(res$type)) stop(noquote("Please indicate if a gillnet or trawl selectivity experiment was peroformed (see 'param$type')."))
  if(res$type == "gillnet"){
    numNet1 <- res$CatchPerNet_mat[,1]
    numNet2 <- res$CatchPerNet_mat[,2]
    msNet1 <- res$meshSizes[1]
    msNet2 <- res$meshSizes[2]

    # log ratios y = ln(numNet2 / numNet1)
    lnNet2_Net1 <- log(numNet2 / numNet1)
    lnNet2_Net1[which(numNet1 == 0 | numNet2 == 0 |
                        lnNet2_Net1 == Inf | lnNet2_Net1 == -Inf)] <- NA

    #regression analysis
    mod <- lm(lnNet2_Net1 ~ classes.num)
    sum.mod <- summary(mod)
    a <- sum.mod$coefficients[1]
    b <- sum.mod$coefficients[2]

    # SF
    SF <- (- 2 * a) / (b * (msNet1 + msNet2))
    LmNet1 <- SF * msNet1
    LmNet2 <- SF * msNet2

    # standard deviation
    s2 <- SF * ((msNet2 - msNet1) / b)
    s <- sqrt(s2)

    # points on selection curves
    SNet1 <- exp(-((classes.num - LmNet1)^2 / (2 * s2)))
    SNet2 <- exp(-((classes.num - LmNet2)^2 / (2 * s2)))

    # numbers in population
    NNet1 <- numNet1 / SNet1
    NNet2 <- numNet2 / SNet2

    res2 <- list(classes.num=classes.num,
                 SNet1=SNet1,
                 SNet2=SNet2,
                 lnNet2_Net1 = lnNet2_Net1,
                 reg.coeffs = sum.mod$coefficients,
                 LmNet1 = LmNet1,  #Lm optimum length
                 LmNet2 = LmNet2,
                 SF = SF,
                 stand.dev = s)
  }

  # Trawl
  if(res$type == "trawl_net"){
    numCodend <- res$CatchPerNet_mat[,2]
    numCover <- res$CatchPerNet_mat[,1]
    meshsizeCodend <- res$meshSizes[2]

    # calculate fraction retained (SL obs)
    SLobs <- numCodend/(numCodend + numCover)

    # ln( 1 / SL - 1)
    lnSL <- log(1/SLobs - 1)

    #excluding point where no (0) or full (1) retention was obtained and all values beyond those points, even if they are between 0 and 1
    lnSL[which(classes.num >= classes.num[which(SLobs == 1)])] <- NA
    lnSL[which(lnSL == Inf | lnSL == -Inf)] <- NA

    #model
    mod <- lm(lnSL ~ classes.num)
    sum.mod <- summary(mod)
    S1 <- sum.mod$coefficients[1]
    S2 <- abs(sum.mod$coefficients[2])

    #L25,L50,L75
    L25 <- (S1 - log(3)) / S2
    L50 <- S1 / S2
    L75 <- (S1 + log(3)) / S2

    #estimated SL (SL est)
    SLest <- 1 / (1 + exp(S1 - S2 * classes.num))

    #Selection factor
    SF <- L50/meshsizeCodend

    res2 <- list(classes.num = classes.num,
                 SLobs = SLobs,
                 SLest = SLest,
                 lnSL = lnSL,
                 reg.coeffs = sum.mod$coefficients,
                 S1 = S1,
                 S2 = S2,
                 L25 = L25,
                 L50 = L50,
                 L75 = L75,
                 SF = SF)
  }
  ret <- c(res,res2)
  class(ret) = 'select'
  if(plot) plot(ret)
  return(ret)
}
