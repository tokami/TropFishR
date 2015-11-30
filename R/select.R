#' @title Selectivity model
#
#' @description  This function estimates the selecitvity of different gillnets and trawl nets
#'    from experimental catches.
#'
#' @param param A list with following parameters: midlengths of size classes
#'    (\code{$midLengths}), a string indicating which type of gear was used (\code{$type}) so far
#'    \code{"gillnet"} or \code{"trawl_net"} are possible types, a vector with increasing mesh sizes (\code{$meshSizes}),
#'    and a matrix with the catches per net in corresponding order of mesh sizes (\code{$CatchPerNet_mat}).
#'
#' @examples
#' # Gillnet selectivity
#' # load data
#' data(tilapia)
#'
#' # run model
#' output <- select(tilapia)
#'
#' # plot results
#' #plot(output)
#'
#' # Trawl selectivity
#' # load data
#' data(bream)
#'
#' # run model
#' output <- select(bream)
#'
#' # plot results
#' #plot(output)
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
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export


select <- function(param){

  res <- param
  classes <- as.character(res$midLengths)

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
  classes.num <- as.numeric(classes.num[,1])

  #HHHHHHHHHHH#
  #  GILLNET  #
  #HHHHHHHHHHH#
  if(res$type == "gillnet"){  #if("msCodend" %in% names(res) == FALSE){
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

    #   #L25,L50,L75
    #   L25 <- (S1 - log(3)) / S2
    #   L50 <- S1 / S2
    #   L75 <- (S1 + log(3)) / S2

    # points on selection curves
    SNet1 <- exp(-((classes.num - LmNet1)^2 / (2 * s2)))
    SNet2 <- exp(-((classes.num - LmNet2)^2 / (2 * s2)))

    # numbers in population
    NNet1 <- numNet1 / SNet1
    NNet2 <- numNet2 / SNet2

    # final plot
    classes.num.plot <- classes.num - 0.5
    xL <- seq(min(classes.num,na.rm = T),max(classes.num,na.rm = T),0.01)


    if(max(numNet2,na.rm=T) > max(numNet1,na.rm=T)){
      XX <- numNet2
    }else  XX <- numNet1

    #create plot
    op <- par(mfrow=c(2,1), mar=c(4,4,2,4),oma=c(1,1,1,1))
    plot(classes.num,lnNet2_Net1)
    abline(a=sum.mod$coefficients[1], b=sum.mod$coefficients[2])

    plot(classes.num.plot, XX, type='n', bty='n',xaxt ='n', ylim = c(0,max(XX,ra.rm = TRUE)),
         ylab="numbers caught", xlab="fish length [cm]")
    lines(classes.num.plot, numNet1, type='s')
    axis(side=1, at=classes.num)
    lines(classes.num.plot, numNet2, type='s', lty=2)
    par(new=TRUE)
    plot(xL, exp(- ((xL - LmNet1)^2 / (2 * s2))), type = "l", col='darkgreen',lwd=2.5,
         axes=F,bty = "n", xlab = "", ylab = "")
    lines(xL, exp(- ((xL - LmNet2)^2 / (2 * s2))), type = "l", col='orange',lwd=2.5,
          bty = "n", xlab = "", ylab = "")
    axis(side=4, at = pretty(range(SNet1)))
    mtext("fractions retained", side=4, line=3)
    text(labels = paste("ms = ",msNet1,"cm"), x = LmNet1, y =
           ((exp(- ((LmNet1 - LmNet1)^2 / (2 * s2))))/1.1),
         col = 'darkgreen')
    text(labels = paste("ms = ",msNet2,"cm"), x = LmNet2, y =
           ((exp(- ((LmNet2 - LmNet2)^2 / (2 * s2))))/1.1),
         col = 'orange')
    par(op)

    res2 <- list(classes.num=classes.num,
                 SNet1=SNet1,
                 SNet2=SNet2,
                 LmNet1 = LmNet1,  #Lm optimum length
                 LmNet2 = LmNet2,
                 SF = SF,
                 stand.dev = s)

    ret <- c(res,res2)
    class(ret) = 'GillnetSelect'
  }


  #HHHHHHHHHHH#
  #   TRAWL   #
  #HHHHHHHHHHH#
  if(res$type == "trawl_net"){  #if("msCodend" %in% names(res) == TRUE){
    numCodend <- res$CatchPerNet_mat[,2]
    numCover <- res$CatchPerNet_mat[,1]
    meshsizeCodend <- res$meshSizes[2]   #check, so far works just for two meshsizes

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

    #plot
    op <- par(mfrow=c(2,1), mar=c(4,4,2,2),oma=c(1,1,1,1))
    plot(classes.num,lnSL)#, main = "Regression analysis")
    abline(a=sum.mod$coefficients[1], b=sum.mod$coefficients[2])

    xL <- seq(min(classes.num,na.rm=TRUE),max(classes.num,na.rm=TRUE),0.1)
    plot(SLobs ~ classes.num,
         xlab = "L [cm]", ylab = "fraction retained", pch = 16 , cex =1.4)
    lines(x = xL, y = 1/(1 + exp(S1 - S2 * xL)), col='blue', lwd = 1.5)
    segments(x0 = L25, x1 = L25, y0 = 0,
             y1 = 1/(1 + exp(S1 - S2 * L25)),col= 'gray40',lty = 2, lwd=1.5)
    segments(x0 = 0, x1 = L25, y0 = 1/(1 + exp(S1 - S2 * L25)),
             y1 = 1/(1 + exp(S1 - S2 * L25)),col= 'gray40',lty = 2, lwd=1.5)
    text(labels = "L25%", x = L25, y = (1/(1 + exp(S1 - S2 * L25)))/3,
         col = 'gray40')
    segments(x0 = L50, x1 = L50, y0 = 0,
             y1 = 1/(1 + exp(S1 - S2 * L50)),col= 'gray40',lty = 2, lwd=1.5)
    segments(x0 = 0, x1 = L50, y0 = 1/(1 + exp(S1 - S2 * L50)),
             y1 = 1/(1 + exp(S1 - S2 * L50)),col= 'gray40',lty = 2, lwd=1.5)
    text(labels = "L50%", x = L50, y = (1/(1 + exp(S1 - S2 * L50)))/3,
         col = 'gray40')
    segments(x0 = L75, x1 = L75, y0 = 0,
             y1 = 1/(1 + exp(S1 - S2 * L75)),col= 'gray40',lty = 2, lwd=1.5)
    segments(x0 = 0, x1 = L75, y0 = 1/(1 + exp(S1 - S2 * L75)),
             y1 = 1/(1 + exp(S1 - S2 * L75)),col= 'gray40',lty = 2, lwd=1.5)
    text(labels = "L75%", x = L75, y = (1/(1 + exp(S1 - S2 * L75)))/3,
         col = 'gray40')
    par(op)

    res2 <- list(classes.num = classes.num,
                 SLobs = SLobs,
                 SLest = SLest,
                 S1 = S1,
                 S2 = S2,
                 L25 = L25,
                 L50 = L50,
                 L75 = L75,
                 SF = SF)

    ret <- c(res,res2)
    class(ret) = 'TrawlSelect'
  }
  return(ret)
}
