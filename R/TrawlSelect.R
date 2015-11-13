#' @title Trawl net selctivity
#
#' @description  This function estimates the selecitvity of trawl nets.
#'
#' @param param List with parameters: midlengths of size classes, number in net, number in codend and meshsize of codend
#'
#' @examples
#' data(data_TrawlSelect)
#' output <- TrawlSelect(data_TrawlSelect)
#' plot(output)
#'
#' @details To calculate selection factor (SF), L25, L50 and L75 for trawl nets /fisheries.
#'
#' @references xx
#'
#'
#' @export

TrawlSelect <- function(param){

  res <- param
  classes <- as.character(res$midLengths)
  numCodend <- res$numCodend
  numCover <- res$numCover
  meshsizeCodend <- res$meshsizeCodend

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
  classes.num <- as.numeric(classes.num[,1])

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

  res2 <- list(classes.num=classes.num,
               SLobs=SLobs,
              SLest = SLest,
              S1 = S1,
              S2 = S2,
              L25 = L25,
              L50 = L50,
              L75 = L75,
              SF = SF)

  ret <- c(res,res2)
  class(ret) = 'TrawlSelect'
  return(ret)
}
