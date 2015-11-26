#' @title Gillnet selctivity
#
#' @description  This function estimates the selecitvity of a gillnet from an
#'    experimental catch with two gillnets with different mesh sizes.
#'
#' @param param A list with following parameters: midlengths of size classes
#'    (\code{$midLengths}), number of fish caught with net 1 (\code{$numNet1}),
#'    number of fish caught with net 2 (\code{$numNet2}), and the meshsizes of both
#'    nets (\code{$msNet1} & \code{$msNet2}).
#'
#' @examples
#' # load data
#' data(tilapia)   ### == data_GillnetSelect
#'
#' # run model
#' output <- GillnetSelect(tilapia)
#'
#' # plot results
#' plot(output)
#'
#' @details This function estimates the fractions retained by each net (SNet1 & SNet2), the
#'   optimum lengths for each net, the selection factor (SF), and the standard deviation
#'   of the factor (stand.dev). Calculations are based on a normal distribution with common spread.
#'   Assumptions of this method are, that (i) the optimum length Lm is proportional to the mesh
#'   size (Lm = SF * m), (ii) the selection curves are normally distributed with a common
#'   standard deviation, (iii) the nets have the same fishing power (same dimensions and material).
#'   Requirements for the experimental set-up are: selection curves corresponding to the two
#'   mesh sizes have to overlap, and the nets have to be set in the same area, during the
#'   same time.
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

GillnetSelect <- function(param){

  res <- param
  classes <- as.character(res$midLengths)

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
  classes.num <- as.numeric(classes.num[,1])

  numNet1 <- res$numNet1
  numNet2 <- res$numNet2
  msNet1 <- res$msNet1
  msNet2 <- res$msNet2

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

  plot(classes.num.plot,XX,type='n', bty='n',xaxt ='n',
       ylab="numbers caught",xlab="fish length [cm]")
  lines(classes.num.plot,numNet1,type='s')
  axis(side=1,at=classes.num)
  lines(classes.num.plot,numNet2,type='s', lty=2)
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
  return(ret)
}
