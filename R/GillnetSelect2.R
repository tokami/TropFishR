#' @title Gillnet selctivity
#
#' @description  This function estimates the selecitvity of gillnets nets.
#'
#' @param param List with parameters: midlengths of size classes, number in net1, number in net2 and meshsizes of both nets
#'
#' @examples
#' data("data_GillnetSelect")
#' output <- GillnetSelect(data_GillnetSelect)
#' plot(output)
#'
#' @details Preconditions: selection curves of two mesh sizes must overlap, nets set in same area, during the same time. To calculate selection factor (SF), L25, L50 and L75 for gillnets. Assumptions: optimum length Lm is proportional to mesh size (Lm = SF * m), two selection curves have the same standard deviation, nets have the same fishing power (same dimensions and material). Assumption that selection curves are normally distributed with a common standard deviation.
#'
#' @references XX
#'
#'
#' @export


GillnetSelect <- function(param){

  res <- param
  classes <- as.character(res$midLengths)
  numNet1 <- res$numCa1
  numNet2 <- res$numCa2
  msNet1 <- res$msNet1
  msNet2 <- res$msNet2

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
  classes.num <- as.numeric(classes.num[,1])

  # log ratios y = ln(numNet2 / numNet1)
  lnNet2_Net1 <- log(numNet2 / numNet1)
  lnNet2_Net1[which(numNet1 == 0 | numNet2 == 0)] <- NA
  lnNet2_Net1[which(lnNet2_Net1 == Inf | lnNet2_Net1 == -Inf)] <- NA

  #regression analysis
  mod <- lm(lnNet2_Net1 ~ classes.num)
  sum.mod <- summary(mod)
  a <- sum.mod$coefficients[1]
  b <- sum.mod$coefficients[2]

  #plot
  plot(classes.num,lnNet2_Net1, main = "Regression analysis")
  abline(a=sum.mod$coefficients[1], b=sum.mod$coefficients[2])

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

  # final plot
  classes.num.plot <- classes.num - 0.5
  xL <- seq(min(classes.num,na.rm = T),max(classes.num,na.rm = T),0.01)


  if(max(numNet2,na.rm=T) > max(numNet1,na.rm=T)){
    XX <- numNet2
  }else  XX <- numNet1

  #create plot
  par(mar = c(5, 4, 4, 4) + 0.3)
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

  res2 <- list(classes.num=classes.num,
               SNet1=SNet1,
               SNet2=SNet2,
               LmNet1 = LmNet1,  #Lm optimum length
               LmNet2 = LmNet2,
               #L25 = L25,
               #L50 = L50,
               #L75 = L75,
               SF = SF,
               stand.dev = s)

  ret <- c(res,res2)
  class(ret) = 'GillnetSelect'
  return(ret)
}
