#' @title Gillnet selctivity
#
#' @description  This function estimates the selecitvity of gillnets nets.
#'
#' @param classes Midpoints of the length class as vector (length frequency data) or ages as vector (age composition data).
#' @param catch Catch as vector, or a matrix with catches of subsequent years if the catch curve with constat time intervals should be applied.
#' @param numNet1 Number of individuals caught using the smaller meshed net
#' @param numNet2 Number of individuals caught using the larger meshed net
#' @param msNet1 mesh size of smaller meshed net
#' @param msNet2 mesh size of larger meshed net
#'
#' @examples
#' data("ex.GillnetSelect")
#' with(ex.GillnetSelect,GillnetSelect(midLengths,numCa1,numCa2,
#'    msNet1 = 8.1,msNet2 = 9.1))
#'
#'
#' @details Preconditions: selection curves of two mesh sizes must overlap, nets set in same area, during the same time. To calculate selection factor (SF), L25, L50 and L75 for gillnets. Assumptions: optimum length Lm is proportional to mesh size (Lm = SF * m), two selection curves have the same standard deviation, nets have the same fishing power (same dimensions and material). Assumption that selection curves are normally distributed with a common standard deviation.
#'
#' @references
#'
#'
#' @export

GillnetSelect <- function(classes, numNet1, numNet2, msNet1, msNet2){

  df.GS <- cbind(classes,numNet1,numNet2)
  df.GS <- as.data.frame(df.GS)
  df.GS$classes <- as.character(df.GS$classes)

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(df.GS$classes, split="\\+"))
  df.GS$classes.num <- as.numeric(classes.num[,1])

  # log ratios y = ln(numNet2 / numNet1)
  df.GS$lnNet2_Net1 <- log(df.GS$numNet2 / df.GS$numNet1)
  df.GS$lnNet2_Net1[which(df.GS$numNet1 == 0 | df.GS$numNet2 == 0)] <- NA
  df.GS$lnNet2_Net1[which(df.GS$lnNet2_Net1 == Inf | df.GS$lnNet2_Net1 == -Inf)] <- NA

  #regression analysis
  mod.GS <- lm(lnNet2_Net1 ~ classes.num, df.GS)
  sum.mod.GS <- summary(mod.GS)
  a.GS <- sum.mod.GS$coefficients[1]
  b.GS <- sum.mod.GS$coefficients[2]

  #plot
  plot(df.GS$classes.num,df.GS$lnNet2_Net1, main = "Regression analysis")
  abline(a=sum.mod.GS$coefficients[1], b=sum.mod.GS$coefficients[2])
  plot1 = recordPlot()

  # SF
  SF.GS <- (- 2 * a.GS) / (b.GS * (msNet1 + msNet2))
  LmNet1 <- SF.GS * msNet1
  LmNet2 <- SF.GS * msNet2

  # standard deviation
  s2.GS <- SF.GS * ((msNet2 - msNet1) / b.GS)
  s.GS <- sqrt(s2.GS)

  # points on selection curves
  df.GS$SNet1 <- exp(-((df.GS$classes.num - LmNet1)^2 / (2 * s2.GS)))
  df.GS$SNet2 <- exp(-((df.GS$classes.num - LmNet2)^2 / (2 * s2.GS)))

  # numbers in population
  df.GS$NNet1 <- df.GS$numNet1 / df.GS$SNet1
  df.GS$NNet2 <- df.GS$numNet2 / df.GS$SNet2

  # final plot
  df.GS$classes.num.plot <- df.GS$classes.num - 0.5
  xL <- seq(min(df.GS$classes.num,na.rm = T),max(df.GS$classes.num,na.rm = T),0.01)


  if(max(df.GS$numNet2,na.rm=T) > max(df.GS$numNet1,na.rm=T)){
    XX <- df.GS$numNet2
  }else  XX <- df.GS$numNet1

  #create plot
  par(mar = c(5, 4, 4, 4) + 0.3)
  plot(df.GS$classes.num.plot,XX,type='n', bty='n',xaxt ='n',
       ylab="numbers caught",xlab="fish length [cm]")
  lines(df.GS$classes.num.plot,df.GS$numNet1,type='s')
  axis(side=1,at=df.GS$classes.num)
  lines(df.GS$classes.num.plot,df.GS$numNet2,type='s', lty=2)
  par(new=TRUE)
  plot(xL, exp(- ((xL - LmNet1)^2 / (2 * s2.GS))), type = "l", col='darkgreen',lwd=2.5,
        axes=F,bty = "n", xlab = "", ylab = "")
  lines(xL, exp(- ((xL - LmNet2)^2 / (2 * s2.GS))), type = "l", col='orange',lwd=2.5,
       bty = "n", xlab = "", ylab = "")
  axis(side=4, at = pretty(range(df.GS$SNet1)))
  mtext("fractions retained", side=4, line=3)
  text(labels = paste("ms = ",msNet1,"cm"), x = LmNet1, y =
         ((exp(- ((LmNet1 - LmNet1)^2 / (2 * s2.GS))))/1.1),
       col = 'darkgreen')
  text(labels = paste("ms = ",msNet2,"cm"), x = LmNet2, y =
         ((exp(- ((LmNet2 - LmNet2)^2 / (2 * s2.GS))))/1.1),
       col = 'orange')
  plot2 <- recordPlot()


  #all results together
  df.results <- data.frame(
    "SF" = SF.GS,
    "Lm net1" = LmNet1,
    "Lm net2" = LmNet2,
    "stand.dev." = s.GS
  )

  #save all in list
  results.GS <- list()
  results.GS[[1]] <- df.GS
  results.GS[[2]] <- df.results
  results.GS[[3]] <- plot1
  results.GS[[4]] <- plot2

  names(results.GS) <- c("Dataframe","Selectivity_Results",
                         "Regression_plot","Gear_selection_ogives")

  return(results.GS)
}
