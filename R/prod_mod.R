#' @title Production models: Schaefer and Fox models
#'
#' @description Schaefer and Fox models are examples of so called surplus production
#'    models, and are used to calculate MSY and BMSY and so on.
#'    First sentence. second sentence.
#'
#' @param Y Yield as vector (catch in weight)
#' @param f Fishing effort as vector
#'
#' @examples
#' # load data
#' data(trawl_fishery_Java)
#'
#' # run model
#' output = prod_mod(Y = trawl_fishery_Java$yield,
#'                      f = trawl_fishery_Java$effort)
#' output[[2]]
#'
#' @details MSY
#'
#' @references
#' Fox, W. W. Jr., 1970. An exponential surplus-yield model for optimizing exploited fish
#' populations. \emph{Trans.Am.Fish.Soc.}, 99:80-88
#'
#' Graham, M., 1935. Modern theory of exploiting a fishery and application to North Sea
#' trawling. \emph{J.Cons.CIEM}, 10(3):264-274
#'
#' Schaefer, M., 1954. Some aspects of the dynamics of populations important to the
#' management of the commercial marine fisheries. \emph{Bull.I-ATTC/Bol. CIAT}, 1(2):27-56
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

prod_mod <- function(Y, f){

  mean.f <- mean(f,na.rm=T)
  sd.f <- sd(f,na.rm=T)
  tn2 <- 2.37   #WHERE IS THIS VALUE FROM????

  #Schaefer
  CPUE.S <- Y/f
  mean.CPUE.S <- mean(CPUE.S,na.rm=T)
  sd.CPUE.S <- sd(CPUE.S,na.rm=T)
  a.S <- summary(lm(CPUE.S ~ f))$coefficients[1]
  b.S <- summary(lm(CPUE.S ~ f))$coefficients[2]
  sb2.S <- ((sd.CPUE.S/sd.f)^2 - b.S^2 )/ (length(f)-2)
  sb.S <- summary(lm(CPUE.S ~ f))$coefficients[4]
  conf.b.S <- c(b.S - tn2 * sb.S, b.S + tn2 * sb.S)
  sa2.S <- sb2.S * (sd.f^2 * (length(f)-1)/length(f) + mean.f^2)
  sa.S <- summary(lm(CPUE.S ~ f))$coefficients[3]
  conf.a.S <- c(a.S - tn2 * sa.S, a.S + tn2 * sa.S)
  MSY.S <- -0.25 * a.S^2/b.S
  fMSY.S <- -0.5 * a.S/b.S

  #Fox
  CPUE.F <- log(Y/f)
  mean.CPUE.F <- mean(CPUE.F,na.rm=T)
  sd.CPUE.F <- sd(CPUE.F,na.rm=T)
  a.F <- summary(lm(CPUE.F ~ f))$coefficients[1]
  b.F <- summary(lm(CPUE.F ~ f))$coefficients[2]
  sb2.F <- ((sd.CPUE.F/sd.f)^2 - b.F^2 )/ (length(f)-2)
  sb.F <- summary(lm(CPUE.F ~ f))$coefficients[4]
  conf.b.F <- c(b.F - tn2 * sb.F, b.F + tn2 * sb.F)
  sa2.F <- sb2.F * (sd.f^2 * (length(f)-1)/length(f) + mean.f^2)
  sa.F <- summary(lm(CPUE.F ~ f))$coefficients[3]
  conf.a.F <- c(a.F - tn2 * sa.F, a.F + tn2 * sa.F)
  MSY.F <- -(1/b.F) * exp(a.F-1)
  fMSY.F <- -1/b.F


  #Plot
  x = seq(min(f),max(f),1)
  y = exp(a.F + b.F*x)
  par(new=F, mar=c(5, 4, 4, 2) + 0.1)
  plot(CPUE.S ~ f, xlab='Fishing effort', ylab='CPUE')
  abline(a = a.S, b=b.S)
  lines(x, y = y, col = 'blue', lty = 5)
  legend("topright", col=c('black','blue'), bty ='n',
         legend = c('Schaefer','Fox'), lty=c(1,5))
  plot1 <- recordPlot()


  x = seq(0,abs(a.S / b.S),1)
  y.S = a.S*x + b.S*x^2
  y.F = x * exp(a.F + b.F*x)
  par(new=F)
  plot(x,rep(max(Y),length(x)),type='n',
       ylim=c(0,max(Y)), xlab = "Fishing effort", ylab = "Yield")
  points(Y ~ f)
  lines(x, y = y.S)
  lines(x, y = y.F,col='blue', lty = 5)
  segments(x0=fMSY.S,x1=fMSY.S,y0=0,y1=MSY.S,lty=2)
  segments(x0=0,x1=fMSY.S,y0=MSY.S,y1=MSY.S,lty=2)
  points(fMSY.S,MSY.S,col='red',pch=16)
  segments(x0=fMSY.F,x1=fMSY.F,y0=0,y1=MSY.F,lty=2,col='blue')
  segments(x0=0,x1=fMSY.F,y0=MSY.F,y1=MSY.F,lty=2,col='blue')
  points(fMSY.F,MSY.F,col='red',pch=16)
  legend("bottomleft", col=c('black','blue'), bty ='n',
         legend = c('Schaefer','Fox'), lty=c(1,5))
  plot2 <- recordPlot()


  results.list <- list()
  results <- data.frame(Schaefer = c(MSY.S,fMSY.S),
                        Fox = c(MSY.F,fMSY.F))

  rownames(results) <- c("MSY","fMSY")
  results.list[[1]] <- results
  results.list[[2]] <- plot1
  results.list[[3]] <- plot2

  return(results.list)
}

