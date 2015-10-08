#' @title Schaefer and Fox models
#'
#' @description Schaefer and Fox models are examples of so called surplus production models, and are used to calculate MSY and BMSY and so on.
#'
#' @param Y Yield as vector (catch in weight)
#' @param f Fishing effort as vector
#'
#' @examples
#' data(ex.SchaeferFox)
#' SchaeferFox(Y = ex.SchaeferFox$yield, f = ex.SchaeferFox$effort)
#'
#' @details MSY
#'
#' @references
#' Graham 1935 Schaefer 1954 Fox 1970
#'
#' @export

SchaeferFox <- function(Y, f){

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

  results <- data.frame(Schaefer = c(MSY.S,fMSY.S),
                        Fox = c(MSY.F,fMSY.F))

  rownames(results) <- c("MSY","fMSY")

  return(results)
}

