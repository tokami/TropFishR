#' @title Non-equilibrium dynamic biomass models
#'
#' @description Non-equilibrium dynamic biomass models
#'
#' @param param a list of parameters
#'
#' @examples
#' \donttest{
#' }
#'
#' @details MSY
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export


prod_mod_dyn <- function(param){

  # Y= Yield(catch) in tons I= abundance index (cpue)kg
  #SAYA Catch and Effort data 1989-2004 DYNAMIC BIOMASS MODEL
  yrs <- 1989:2004
  Y <- c( 2177, 1410, 1782, 2825, 3173, 3142, 2957, 2283, 1798, 2054, 2107, 2099, 1283, 2090, 2354, 1689) I <- c( 74.3, 73.0, 88.0, 67.7, 69.1, 66.6, 67.0, 57.8, 71.8, 75.9, 70.1, 77.8, 124.1, 83.3, 80.1, 71.2)
  #NAZARETH catch and effort data DYNAMIC BIOMASS MODEL 1989-2004
  yrs <- 1989:2004
  Y <- c( 837,914,793,952,1358,1494,1533,1253,1720,1086,1121,1080,1366,918,468,855) I <- c( 75.7,78.5,81.3,78.2,66.1,66.5,64.3,52.5,66.1,81.4,76.2,90.5,99.2,93.3,72.8,84.2)
  # FORMULA== By+1 = By + r By(1 - By/K) - Yy
  # Yy=yield/catch
  # Initial parameters of starting biomass B0, Carrying capacity K and the rate of pop. growth r
  # q= coefficient catchability
  # There are four parameters to be estimated K, B0,r and q
  #It is not advisable to estimate all- first estimate K, B0 and r with fixed q, then estimate r and q, then all
  B0<-2*mean(Y) K<-B0*1.3
  r<-1 q<-mean(I)/B0
  # biomass greater than catch # K>B0- carrying capacity
  # rate of population growth # as I=qB
  input<-c(K,B0,r) B<-B0 par(mfrow=c(3,3))
  ### ssefn<-function(input){
  K<-input[1] B0<-input[2] r<-input[3]
  B<-B0 Yvec<-NULL Bvec<-NULL Ihat<-NULL yrs<-1:length(Y) for (y in yrs){
    SY<-r*B*(1-B/K) Bvec<-c(Bvec,B) Ihat<-c(Ihat,q*B) B<-B+SY-Y[y] B<-ifelse(B<0,0,B)
  } SSE<-sum((I-Ihat)^2) return(SSE)



  #### to optimise estA<-nlm(ssefn,input,typsize=input,iterlim=1000) estA<-nlm(ssefn,estA$est,typsize=input,iterlim=1000) estA
  # nlm- non linear minimization
  # using the result of the first estimate we estimate again
  #### estimate r and q K <- estA$est[1]
  B0 <- estA$est[2]
  r <- estA$est[3]
  q <- mean(I)/B0
  input2 <- c(r,q)
  B <- B0 ssefn2<-function(input){
    r<-input[1] q<-input[2]
    B<-B0 Yvec<-NULL
    # now we estimate the two parameters
    UNU-Fisheries Training Programme 54
    Dharmendra
    Bvec<-NULL Ihat<-NULL yrs<-1:length(Y) for (y in yrs){
      SY<-r*B*(1-B/K) Bvec<-c(Bvec,B) Ihat<-c(Ihat,q*B) B<-B+SY-Y[y] B<-ifelse(B<0,0,B)
    } SSE<-sum((I-Ihat)^2)
    # plot(yrs,I, type="b", xlab="year", ylab="CPUE") # lines(Ihat, col="red")
    return(SSE) }
  estB <- nlm(ssefn2,input2,typsize=input2,iterlim=1000)
  estB <- nlm(ssefn2,estB$est,typsize=estB$est,iterlim=1000)
  estB
  #####
  # plotting these data
  # to calculate the predicted biomass for all ages and the predicted index K<-estA$est[1]
  B0<-estA$est[2]
  r<-estB$est[1]
  q<-estB$est[2]
  B<-B0 Yvec<-NULL Bvec<-NULL Ihat<-NULL yrs<-1:length(Y) for (y in yrs){
    SY<-r*B*(1-B/K) Bvec<-c(Bvec,B) Ihat<-c(Ihat,q*B) B<-B+SY-Y[y]
  }
  par(mfrow=c(3,3))
  plot(yrs, Y, type="b", xlab="year", ylab="yield", main="Yield Trajectory NAZARETH", ylim=c(0, max(Y)*1.05)) plot(yrs, Bvec, type="b", xlab="year", ylab="biomass", main="Biomass Trajectory", ylim=c(0, max(Bvec)*1.05)) plot(Y, Bvec, xlab="yield", ylab="predicted Biomass", main="Corr. between yield and biomass")
  cor(Bvec,Y)
  plot(Bvec, I, xlab="Bvec-biomass", ylab="I-cpue", main="Corr. between biomass and CPUE") lines(lowess(Bvec, I), col=2)
  plot(log(Bvec), log(I), xlab="Bvec-biomass", ylab="I-cpue", main="Corr. between biomass and CPUE") plot(yrs,I, type="b",xlab="Year", ylab="CPUE", ylim=c(0, max(I)*1.05))
  plot(I, Ihat, xlab="observed CPUE", ylab="Ihat", xlim=c(0, max(c(I,Ihat))), ylim=c(0, max(c(I,Ihat)))) plot(yrs, I, xlab="years", ylab="CPUE", ylim=c(0, max(I)*1.05), type="b")
  lines(yrs, Ihat, col=2, type="b")
  cor(Bvec, I)
  ## the equilibrium yield
  Blevels <- seq(0,3000,10)
  EYlevels <- r*Blevels*(1-Blevels/K)
  EYlevels <- ifelse(EYlevels<0, 0, EYlevels)
  plot(Blevels, EYlevels, type="l" , xlab="Biomass", ylab="Yield", main="Yield Curve NAZARETH", ylim=c(0,max(c(EYlevels,Yvec))))
  lines(Bvec,Yvec,type="b", col=2)
  Elevels <- seq(0,10000,1000) plot(Elevels, EYlevels, type="l")
  MSY <- r*K/4
  Bmsy <- K/2
  Emsy <- r/(2*q)*1000
  Fmsy <- r/2
  MSY
  Emsy
  abline(h=MSY, type="b",col=2) abline(v=Bmsy, col=2)




}


