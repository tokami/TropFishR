## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  fig.width=6, fig.height=5
)

## ------------------------------------------------------------------------
library(TropFishR)
t <- seq(-0.2, 3, length.out = 200)
K <- c(2, 1, 0.5)
COL <- rep(1,3)
LTY <- 1:3
for(i in seq(K)){
  Lt <- VBGF(param = list(Linf = 20, K = K[i], t0 = -0.1), t = t)
  if(i == 1){
    plot(t, Lt, t="l", ylim = c(0,22), yaxs="i", col = COL[i], lty = LTY[i])
    abline(v = 0, col = 8, lty = 3)
    abline(h = 20, col = 3, lty = 3)
    points(x = -0.1, y = 0, pch = 16, col = 4)
    text(x = -0.1, y = 0, labels = expression(italic(t[0])), adj=c(1,-0.5), col=4)
    text(x = -0.1, y = 20, labels = expression(italic(L[infinity])), adj=c(1,-0.5), col=3)
    legend("bottomright", legend = paste("K =", K), lty=LTY, col=COL, bty="n")
  }else{
    lines(t, Lt, col = COL[i], lty = LTY[i])
  }
}

## ------------------------------------------------------------------------
library(TropFishR)
t <- seq(-0.2, 3, length.out = 200)
Lt <- VBGF(param = list(Linf = 20, K = 1, t0 = -0.1, ts = 0, C=0), t = t)
Cs <- seq(0.25,1,0.25)
COLs <- 1:5
plot(t, Lt, t="l", ylim=c(0,22), yaxs="i", col=1)
for(i in seq(Cs)){
  lines(t, VBGF(param = list(Linf = 20, K = 1, t0 = -0.1, ts = 0, C=Cs[i]), t = t), col=COLs[i+1])
}
legend("bottomright", legend=paste("C =", c(0,Cs)), lty=1, col=COLs, bty="n")
abline(v=0, col = 8, lty = 3)
abline(h = 20, col = 8, lty=3)

## ------------------------------------------------------------------------
data("alba")
tmplfq <- list(
  midLengths = alba$midLengths,
  dates = alba$dates,
  catch = alba$catch
)

## ------------------------------------------------------------------------
class(tmplfq) <- "lfq"
plot(tmplfq, Fname="catch", hist.sc = 1)

## ------------------------------------------------------------------------
alba <- lfqRestructure(alba, MA = 7)
plot(alba, Fname="rcounts", hist.sc = 0.75)

## ------------------------------------------------------------------------
alba <- lfqRestructure(alba, MA=7)
plot(alba, Fname="rcounts", hist.col = c("white", "black"), 
  image.col = c(rep(rgb(1,0.8,0.8),1000), "white", rep(rgb(0.8,0.8,1),1000)),
  ylim = c(0,max(alba$midLengths+0.5)))
tmp <- lfqFitCurves(alba, par = list(Linf=11, K=2.5, t_anchor=0.5), 
  draw = TRUE, col=4, lty=2)

## ------------------------------------------------------------------------
PW <- powell_wetherall(alba, catch_columns = 1:7, reg_int = c(2,9) )
PW$Linf_est
PW$confidenceInt_Linf

## ------------------------------------------------------------------------
alba2 <- ELEFAN(
  x = alba,  MA = 7,
  Linf_range = seq(7, 20, length.out = 30),
  K_range = exp(seq(log(0.1),log(5), length.out = 30)),
  method = "cross",
  cross.date = alba$dates[3], 
  cross.midLength = alba$midLengths[5], 
  contour = TRUE, add.values = FALSE,
  hide.progressbar = TRUE # change to 'TRUE' to follow algorithm's progression
)
points(alba2$par["Linf"], alba2$par["K"], pch="*", cex=2, col=2)
unlist(alba2$par)
alba2$Rn_max

## ------------------------------------------------------------------------
plot(alba2, Fname="rcounts")
points(alba$dates[3], alba$midLengths[5], pch="*", cex=2, col=2)

## ------------------------------------------------------------------------
set.seed(1)
alba3 <- ELEFAN_SA(
  x = alba, 
  seasonalised = FALSE, 
  init_par = alba2$par[1:5], 
  low_par = list(Linf=PW$confidenceInt_Linf[1], K=1, t_anchor=0, ts=0, C=0),
  up_par = list(Linf=PW$confidenceInt_Linf[2], K=4, t_anchor=1, ts=1, C=1), 
  SA_temp = 2e5,
  SA_time = 60,
  maxit = 400,
  MA = 7,
  plot.score = TRUE, 
  verbose = FALSE
)
unlist(alba3$par)
alba3$Rn_max

## ------------------------------------------------------------------------
plot(alba3, Fname="rcounts")

## ------------------------------------------------------------------------
set.seed(1)
alba4 <- ELEFAN_GA(
  x = alba, 
  seasonalised = FALSE, 
  low_par = list(Linf=PW$confidenceInt_Linf[1], K=1, t_anchor=0, ts=0, C=0),
  up_par = list(Linf=PW$confidenceInt_Linf[2], K=4, t_anchor=1, ts=1, C=1), 
  popSize = 60,
  pmutation = 0.2,
  maxiter = 100,
  run = 20,
  MA = 7,
  plot.score = TRUE, 
  monitor = FALSE,
  parallel = FALSE
)
unlist(alba4$par)
alba4$Rn_max

## ------------------------------------------------------------------------
plot(alba4, Fname="rcounts")

## ------------------------------------------------------------------------
set.seed(1)
alba5 <- ELEFAN_GA(
  x = alba, 
  seasonalised = TRUE, 
  low_par = list(Linf=PW$confidenceInt_Linf[1], K=0.1, t_anchor=0, ts=0, C=0),
  up_par = list(Linf=PW$confidenceInt_Linf[2], K=4, t_anchor=1, ts=1, C=1), 
  popSize = 60,
  pmutation = 0.2,
  maxiter = 100,
  run = 20,
  MA = 7,
  plot.score = TRUE, 
  monitor = FALSE,
  parallel = FALSE
)
unlist(alba5$par)
alba5$Rn_max
plot(alba5, Fname="rcounts")

## ------------------------------------------------------------------------
obj <- alba5
time <- date2yeardec(obj$dates)
t <- seq(floor(min(time)), ceiling(max(time)),0.01)
PAR <- obj$par
PAR$t_anchor <- PAR$t_anchor + floor(min(time))

Lt <- VBGF(param = PAR, t = t)
Kt <- (log(PAR$Linf-Lt[-length(t)]) - log(PAR$Linf-Lt[-1]))/0.01
plot((t[-1]+t[-length(t)])/2, Kt, t="l", ylim=range(Kt)+(c(-1,1)*dist(range(Kt)*0.1)), xlab = "", ylab=bquote(K[t]))
abline(h=PAR$K, lty=2, col=8)
tmin <- (PAR$ts+1.5)%%1
tmax <- (PAR$ts+1)%%1
points(x = floor(min(t))+1+c(tmin,tmax), y = range(Kt))
abline(v = floor(min(t))+1+c(tmin,tmax), lty=2, col=c(2,3))
text(x = floor(min(t))+1 + tmax, y = max(Kt), labels = bquote(t[max] == .(sprintf(fmt = "%.2f",round(tmax,2)))~.(paste0("(",format(yeardec2date(tmax), "%b-%d"),")"))), pos=3, srt=90, offset = -5)
text(x = floor(min(t))+1 + tmin, y = min(Kt), labels = bquote(t[min] == .(sprintf(fmt = "%.2f",round(tmin,2)))~.(paste0("(",format(yeardec2date(tmin), "%b-%d"),")"))), pos=3, srt=90, offset = 05)
mtext("Instantaneous growth rate (K)", side=3, line=0.25)

## ------------------------------------------------------------------------
true_par <- list(Linf = 80, K = 0.5, t_anchor = 0.25,C = 0.75,  ts = 0.5, phiL = 3.51)

## ------------------------------------------------------------------------
set.seed(1)
synLFQ4 <- ELEFAN_GA(
  x = synLFQ4, 
  seasonalised = TRUE, 
  low_par = list(Linf=70, K=0.1, t_anchor=0, ts=0, C=0),
  up_par = list(Linf=110, K=1, t_anchor=1, ts=1, C=1), 
  popSize = 60,
  pmutation = 0.2,
  maxiter = 100,
  run = 20,
  MA = 11,
  plot.score = TRUE, 
  monitor = FALSE,
  parallel = FALSE
)

## ------------------------------------------------------------------------
tmp <- as.data.frame(rbind(unlist(true_par), unlist(synLFQ4$par)))
rownames(tmp) <- c("true", "estimated")
tmp$Rn <- c(lfqFitCurves(synLFQ4, par = true_par)$Rn_max, synLFQ4$Rn_max)
tmp <- round(tmp,3)
tmp

## ------------------------------------------------------------------------
plot(synLFQ4, Fname="rcounts", draw = FALSE, 
  image.col = NA, region.col = "grey90") 
tmp <- lfqFitCurves(synLFQ4, par = true_par, col=5,
  lty=1, draw = TRUE)
tmp <- lfqFitCurves(synLFQ4, par = synLFQ4$par, col=4,
  lty=2, draw = TRUE)
legend("top", ncol=2, legend = c("true", "estimated"),
  col=c(5,4), lty=c(1,2), bg = "grey90", cex = 0.75)

## ------------------------------------------------------------------------
synLFQ4 <- ELEFAN_GA(
  x = synLFQ4, 
  seasonalised = TRUE, 
  low_par = list(Linf=70, K=0.1, t_anchor=0, ts=0, C=0),
  up_par = list(Linf=110, K=1, t_anchor=1, ts=1, C=1), 
  popSize = 60,
  pmutation = function(...) GA::ga_pmutation(..., p0=0.5, p=0.1),
  maxiter = 100,
  run = 20,
  MA = 11,
  plot.score = TRUE, 
  monitor = FALSE,
  parallel = FALSE
)

