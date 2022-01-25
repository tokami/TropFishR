## ----ReaddataLoadLibraries, message=FALSE, include=FALSE, echo=FALSE----------
knitr::opts_chunk$set(echo = TRUE,
                      cache = FALSE,
                      warning = FALSE,
                      eval = TRUE,
                      error = FALSE,
                      warning = FALSE,
                      message = FALSE,
                      include = TRUE,
                      collapse = TRUE,
                      comment = "#>",
                      fig.show = "hold",
                      fig.width=8, fig.height=7)

## ---- echo = TRUE, include = TRUE, eval = FALSE-------------------------------
#  install.packages("TropFishR", repos = "https://cran.rstudio.com/")

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  install.packages("remotes")
#  remotes::install_github("tokami/TropFishR")

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
library(TropFishR)

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
data("synLFQ7")
lfq <- synLFQ7

## ----Figure 1, echo=TRUE, eval=TRUE, fig.cap="Length frequency data visualised in terms of (a) catches and (b) restructured data with MA = 7."----
## set seed value for reproducible results
set.seed(1)

## adjust bin size
lfq_bin2 <- lfqModify(lfq, bin_size = 2)

## plot raw and restructured LFQ data
ma <- 7
lfq_bin2_res <- lfqRestructure(lfq_bin2, MA = 7, addl.sqrt = FALSE)

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq_bin2_res, Fname = "catch", date.axis = "modern")
plot(lfq_bin2_res, Fname = "rcounts", date.axis = "modern")
par(opar)

## ----elefan_search_space, echo=TRUE, eval=TRUE--------------------------------
## coarse estimate of Linf
linf_guess <- max(lfq_bin2$midLengths) / 0.95

## lower search space bounds
low_par <- list(Linf = 0.8 * linf_guess,
                K = 0.01,
                t_anchor = 0,
                C = 0,
                ts = 0)

## upper search space bounds
up_par <- list(Linf = 1.2 * linf_guess,
               K = 1,
               t_anchor = 1,
               C = 1,
               ts = 1)

## ----Figure 4, fig.height=5, fig.width=5, echo=TRUE, eval=TRUE, results="hide", fig.cap="Score graph of the ELEFAN method with simulated annealing. Green dots indicate the runnning minimum value of the cost function, while blue dots indicate the mean score of each iteration. The red line shows the decline of the 'temperature' value, which describes the probability of accepting worse solutions as the parameter space is explored."----
## run ELEFAN with simulated annealing
res_SA <- ELEFAN_SA(lfq_bin2, SA_time = 60*0.5, SA_temp = 6e5,
                   MA = ma, seasonalised = TRUE, addl.sqrt = FALSE,
                   init_par = list(Linf = linf_guess,
                                   K = 0.5,
                                   t_anchor = 0.5,
                                   C=0.5,
                                   ts = 0.5),
                   low_par = low_par,
                   up_par = up_par)

## show results
res_SA$par
res_SA$Rn_max

## ----Figure 5, fig.height=5, fig.width=5, eval=TRUE, results="hide", fig.cap="Score graph of the ELEFAN method with genetic algorithms. Green dots indicate the runnning maximum value of the fitness function, while blue dots indicate the mean score of each iteration."----
## run ELEFAN with genetic algorithm
res_GA <- ELEFAN_GA(lfq_bin2, MA = ma, seasonalised = TRUE,
                    maxiter = 50, addl.sqrt = FALSE,
                    low_par = low_par,
                    up_par = up_par,
                    monitor = FALSE)

## show results
res_GA$par
res_GA$Rn_max

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  ## list for results
#  JK <- vector("list", length(lfq_bin2$dates))
#  
#  ## loop
#  for(i in 1:length(lfq_bin2$dates)){
#    loop_data <- list(dates = lfq_bin2$dates[-i],
#                    midLengths = lfq_bin2$midLengths,
#                    catch = lfq_bin2$catch[,-i])
#    tmp <- ELEFAN_GA(loop_data, MA = ma, seasonalised = TRUE,
#                      maxiter = 50, addl.sqrt = FALSE,
#                      low_par = low_par,
#                      up_par = up_par,
#                      monitor = FALSE, plot = FALSE)
#    JK[[i]] <- unlist(c(tmp$par, list(Rn_max=tmp$Rn_max)))
#  }
#  
#  ## bind list into dataframe
#  JKres <- do.call(cbind, JK)
#  
#  ## mean
#  JKmeans <- apply(as.matrix(JKres), MARGIN = 1, FUN = mean)
#  
#  ## confidence intervals
#  JKconf <- apply(as.matrix(JKres), MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975)))
#  JKconf <- t(JKconf)
#  colnames(JKconf) <- c("lower","upper")
#  
#  ## show results
#  JKconf

## ----Figure 6, echo = TRUE, fig.cap="Graphical fit of estimated and true growth curves plotted through the length frequency data. The growth curves with the true values are displayed in grey, while the blue and green curves represent the curves of ELEFAN_SA and ELEFAN_GA, respectively."----
## plot LFQ and growth curves
plot(lfq_bin2_res, Fname = "rcounts",date.axis = "modern", ylim=c(0,130))
lt <- lfqFitCurves(lfq_bin2, par = list(Linf=123, K=0.2, t_anchor=0.25, C=0.3, ts=0),
                   draw = TRUE, col = "grey", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq_bin2, par = res_SA$par,
                  draw = TRUE, col = "darkblue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(lfq_bin2, par = res_GA$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)

## -----------------------------------------------------------------------------
## assign estimates to the data list
lfq_bin2 <- lfqModify(lfq_bin2, par = res_GA$par)

## ---- echo=TRUE---------------------------------------------------------------
## estimation of M
Ms <- M_empirical(Linf = lfq_bin2$par$Linf, K_l = lfq_bin2$par$K, method = "Then_growth")
lfq_bin2$par$M <- as.numeric(Ms)

## show results
paste("M =", as.numeric(Ms))

## -----------------------------------------------------------------------------

## define plus group as largest length class smaller than Linf
plus_group <- lfq_bin2$midLengths[max(which(lfq_bin2$midLengths < lfq_bin2$par$Linf))]

## summarise catch matrix into vector and add plus group
lfq_catch_vec <- lfqModify(lfq_bin2, vectorise_catch = TRUE, plus_group = plus_group)


## ----Figure 7,echo=TRUE, fig.width=6, fig.height=5, fig.cap="Catch curve with selected points for the regression analysis and in the second panel the selection ogive with age at first capture.", message = FALSE, warning=FALSE----

## run catch curve
res_cc <- catchCurve(lfq_catch_vec, reg_int = c(18,55), calc_ogive = TRUE)

## assign estimates to the data list
lfq_catch_vec$par$Z <- res_cc$Z
lfq_catch_vec$par$FM <- as.numeric(lfq_catch_vec$par$Z - lfq_catch_vec$par$M)


## ---- echo=FALSE, eval=TRUE---------------------------------------------------
paste("Z =",round(lfq_catch_vec$par$Z,2))
paste("FM =",round(lfq_catch_vec$par$FM,2))
paste("L50 =",round(res_cc$L50,2))

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
lfq_catch_vec$par$E <- lfq_catch_vec$par$FM / lfq_catch_vec$par$Z
paste("E =",round(lfq_catch_vec$par$E,2))

## -----------------------------------------------------------------------------
## assign length-weight parameters to the data list
lfq_catch_vec$par$a <- 0.015
lfq_catch_vec$par$b <- 3

## assign maturity parameters
lfq_catch_vec$par$Lmat <- 35
lfq_catch_vec$par$wmat <- 8


## list with selectivity parameters
selectivity_list <- list(selecType = "trawl_ogive",
                         L50 = res_cc$L50, L75 = res_cc$L75)

## ----Figure 9, echo=TRUE, eval=TRUE, fig.cap="Results of the Thompson and Bell model: (a) Curves of yield and biomass per recruit. The black dot represents yield and biomass under current fishing pressure. The yellow and red dashed lines represent fishing mortality for maximum sustainable yield (Fmax) and fishing mortality to fish the stock at 50% of the virgin biomass (F0.5). (b) exploration of impact of different exploitation rates and Lc values on the relative yield per recruit."----
## Thompson and Bell model with changes in F
TB1 <- predict_mod(lfq_catch_vec, type = "ThompBell",
                   FM_change = seq(0,1.5,0.05),
                   stock_size_1 = 1,
                   curr.E = lfq_catch_vec$par$E,
                   s_list = selectivity_list,
                   plot = FALSE, hide.progressbar = TRUE)

## Thompson and Bell model with changes in F and Lc
TB2 <- predict_mod(lfq_catch_vec, type = "ThompBell",
                   FM_change = seq(0,1.5,0.1),
                   Lc_change = seq(25,50,0.1),
                   stock_size_1 = 1,
                   curr.E = lfq_catch_vec$par$E,
                   curr.Lc = res_cc$L50,
                   s_list = selectivity_list,
                   plot = FALSE, hide.progressbar = TRUE)

## plot results
par(mfrow = c(2,1), mar = c(4,5,2,4.5), oma = c(1,0,0,0))
plot(TB1, mark = TRUE)
mtext("(a)", side = 3, at = -0.1, line = 0.6)
plot(TB2, type = "Isopleth", xaxis1 = "FM", mark = TRUE, contour = 6)
mtext("(b)", side = 3, at = -0.1, line = 0.6)

## Biological reference levels
TB1$df_Es

## Current yield and biomass levels
TB1$currents

