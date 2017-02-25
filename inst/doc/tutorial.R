## ----ReaddataLoadLibraries, message=FALSE, include=FALSE, echo=FALSE-----
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

## ---- echo = TRUE, include = TRUE, eval = FALSE--------------------------
#  install.packages("TropFishR", repos = "https://cran.rstudio.com/")

## ---- echo = TRUE, include = TRUE, eval = FALSE--------------------------
#  install.packages("devtools")
#  devtools::install_github("tokami/TropFishR")

## ---- eval=TRUE,echo=TRUE------------------------------------------------
library(TropFishR)

## ------------------------------------------------------------------------
data("synLFQ7")

## ----Figure 1, echo=TRUE, eval=TRUE, fig.cap="Length frequency data visualised in terms of (a) catches and (b) restructured data with MA = 7."----
# set seed value for reproducible results
set.seed(1)

# adjust bin size
synLFQ7 <- lfqModify(synLFQ7, bin_size = 4)

# plot raw and restructured LFQ data
lfqbin <- lfqRestructure(synLFQ7, MA = 5, addl.sqrt = TRUE)
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfqbin, Fname = "catch", date.axis = "modern")
plot(lfqbin, Fname = "rcounts", date.axis = "modern")
par(opar)

## ----Figure 2, fig.width=6, fig.height=5, echo=TRUE, eval=TRUE, fig.cap="Powell-Wetherall plot to derive an estimate of Linf."----
# Powell Wetherall plot
res_PW <- powell_wetherall(param = synLFQ7,
                           catch_columns = 1:ncol(synLFQ7$catch),
                           reg_int = c(10,30))
# show results
paste("Linf =",round(res_PW$Linf_est), "±", round(res_PW$se_Linf))

## ---- include=TRUE, fig.width=6, fig.height=6, eval = FALSE, echo=TRUE, fig.cap="Example of a K-Scan application for the estimation of K to a corresponding Linf value. Graph shows the score function for different K values."----
#  # ELEFAN with K-Scan
#  res_KScan <- ELEFAN(synLFQ7, Linf_fix = res_PW$Linf_est,
#                      MA=5, hide.progressbar = TRUE, addl.sqrt = TRUE)
#  
#  # show results
#  res_KScan$par; res_KScan$Rn_max

## ----Figure 3, fig.width=8, eval = FALSE,  fig.cap="Results of response surface analysis of the Thumbprint Emperor data. Red colours indicate higher scoring parameter combinations indicative of a better fit."----
#  # Response surface analyss
#  res_RSA <- ELEFAN(synLFQ7, Linf_range = seq(119,139,1), MA = 5,
#                    K_range = seq(0.01,2,0.1), addl.sqrt = TRUE,
#                    hide.progressbar = TRUE, contour=5)
#  
#  # show results
#  res_RSA$par; res_RSA$Rn_max

## ---- eval = FALSE, echo=TRUE, include=TRUE, fig.cap="Analysing three highest local maxima with finer resolution. ResSurAna4 produced highest score and best visual fit through restructured data and was therefore chosen as the best solution of the traditional ELEFAN method, growth parameter are displayed in Table 2, column ELEFAN."----
#  # find 3 highest score values
#  n <- length(res_RSA$score_mat)
#  best_scores <- sort(res_RSA$score_mat,partial=n-0:2)[n-0:2]
#  ind <- arrayInd(which(res_RSA$score_mat %in% best_scores),
#                  dim(res_RSA$score_mat))
#  Ks <- as.numeric(rownames(res_RSA$score_mat)[ind[,1]])
#  Linfs <- as.numeric(colnames(res_RSA$score_mat)[ind[,2]])
#  
#  res_loop <- vector("list", 3)
#  for(i in 1:3){
#    tmp <- ELEFAN(synLFQ7,
#                  Linf_range = seq(Linfs[i]-2, Linfs[i]+2, 0.2),
#                  K_range = seq(Ks[i]-0.1, Ks[i]+0.1, 0.05),
#                  MA = 11,
#                  addl.sqrt = TRUE,
#                  hide.progressbar = TRUE,
#                  contour=5)
#    res_loop[[i]] <- cbind(Rn_max=tmp$Rn_max, t(as.matrix(tmp$par)))
#  }
#  results <- do.call(rbind, res_loop)

## ----Figure 4,  fig.height=5, fig.width=5, eval=TRUE, fig.cap="Score graph of the ELEFAN method with simulated annealing. Green dots indicate the runnning minimum value of the cost function, while blue dots indicate the mean score of each iteration. The red line shows the decline of the 'temperature' value, which describes the probability of accepting worse solutions as the parameter space is explored."----
# run ELEFAN with simulated annealing
res_SA <- ELEFAN_SA(synLFQ7, SA_time = 60*0.5, SA_temp = 6e5,
                    MA = 5, seasonalised = TRUE, addl.sqrt = TRUE,
                    init_par = list(Linf = 129, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                    low_par = list(Linf = 119, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                    up_par = list(Linf = 139, K = 1, t_anchor = 1, C = 1, ts = 1))
# show results
res_SA$par; res_SA$Rn_max

## ---- eval = FALSE, echo = TRUE------------------------------------------
#  JK <- vector("list", length(synLFQ7$dates))
#  for(i in 1:length(synLFQ7$dates)){
#    loop_data <- list(dates = synLFQ7$dates[-i],
#                    midLengths = synLFQ7$midLengths,
#                    catch = synLFQ7$catch[,-i])
#    tmp <- ELEFAN_SA(loop_data, SA_time = 60*0.5, SA_temp = 6e5,
#                     MA = 5, addl.sqrt = TRUE,
#                     init_par = list(Linf = 129, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
#                     low_par = list(Linf = 119, K = 0.01, t_anchor = 0, C = 0, ts = 0),
#                     up_par = list(Linf = 139, K = 1, t_anchor = 1, C = 1, ts = 1),
#                     plot = FALSE)
#    JK[[i]] <- unlist(c(tmp$par,list(Rn_max=tmp$Rn_max)))
#  }
#  JKres <- do.call(cbind, JK)
#  # mean
#  JKmeans <- apply(as.matrix(JKres), MARGIN = 1, FUN = mean)
#  # confidence intervals
#  JKconf <- apply(as.matrix(JKres), MARGIN = 1, FUN = function(x) t.test(x)$conf.int[c(1,2)])
#  JKconf <- t(JKconf)
#  colnames(JKconf) <- c("lower","upper")
#  
#  # show results
#  JKconf

## ----Figure 5, fig.height=5, fig.width=5, eval=TRUE, fig.cap="Score graph of the ELEFAN method with genetic algorithm. Green dots indicate the runnning maximum value of the fitness function, while blue dots indicate the mean score of each iteration."----
# run ELEFAN with genetic algorithm
res_GA <- ELEFAN_GA(synLFQ7, MA = 5, seasonalised = TRUE, maxiter = 10, addl.sqrt = TRUE,
                    low_par = list(Linf = 119, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                    up_par = list(Linf = 129, K = 1, t_anchor = 1, C = 1, ts = 1),
                    monitor = FALSE)
# show results
res_GA$par; res_GA$Rn_max

## ----Figure 6, echo = TRUE, fig.cap="Graphical fit of estimated and true growth curves plotted through the length frequency data. The growth curves with the true values are displayed in grey, while the yellow and orange curves represent the curves of the traditional ELEFAN method and ELEFAN with simulated annealing, respectively."----
# plot LFQ and growth curves
plot(lfqbin, Fname = "rcounts",date.axis = "modern", ylim=c(0,130))
lt <- lfqFitCurves(synLFQ7, par = list(Linf=123, K=0.2, t_anchor=0.25, C=0.3, ts=0),
                   draw = TRUE, col = "blue", lty = 1, lwd=1.5)
# lt <- lfqFitCurves(synLFQ7, par = res_RSA$par,
#                    draw = TRUE, col = "goldenrod1", lty = 1, lwd=1.5)
lt <- lfqFitCurves(synLFQ7, par = res_SA$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)
lt <- lfqFitCurves(synLFQ7, par = res_GA$par,
                   draw = TRUE, col = "darkred", lty = 1, lwd=1.5)

## ------------------------------------------------------------------------
# assign estimates to the data list
synLFQ7 <- c(synLFQ7, res_GA$par)

## ---- echo=TRUE----------------------------------------------------------
# estimation of M
Ms <- M_empirical(Linf = res_GA$par$Linf, K_l = res_GA$par$K, method = "Then_growth")
synLFQ7$M <- as.numeric(Ms)
# show results
paste("M =", as.numeric(Ms))

## ----Figure 7,echo=TRUE, fig.width=6, fig.height=5, fig.cap="Catch curve with selected points for the regression analysis and in the second panel the selection ogive with age at first capture.", message = FALSE, warning=FALSE----
# modify the data list
synLFQ7 <- lfqModify(synLFQ7)
# run catch curve
res_cc <- catchCurve(synLFQ7,reg_int = c(14,50), calc_ogive = TRUE)
# assign estimates to the data list
synLFQ7$Z <- res_cc$Z
synLFQ7$FM <- as.numeric(synLFQ7$Z - synLFQ7$M)
synLFQ7$E <- synLFQ7$FM/synLFQ7$Z
paste("Z =",round(synLFQ7$Z,2))
paste("FM =",round(synLFQ7$FM,2))
paste("E =",round(synLFQ7$E,2))
paste("L50 =",round(res_cc$L50,2))

## ----Figure 8, echo=TRUE, fig.cap="Results of Jones' cohort analysis (CA).", message=FALSE,warning=FALSE----
# add plus group which is smaller than Linf
synLFQ7 <- lfqModify(synLFQ7, plus_group = c(TRUE,122))

# assign length-weight parameters to the data list
synLFQ7$a <- 0.01
synLFQ7$b <- 3
# run CA
vpa_res <- VPA(param = synLFQ7, terminalF = synLFQ7$FM,
               analysis_type = "CA",
               plot=TRUE, catch_corFac = (1+4/12))
# stock size
sum(vpa_res$annualMeanNr, na.rm =TRUE) / 1e3
# stock biomass
sum(vpa_res$meanBiomassTon, na.rm = TRUE)
# assign F per length class to the data list
synLFQ7$FM <- vpa_res$FM_calc

## ----Figure 9, echo=TRUE, eval=TRUE, fig.cap="Thompson and Bell model for the Thumbprint Emperor data: (a) Curves of yield and biomass per recruit. The black dot represents yield and biomass under current fishing pressure. The yellow and red dashed lines represent fishing mortality for maximum sustainable yield (Fmsy) and fishing mortality to fish the stock at 50% of the virgin biomass (F0.5). (b) exploration of impact of different exploitation rates and Lc values on the relative yield per recruit."----
# Thompson and Bell model with changes in F
TB1 <- predict_mod(synLFQ7, type = "ThompBell",
                   FM_change = seq(0,1.5,0.05),  stock_size_1 = 1,
                   curr.E = synLFQ7$E, plot = FALSE, hide.progressbar = TRUE)
# Thompson and Bell model with changes in F and Lc
TB2 <- predict_mod(synLFQ7, type = "ThompBell",
                   FM_change = seq(0,1.5,0.1), Lc_change = seq(20,50,0.1),
                   stock_size_1 = 1,
                   curr.E = synLFQ7$E, curr.Lc = res_cc$L50,
                   s_list = list(selecType = "trawl_ogive",
                                 L50 = res_cc$L50, L75 = res_cc$L75),
                   plot = FALSE, hide.progressbar = TRUE)
# plot results
par(mfrow = c(2,1), mar = c(4,5,2,4.5), oma = c(1,0,0,0))
plot(TB1, mark = TRUE)
mtext("(a)", side = 3, at = -1, line = 0.6)
plot(TB2, type = "Isopleth", xaxis1 = "FM", mark = TRUE, contour = 6)
mtext("(b)", side = 3, at = -0.1, line = 0.6)
# Biological reference levels
TB1$df_Es
# Current yield and biomass levels
TB1$currents

