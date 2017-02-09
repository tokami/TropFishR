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
#  install.packages("TropFishR")

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
lfqbin <- lfqRestructure(synLFQ7, MA = 11, addl.sqrt = TRUE)
par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfqbin, Fname = "catch", date.axis = "modern")
plot(lfqbin, Fname = "rcounts", date.axis = "modern")

## ----Figure 2, fig.width=6, fig.height=5, echo=TRUE, eval=TRUE, fig.cap="Powell-Wetherall plot to derive an estimate of Linf."----
# Powell Wetherall plot
res_PW <- powell_wetherall(param = synLFQ7,
                           catch_columns = 1:ncol(synLFQ7$catch),
                           reg_int = c(9,29))
# show results
paste("Linf =",round(res_PW$Linf_est), "±", round(res_PW$se_Linf))

## ---- include=TRUE, fig.width=6, fig.height=6, eval = FALSE, echo=TRUE, fig.cap="Example of a K-Scan application for the estimation of K to a corresponding Linf value. Graph shows the score function for different K values."----
#  # ELEFAN with K-Scan
#  res_KScan <- ELEFAN(synLFQ7, Linf_fix = res_PW$Linf_est,
#                      MA=11, hide.progressbar = TRUE, addl.sqrt = TRUE)
#  
#  # show results
#  res_KScan$par; res_KScan$Rn_max

## ----Figure 3, fig.width=8, eval = FALSE,  fig.cap="Results of response surface analysis of the Thumbprint Emperor data. Red colours indicate higher scoring parameter combinations indicative of a better fit."----
#  # Response surface analyss
#  res_RSA <- ELEFAN(synLFQ7, Linf_range = seq(114,134,1), MA = 11,
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
res_SA <- ELEFAN_SA(synLFQ7, SA_time = 60*1, SA_temp = 6e5,
                    MA = 11, seasonalised = TRUE, addl.sqrt = TRUE,
                    init_par = list(Linf = 124, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                    low_par = list(Linf = 119, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                    up_par = list(Linf = 129, K = 1, t_anchor = 1, C = 1, ts = 1))
# show results
res_SA$par; res_SA$Rn_max

