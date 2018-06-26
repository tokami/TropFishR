## ----ReaddataLoadLibraries, message=FALSE, include=FALSE, echo=FALSE-----
knitr::opts_chunk$set(echo = TRUE,
                      cache = FALSE,
                      warning = FALSE,
                      eval = TRUE,
                      error = FALSE,
                      message = FALSE,
                      include = TRUE,
                      collapse = TRUE,
                      comment = "#>",
                      fig.show = "hold",
                      fig.width=8, fig.height=6,
                      fig.align="center")

## ---- echo = TRUE, include = TRUE, eval = FALSE--------------------------
#  install.packages("TropFishR", repos = "https://cran.rstudio.com/")

## ----plot, eval=TRUE, echo=TRUE------------------------------------------
library(TropFishR)

## ---- eval=TRUE,echo=TRUE, fig.width=7, fig.height=4---------------------
## load data set into R environment
data("synLFQ8")
## plot length-frequency data
plot(synLFQ8, Fname = "catch")

## ---- eval=TRUE,echo=TRUE------------------------------------------------
## aggregate lfq data per year
lfq <- lfqModify(synLFQ8, aggregate = "year")

## ---- eval=TRUE,echo=TRUE, fig.width=7, fig.height=4---------------------
## plot lfq data in standard TropFishR manner
plot(lfq, Fname = "catch")

## ---- eval=TRUE,echo=TRUE, fig.width=7, fig.height=6---------------------
## plot data in LBB manner
plotLBB.data(lfq)

## ---- eval=TRUE,echo=TRUE------------------------------------------------
## add length at maturity to lfq data
lfq$Lm50 <- 38  

## ---- eval=FALSE,echo=TRUE-----------------------------------------------
#  ## run LBB model
#  res <- LBB(lfq, plot = TRUE)

## ---- eval=TRUE,echo=FALSE, fig.width=7, fig.height=6, tidy=FALSE, size="\\tiny", out.width="0.8\\linewidth"----
if(require(rjags)){
    res <- LBB(lfq, n.cluster = 1,
               plot = TRUE)
}else{
    data("synLFQ8res")
    res <- synLFQ8res
}

## ---- eval=TRUE,echo=TRUE, fig.width=4, fig.height=3.2-------------------
par(cex=0.7)
## plot results
plotLBB(res)

## ---- eval=TRUE,echo=TRUE, fig.width=8, fig.height=3---------------------
## plot results as time series
plotLBB.ts(res)

