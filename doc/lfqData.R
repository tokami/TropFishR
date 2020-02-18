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

## ----eval=FALSE, echo = TRUE--------------------------------------------------
#  lfq1 <- read.csv2("lfq1.csv")
#  
#  lfq1$date <- as.Date(lfq1$date, format = "%d.%m.%Y")
#  
#  lfq1new <- lfqCreate(data = lfq1, Lname = "length", Dname = "date")
#  
#  plot(lfq1new, Fname = "catch")

## ----eval=FALSE, echo = TRUE--------------------------------------------------
#  lfq2 <- read.csv2("lfq2.csv")
#  
#  lfq2$date <- as.Date(lfq2$date, format = "%d/%m/%y")
#  
#  lfq2new <- lfqCreate(data = lfq2, Lname = "Length_CM", Dname = "DATE", Fname = "Frequency")
#  
#  plot(lfq2new, Fname = "catch")

## ----eval=FALSE, echo = TRUE--------------------------------------------------
#  lfq3 <- read.csv2("lfq3.csv")
#  
#  dates <- colnames(lfq3)[-1]
#  dates <- strsplit(dates, "X")
#  dates <- unlist(lapply(dates, function(x) x[2]))
#  dates <- as.Date(dates, "%Y-%d-%m")
#  
#  lfq3new <- list(dates = dates,
#                  midLengths = lfq3$lengthClass,
#                  catch = as.matrix(lfq3[,-1]))
#  class(lfq3new) <- "lfq"
#  lfq3new$catch[is.na(lfq3new$catch)] <- 0
#  
#  plot(lfq3new, Fname = "catch")

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  lfq2a <- read.table("lfq2.txt", sep="\t", dec=',', fileEncoding="latin1", skipNul=TRUE)
#  
#  install.packages("openxlsx")
#  require(openxlsx)
#  lfq3a <- read.xlsx("lfq3.xlsx", sheet = 1)
#  

