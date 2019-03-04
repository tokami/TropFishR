#' @title Bhattacharya's method
#'
#' @description Find relative frequencies and frequency distribution of distinct cohorts
#'      in the observed length frequency distribution by resolving it into
#'      Gaussian components.
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \code{midLengths} midpoints of the length class as vector,
#'   \item \code{catch} a vector with the catch per length class or a matrix with
#'      catches per length class of subsequent years;
#' }
#' @param n_rnorm number of observations for the function \code{\link{rnorm}}. The
#'    default is 1000.
#' @param savePlots logical; indicating whether the analyis graphs should be recorded
#'
#' @keywords function Bhattacharya length-frequency lfq
#'
#' @examples
#' # The following example requires to choose certain values for the regression analyses:
#' # first cohort:   point  2 to  8
#' # second cohort:  point 12 to 17
#' # third cohort:   point 19 to 23
#' # fourth cohort:  point 26 to 30
#'
#' \donttest{
#'  data(synLFQ1)
#'  Bhattacharya(param = synLFQ1)
#' }
#'
#' @details This method includes the \code{\link{identify}} function, which allows to
#'      choose points on a graphical device manually. To stop this process please press
#'      right mouse click on the graph, and in case you are working on a windows
#'      maschine click on "Stop". An error will be caused if the graphical device
#'      is closed manually. After you have selected the points for regression analysis
#'      you will be asked if you want to redo the selection or if you want to continue.
#'      Please enter in the Console "y" for continuation if you are satisfied
#'      with your selection and the corresponding Gaussian distribution or enter
#'      "redo" if you want to repeat the selection procedure. This function
#'      allows a maximum of 12 cohorts or seperate distributions in one sample. Please
#'      find more details in the Vignette of this package or in the FAO manual by
#'      Sparre and Venema (1998) (\strong{Section 3.4.1, p. 80}).
#'
#' @return A list with the input parameters and
#'    \itemize{
#'    \item \strong{regressionLines} dataframe with intercept, slope, start and end points
#'       of the regression lines,
#'    \item \strong{Lmean_SD_list} dataframe with the mean length (Lmean), standard deviation (SD),
#'      and seperation index (SI) for each cohort,
#'    \item \strong{bhat_results} dataframe with the results of the Bhattacharya method,
#'    \item \strong{distributions} list with the x and y values of selected distributions,
#'    \item \strong{cohort_plots} list with analysis plots (when savePlots = TRUE).
#'    }
#'
#' @importFrom grDevices dev.new dev.off recordPlot
#' @importFrom graphics abline hist identify layout lines par plot points text title
#' @importFrom stats dnorm lm rnorm sd
#' @importFrom utils flush.console
#'
#' @references
#' Bhattacharya, C.G., 1967. A simple method of resolution of a distribution into Gaussian
#' components, \emph{Biometrics}, 23:115-135
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export


Bhattacharya <- function(param, n_rnorm = 1000, savePlots = FALSE){

  res <- param
  if("midLengths" %in% names(res) == TRUE){
    midLengths <- as.character(res$midLengths)
    # create column without plus group (sign) if present
    midLengths.num <- do.call(rbind,strsplit(midLengths, split="\\+"))
    midLengths.num <- as.numeric(midLengths.num[,1])
  }else stop(
    noquote('The parameter list does not contain an object with name "midLengths"'))
  if("catch" %in% names(res) == TRUE){
    catch <- res$catch
  }else stop(
    noquote('The parameter list does not contain an object with name "catch"'))

  #Transform matrix into vector if provided
  if(class(catch) == 'matrix'){
    catch.vec <- rowSums(catch, na.rm = TRUE)
  }else catch.vec <- catch
  if(length(midLengths.num) != length(catch.vec)) stop(
    noquote("midLengths and catch do not have the same length!"))

  #calculate size class interval
  interval <- midLengths.num[2] - midLengths.num[1]

  #STEP 1: fills second column of bhat with counts per size class
  bhat.table <- data.frame('mean.length.classes' = midLengths.num,
                           'N1.plus' = catch.vec,
                           'log.N1.plus' = NA,
                           'delta.log.N1.plus' = NA,
                           'L' = NA,
                           'delta.log.N' = NA,
                           'log.N1' = NA,
                           'N1' = NA,
                           'N2.plus' = NA)

  bhat.table.list <- vector("list", 13)
  a.b.list <- vector("list", 12)
  l.s.list <- vector("list", 12)
  temp.list <- vector("list", 12)
  coho_plot <- vector("list",12)


  colour.xy <- c('blue','darkgreen','red','goldenrod2','purple','orange',
                 'lightgreen','skyblue','brown','darkblue','darkorange','darkred')
  colour.vec <- rep('black', length(bhat.table$L))

  writeLines("Starting on the left, please choose from black points which \nlie on a straight line! Do not include points which might be \naffected by the next distribution! \nTo stop the process press right click (and choose \n'Stop' if necessary)")
  flush.console()

  for(xy in 1:12){

    #STEP 2: fills third column
    bhat.table$log.N1.plus <- round(log(bhat.table$N1.plus),digits=3)

    #STEP 3: fills fourth coulmn
    if(xy > 1) bhat.table$log.N1.plus[last_choices[2]] <- 0
    for(i in 2:length(bhat.table$log.N1.plus)){
      delta.value <- bhat.table$log.N1.plus[i] - bhat.table$log.N1.plus[i-1]
      bhat.table$delta.log.N1.plus[i] <- delta.value
    }

    #STEP 4: fills fifth coulmn
    checki <- !is.na(bhat.table$delta.log.N1.plus) &
      bhat.table$delta.log.N1.plus != 'Inf'
    bhat.table$L[checki] <- bhat.table$mean.length.classes[checki] - (interval/2)

    #STEP 5: plot of fifth against fourth column
    checki2 <- !is.na(bhat.table$delta.log.N1.plus) &
      bhat.table$delta.log.N1.plus != 'Inf'
    bhat.table$L[checki2] <- (bhat.table$mean.length.classes[checki2] - (interval/2))

    #STEP 6: select points for regression line
    if(xy == 1){bhat.table.list[[1]] <- bhat.table}

    repeat {

      if("id.co1" %in% ls()){
        if(length(id.co1$ind) == 0) break
      }

      repeat {
        # histogramm
        dev.new()#noRStudioGD = TRUE)
        par(mfrow = c(2,1),
            oma = c(6.5, 0.5, 2, 1) + 0.1,
            mar = c(0, 4, 0.5, 0) + 0.1)
        layout(matrix(c(1,2),nrow=2), heights = c(1,2.5))

        freqis <- rep(bhat.table.list[[1]]$mean.length.classes, bhat.table.list[[1]]$N1.plus)
        if(xy == 1){
          hist(freqis, breaks = 50, main = '', xaxt = 'n', probability = TRUE)
          mtext(side = 3, "Click on two points. Escape to Quit.",
                      xpd = NA, cex = 1.25)
        }
        if(xy > 1){
          hist(freqis, breaks = 50, main = '', xaxt = 'n',
                        probability = TRUE, ylim = c(0,maxlim))
          mtext(side = 3, "Click on two points. Escape to Quit.",
                      xpd = NA, cex = 1.25)
        }
        if(xy > 1){
          for(tempRep in 1:(xy-1)){
            lines(temp.list[[tempRep]]$xfit, temp.list[[tempRep]]$yyfit,
                  col=colour.xy[tempRep], lwd=1.5)
          }
        }
        # bhatta plot
        if(xy == 1){
          plot(bhat.table.list[[1]]$L, bhat.table.list[[1]]$delta.log.N1.plus,
               col = colour.vec, pch = 16,
               ylab = expression(paste(delta," log N+")))
          abline(h = 0)
          text(bhat.table.list[[1]]$L, bhat.table.list[[1]]$delta.log.N1.plus,
               labels = row.names(bhat.table.list[[1]]), cex= 0.7, pos=3)
          title(xlab = "L", outer = TRUE, line = 2.5)
        }
        if(xy > 1){
          plot(bhat.table.list[[1]]$L, bhat.table.list[[1]]$delta.log.N1.plus, pch = 4,
               col = colour.vec,
               ylab = expression(paste(delta," log N+")))
          abline(h = 0)
          abline(a = a.b.list[[xy-1]][1], b=a.b.list[[xy-1]][2],
                 col=colour.xy[xy-1], lwd = 1.5)
          points(bhat.table$L[(last_choices[2]+1):length(bhat.table$L)],
                 bhat.table$delta.log.N1.plus[(last_choices[2]+1):length(bhat.table$L)],
                 col = 'black', pch = 16)
          text(bhat.table$L[(last_choices[2]+1):length(bhat.table$L)],
               bhat.table$delta.log.N1.plus[(last_choices[2]+1):length(bhat.table$L)],
               labels=row.names(bhat.table.list[[1]][
                 (last_choices[2]+1):length(bhat.table$delta.log.N1.plus),]),
               cex = 0.7, pos = 3)
          title(xlab = "L", outer = TRUE, line = 2.5)
        }
        # identify regression line
        if(xy > 1) rm(id.co1)
        id.co1 <- identify(bhat.table$L, bhat.table$delta.log.N1.plus, n = 2, pos = TRUE)
        if(length(id.co1$ind) == 0) break

        colour.vec[id.co1$ind[1]:id.co1$ind[2]] <- colour.xy[xy]
        dev.off()

        #STEP 7: calculate mean length and standard deviation of regression line
        x.co1 <- bhat.table$L[id.co1$ind[1]:id.co1$ind[2]]

        # Break loop if selection does not embrace at least 3 points
        if(length(x.co1) < 3){
          writeLines("Your selection is not possible. You have to choose two \npoints which include at least one other point. At least \nthree points are required for a regression line. Please choose again!")
          flush.console()
        }
        if(length(x.co1) >= 3){
          break
        }
      }

      if(xy == 1 & length(id.co1$ind) == 0){
        writeLines("You did not choose any points! Please run the function again \nand choose points to include in the calculation of the \nGaussian components of the observed frequency.")
        flush.console()
        cancel <- TRUE
        break
      }
      if(xy > 1 & length(id.co1$ind) == 0){
        cancel <- FALSE
        break
      }

      y.co1 <- bhat.table$delta.log.N1.plus[id.co1$ind[1]:id.co1$ind[2]]
      m.co1 <- lm(y.co1 ~ x.co1)
      sum.m.co1 <- summary(m.co1)
      a.co1 <- sum.m.co1$coefficients[1]
      b.co1 <- sum.m.co1$coefficients[2]
      l.mean.co1 <- -a.co1/b.co1   #mean length: L(mean)(N1) = -a/b
      s.co1 <- sqrt(-1/b.co1)      #standard deviation: s(N1) = sqrt(-1/b)

      a.b.list[[xy]] <- c(a.co1,b.co1)
      l.s.list[[xy]] <- c(l.mean.co1,s.co1)


      #STEP 8: fill sixth column
      #set.seed ????
      normal.dis.co1 <- rnorm(n=n_rnorm, mean = l.mean.co1, sd=s.co1)
      max.class.ind.co1 <- which(round(max(normal.dis.co1), digits = 0) >= (
        bhat.table$mean.length.classes - (interval/2)) &
          round(max(normal.dis.co1),digits = 0) < (
            bhat.table$mean.length.classes + (interval/2)))
      if(length(max.class.ind.co1) == 0) max.class.ind.co1 <- length(bhat.table$mean.length.classes)
      max.class.co1 <- bhat.table$mean.length.classes[max.class.ind.co1]
      for(i in 1:max.class.ind.co1){
        delta.N <- round(a.co1 + b.co1 * bhat.table$L[i],digits=3)
        bhat.table$delta.log.N[i] <- delta.N
      }


      # Draw plot for checking
      if(xy == 1) min.class.ind.co1 <- 1
      if(xy > 1){
        min.class.ind.co1 <- which(round(min(normal.dis.co1), digits = 0) >= (
          bhat.table$mean.length.classes - (interval/2)) &
            round(min(normal.dis.co1),digits = 0) < (
              bhat.table$mean.length.classes + (interval/2)))
        if(length(min.class.ind.co1) == 0) min.class.ind.co1 <- 1
      }
      min.class.co1 <- bhat.table$mean.length.classes[min.class.ind.co1]

      xfit <- seq(min.class.co1, max.class.co1, length = 100)
      yfit <- dnorm(xfit, mean=l.mean.co1, sd=s.co1)

      freqis <- rep(bhat.table.list[[1]]$mean.length.classes,
                    bhat.table.list[[1]]$N1.plus)
      h <- hist(freqis, breaks = 50, plot = FALSE)

      dev.new()#noRStudioGD = TRUE)
      par(mfrow = c(2,1),
          oma = c(6.5, 0.5, 2, 1) + 0.1,
          mar = c(0, 4, 0.5, 0) + 0.1)
      layout(matrix(c(1,2),nrow=2), heights = c(1,2.5))
      # histogram
      yyfit <-  yfit * abs((max(yfit) - max(h$density)) / max(yfit))
      if(xy == 1) maxlim <- ifelse(max(yyfit) > max(h$density), max(yyfit), max(h$density))
      hist(freqis, breaks = 50, main = "", xaxt = 'n',
           probability = TRUE, ylim = c(0,maxlim))
      if(xy > 1){
        for(tempRep in 1:(xy-1)){
          lines(temp.list[[tempRep]]$xfit, temp.list[[tempRep]]$yyfit,
                col=colour.xy[tempRep], lwd=1.5)
        }
      }
      lines(xfit, yyfit, col=colour.xy[xy], lwd=1.5)
      temp.list[[xy]] <- list(xfit = xfit,
                              yyfit = yyfit)
      # bhatta plot
      plot(bhat.table.list[[1]]$L, bhat.table.list[[1]]$delta.log.N1.plus, pch = 4,
           col = colour.vec,
           ylab = expression(paste(delta," log N+")))
      abline(h = 0)
      abline(a = a.co1, b=b.co1, col=colour.xy[xy], lwd = 1.5)
      points(bhat.table$L[(id.co1$ind[2]+1):length(bhat.table$L)],
             bhat.table$delta.log.N1.plus[(id.co1$ind[2]+1):length(bhat.table$L)],
             col = 'black', pch = 16)
      text(bhat.table$L[(id.co1$ind[1]):(id.co1$ind[2])],
           bhat.table$delta.log.N1.plus[(id.co1$ind[1]):(id.co1$ind[2])],
           labels=row.names(bhat.table.list[[1]][(id.co1$ind[1]):
             (id.co1$ind[2]),]),
           cex= 0.7, pos=3)
      points(bhat.table$L[(id.co1$ind[1]):(id.co1$ind[2])],
             bhat.table$delta.log.N1.plus[(id.co1$ind[1]):(id.co1$ind[2])],
             col = colour.xy[xy], pch = 16)
      text(bhat.table$L[(id.co1$ind[2]+1):length(bhat.table$L)],
           bhat.table$delta.log.N1.plus[(id.co1$ind[2]+1):length(bhat.table$L)],
           labels=row.names(bhat.table.list[[1]][
             (id.co1$ind[2]+1):length(bhat.table$delta.log.N1.plus),]),
           cex= 0.7, pos=3)
      title(xlab = "L", outer = TRUE, line = 2.5)

      if(savePlots == TRUE) coho_plot[[xy]] <- recordPlot() ##coho_plot[xy] <- quartz.save(file = paste("Cohort",xy,"plot", sep='_'),type = "jpeg")

      repSel <- readline(prompt = writeLines("Are you satisfied with your selection and want to continue? \n'y' or 'redo':"))

      if(repSel == '') repSel <- readline(prompt = writeLines("Please type 'y' or 'redo' into the console! \nAre you satisfied with your selection and want to continue? \n'y' or 'redo':"))

      dev.off()

      if(repSel == 'y') break
    }

    #STEP 9: fills one value in seventh column and one in eigth column, get clean starting value
    if(length(id.co1$ind) == 0) break
    choosen <- seq(id.co1$ind[1],id.co1$ind[2],interval)
    # if odd, take midpoint:
    if(length(choosen) %% 2 == 1){
      mid_choice <- id.co1$ind[1] + (id.co1$ind[2]-id.co1$ind[1])/2
    }else{
      # if even, take one after midpoint:
      mid_choice <- (id.co1$ind[1] + (id.co1$ind[2]-id.co1$ind[1])/2) + (interval/2)
    }
    bhat.table$log.N1[mid_choice] <- log(bhat.table$N1.plus[mid_choice])
    bhat.table$N1[mid_choice] <- bhat.table$N1.plus[mid_choice]


    #STEP 10: fills rest of seventh column
    for(i in (mid_choice+1):length(bhat.table$delta.log.N)){
      log.N1 <- bhat.table$log.N1[i-1] + bhat.table$delta.log.N[i]
      bhat.table$log.N1[i] <- log.N1
    }

    #STEP 11: fills rest of eigth column
    for(i in (mid_choice+1):length(bhat.table$delta.log.N)){
      N1i <- round(exp(bhat.table$log.N1[i]),digits=2)
      bhat.table$N1[i] <- N1i
    }

    #STEP 12: fills ninth column
    bhat.table$N1[which(which(is.na(bhat.table$N1)) <
                          max.class.ind.co1)] <- bhat.table$N1.plus[
                            which(which(is.na(bhat.table$N1)) <
                                    max.class.ind.co1)]

    bhat.table$N2.plus <- bhat.table$N1.plus - bhat.table$N1

    bhat.table$N2.plus[which(is.na(bhat.table$N2.plus))] <-
      bhat.table$N1.plus[which(is.na(bhat.table$N2.plus))]


    #for creation of new bhat table
    bhat.table2 <- bhat.table
    bhat.table2$N2.plus[which(bhat.table2$N2.plus < 0)] <- 0
    bhat.table2$N1[which(bhat.table2$N2.plus < 0)] <- bhat.table2$N1.plus

    #creation of new table for continuation
    bhat.table <- data.frame('mean.length.classes' = bhat.table2$mean.length.classes,
                             'N1.plus' = bhat.table2$N2.plus,
                             'log.N1.plus' = NA,
                             'delta.log.N1.plus' = NA,
                             'L' = NA,
                             'delta.log.N' = NA,
                             'log.N1' = NA,
                             'N1' = NA,
                             'N2.plus' = NA)

    #arrange dataframe for saving
    colnames(bhat.table2) <- c('mean.length.classes',paste('N',xy,'.plus',sep=''),
                               paste('log.N',xy,'.plus',sep=''),
                               paste('delta.log.N',xy,'.plus',sep=''),'L',
                               'delta.log.N',paste('log.N',xy,sep=''),
                               paste('N',xy,sep=''),paste('N',xy+1,'.plus',sep=''))

    #save data
    bhat.table.list[[xy+1]] <- bhat.table2
    last_choices <- id.co1$ind
    a.b.list[[xy]] <- c(a.b.list[[xy]],id.co1$ind[1],id.co1$ind[2])
  }

  # Bhatta table
  bhat.table.list <- bhat.table.list[!sapply(bhat.table.list, is.null)]
  bhat.table.list2 <- bhat.table.list[-1]
  bhat.table.list_new <- lapply(bhat.table.list2, function(x)  return(x[,3:9]))
  bhat.results <- do.call(cbind, bhat.table.list_new)
  if(!is.null(bhat.results)){
    ret_pri <- cbind(bhat.table.list[[1]][,1:2], bhat.results)
  }else ret_pri <- bhat.table.list[[1]]

  # Intercept & Slope
  a.b.list <- a.b.list[!sapply(a.b.list,is.null)]
  if(length(a.b.list) > 0){
    a.b.df <- do.call(rbind,a.b.list)
    a.b.df <- cbind(1:length(a.b.df[,1]),a.b.df)
    colnames(a.b.df) <- c("Cohort","Intercept","Slope","Start_point","End_point")
  }else a.b.df = NA

  # Lmean and SD
  l.s.list <- l.s.list[!sapply(l.s.list,is.null)]
  if(length(l.s.list) > 1){
    l.s.df <- do.call(rbind,l.s.list)
    l.s.df <- cbind(1:length(l.s.df[,1]),l.s.df)
    # seperation index (SI)
    SIs <- rep(NA,length(l.s.df[,1]))
    for(i in 1:(length(l.s.df[,1])-1)){
      SIs[i] <- abs(l.s.df[i+1,2]-l.s.df[i,2]) / abs(l.s.df[i+1,3]-l.s.df[i,3])
    }
    l.s.df <- cbind(l.s.df,SIs)
    colnames(l.s.df) <- c("Cohort","Lmean","SD","SI")
  }else l.s.df <- NA

  temp.list <- temp.list[!sapply(temp.list, is.null)]
  coho_plot <- coho_plot[!sapply(coho_plot,is.null)]

  ret <- list(res,
              regression_lines = a.b.df,
              Lmean_SD_list = l.s.df,
              bhat_results = ret_pri,
              distributions = temp.list,
              cohort_plots = coho_plot)
  class(ret) <- "Bhattacharya"
  if(cancel != TRUE){
    plot(ret)
    return(ret)
  }
}
#' @title Powell-Wetherall method
#
#' @description A method to estimate the instantaneous total mortality rate (Z) and
#'    the infinite length of the von Bertalanffy growth equation
#'    (Powell, 1979; Wetherall et al., 1987).
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item  \code{midLengths}: midpoints of the length groups,
#'   \item \code{Linf}: infinite length for investigated species [cm],
#'   \item \code{K}: growth coefficent for investigated species [1/year],
#'   \item \code{t0}: theoretical time zero, at which individuals of this species hatch,
#'   \item \code{catch}: catch as vector, or a matrix with catches of subsequent years;
#' }
#' @param catch_columns optional; in case catch is a matrix or data.frame, a number
#'    indicating which column of the matrix should be analysed (Default: \code{NA}).
#' @param savePlots logical; if TRUE the plot is recorded. Default is FALSE.
#' @param reg_int instead of using the identity method a range can be determined,
#'    which is to be used for the regression analysis. If equal to NULL identity method
#'    is applied (default).
#' @param main title of plot (Default is "Powell-Wetherall plot")
#'
#' @keywords function mortality Z/K Linf
#'
#' @examples
#' \donttest{
#' data(synLFQ3)
#' powell_wetherall(synLFQ3)
#'
#' data(synLFQ5)
#' powell_wetherall(synLFQ5, catch_columns = 1:12)
#' }
#'
#' @details  The first length group or age class within the list object \code{midLengths} or
#'    \code{age} will be used as the Lprim or tprime (length of recruitment to fishery).
#'    This function includes the
#'    \code{identify} function, which asks you to choose two points from a graph manually. The
#'    two points which you choose by clicking on the plot in the graphical device represent
#'    the start and end of the data points, which should be used for the analysis. Based
#'    on these points the regression line is calculated. The Powell and Wetherall method
#'    only works with length-frequency data.
#'
#' @return A list with the input parameters and follwing objects:
#' \itemize{
#'   \item \strong{tmean} or \strong{Lmean}: mean age or length of fish,

#'   \item \strong{Z}: total mortality;}
#' and/or following objects when applying the Powell and Wetherall method:
#' \itemize{
#'   \item \strong{Lmean_Lprime}: dependent variable for regression analysis,
#'   \item \strong{Lprime}: some length for which all fish of that length and
#'      longer are under full exploitation,
#'   \item \strong{Linf_est}: infinite length in [cm] (Linf),
#'   \item \strong{se_Linf}: standard error of Linf,
#'   \item \strong{confidenceInt_Linf}: confidence interval for Linf,
#'   \item \strong{ZK}: total mortality divided by K (Z/K),
#'   \item \strong{se_ZK}: standard error of Z/K,
#'   \item \strong{confidenceInt_ZK}: confidence interval of Z/K;
#' }
#'
#' @importFrom grDevices dev.new recordPlot
#' @importFrom graphics abline identify par plot points
#' @importFrom stats lm qt
#' @importFrom utils flush.console
#'
#' @references
#' Powell, D.G., 1979. Estimation of mortality and growth parameters from the length-
#' frequency of a catch. \emph{Rapp.P.-v.Reun.CIEM}, 175:167-169
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' Wetherall, J.A., J.J. Polovina and S. Ralston, 1987. Estimating growth and mortality
#' in steady-state fish stocks from length-frequency data. \emph{ICLARM Conf.Proc.},
#' (13):53-74
#'
#' @export

powell_wetherall <- function(param, catch_columns = NA,
                             savePlots = FALSE, reg_int = NULL,
                             main = "Powell-Wetherall plot"){

  res <- param
  catch <- res$catch

  if(class(catch) == "data.frame" | class(catch) == "matrix"){
    if(is.na(catch_columns[1])){
      writeLines("By default the whole catch matrix is considered for this analysis. Please be aware that this \n method requires catches representing one year. You can choose separate columns of the catch \n matrix with 'catch_columns'.")
    }else{
      catchmat <- catch[,(catch_columns)]
      if(length(catch_columns) > 1){
        catch <- rowSums(catchmat, na.rm = TRUE)
      }else catch <- catchmat
    } # stop("Please provide a number indicating which column of the catch matrix should be analysed!")
  }

  if("midLengths" %in% names(res)){

    classes <- as.character(res$midLengths)
    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    Linf <- res$Linf
    K <- res$K

    # Error message if catch and age do not have same length
    if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
      if(length(classes) != length(catch[,1])) stop("Ages and catch do not have the same length!")
    }else if(class(catch) == 'numeric'){
      if(length(classes) != length(catch)) stop("Ages and catch do not have the same length!")
    }

    # calculate cumulative catch
    cumCatch <- rev(cumsum(rev(catch)))

    # calculate  C * (L1 + L2) / 2
    c_midlength <- catch * classes.num

    # calculate L prime   ===   x
    interval <-  (classes.num[2] - classes.num[1]) / 2
    Lprime <- classes.num - interval

    # calculate L mean
    sum_midL_c <- rep(NA,length(classes.num))
    Lmean <- rep(NA,length(classes.num))
    for(i in 1:length(c_midlength)){
      sum_midL_c[i] <- sum(c_midlength[i:length(c_midlength)])
      Lmean[i] <- sum(c_midlength[i:length(c_midlength)]) /
        sum(catch[i:length(catch)])
    }

    # calculate Lmean - Lprime
    Lmean_Lprime <- Lmean - Lprime

    #identify plot
    if(is.null(reg_int)){
      repeat{
        writeLines("Please choose the minimum and maximum point in the \ngraph to include for the regression line!")
        flush.console()
        dev.new()#noRStudioGD = TRUE)
        plot(x = Lprime,y = Lmean_Lprime,
             xlab = "Lprime", ylab = "Lmean - Lprime")
        text(Lprime+0.5, Lmean_Lprime+0.5, labels=as.character(order(Lprime)), cex= 0.7)
        mtext(side = 3, "Click on two points. Escape to Quit.",
              xpd = NA, cex = 1.25)
       cutter <- identify(x = Lprime,y = Lmean_Lprime,
                           labels = order(Lprime), n=2)

        if(length(cutter) == 0){
          stop(noquote("You did not choose any points! Please run the function again \nand choose points to include into the estimation of Z."))

        }

        length.cutter <- length(cutter[1]:cutter[2])
        # Break loop if selection does not embrace at least 3 points
        if(length.cutter < 3){
          writeLines("Your selection is not possible. You have to choose two \npoints which include at least one other point. At least \nthree points are required for a regression line. Please choose again!")
          flush.console()
        }
        if(length.cutter >= 3){
          break
        }
      }
    }
    if(!is.null(reg_int)){
      cutter <- reg_int
    }
    if(length(cutter) != 2) stop("You have to provide 2 numbers in reg_int.")

    #calculations + model
    df.BH <- as.data.frame(cbind(classes.num,Lmean_Lprime,Lprime))
    df.BH.cut <- df.BH[cutter[1]:cutter[2],]
    lm1 <- lm(Lmean_Lprime ~ Lprime, data = df.BH.cut)
    sum_lm1 <- summary(lm1)
    r_lm1 <- sum_lm1$r.squared
    intercept_lm1 <- sum_lm1$coefficients[1]
    slope_lm1 <- sum_lm1$coefficients[2]
    se_slope_lm1 <- sum_lm1$coefficients[4]
    se_intercept_lm1 <- sum_lm1$coefficients[3]

    #fit of regression line
    lm1.fit <- sum_lm1$r.squared

    SE_slope <- abs(se_slope_lm1)
    confi_slope <- abs(se_slope_lm1) * qt(0.975,sum_lm1$df[2])
    conf_slope <- slope_lm1 + c(-confi_slope,confi_slope)

    SE_intercept <- abs(se_intercept_lm1)
    confi_intercept <- abs(se_intercept_lm1) * qt(0.975,sum_lm1$df[1])
    conf_intercept <- intercept_lm1 + c(-confi_intercept,confi_intercept)

    # Linf with SE and confidence interval
    Linf.BH <- (-intercept_lm1 / slope_lm1)
    se_Linf.BH <- (abs(SE_intercept)/abs(intercept_lm1) +
                     abs(SE_slope)/abs(slope_lm1)) *
      (abs(intercept_lm1) / abs(slope_lm1))
    confi_Linf <- (abs(SE_intercept)/abs(intercept_lm1) +
                     abs(SE_slope)/abs(slope_lm1)) *
      (abs(intercept_lm1) / abs(slope_lm1)) * qt(0.975,sum_lm1$df[2])
    conf_Linf.BH <- Linf.BH + c(-confi_Linf,confi_Linf)

    # Z/K with SE and confidence interval
    ZK.BH <- (-(1+slope_lm1)/slope_lm1)
    se_ZK.BH <- abs(SE_slope)
    confi_ZK <- se_ZK.BH * qt(0.975,sum_lm1$df[2])
    conf_ZK.BH <- ZK.BH + c(-confi_ZK,confi_ZK)

    #final plot
    plot(x = Lprime,y = Lmean_Lprime,
         xlab = "Lprime", ylab = "Lmean - Lprime",
         cex = 1.5, main = main)
    par(new=T)
    points(x = df.BH.cut$Lprime,y = df.BH.cut$Lmean_Lprime,
           pch = 19, col = 'blue', cex = 1.5)
    abline(a=intercept_lm1,b=slope_lm1,col="blue",lwd = 1.7)
    if(savePlots == TRUE){
      ploti <- recordPlot()
    }else ploti = NA


    #save all in list
    ret <- c(res,list(
      Lmean_Lprime = Lmean_Lprime,
      Lprime = Lprime,
      Linf_est = Linf.BH,
      se_Linf = se_Linf.BH,
      confidenceInt_Linf = conf_Linf.BH,
      ZK = ZK.BH,
      se_ZK = se_ZK.BH,
      confidenceInt_ZK = conf_ZK.BH,
      plot = ploti
    ))
    return(ret)
  }
}
#' @title Von Bertalanffy Growth function (VBGF)
#
#' @description  This function applies the von Bertalanffy growth function (VBGF).
#'    It allows to calculate ages from lengths or lengths from ages based on the special,
#'    generalised or seasonalised VBGF.
#'
#' @param t ages for which to calculate corresponding lengths, or
#' @param L lengths for which to calculate corresponding ages
#' @param param a list with following potential objects:
#' \itemize{
#'   \item \code{Linf}: infinite length for investigated species in cm, or
#'   \item \code{Winf}: infinite weight for investigated species in gramm
#'   \item \code{K}: growth coefficent for investigated species per year
#'   \item \code{t0}: theoretical time zero, at which individuals of this species hatch (default: 0)
#'   \item \code{b}: exponent of weight length relationship (default: 3)
#'   \item \code{D}: surface factor (default: 1)
#'   \item \code{L0}: length at hatching for VBGF with L0
#'   \item \code{ts}: onset of the first oscillation relative to t0
#'   \item \code{C}: intensity of (sinusoid) growth oscillations. Default is no oscillation (C = 0)
#' }
#'
#' @keywords function growth VBGF
#'
#' @examples
#' # calculation of lengths
#' # with t0
#' t <- seq(0,6,0.1)
#' Lt <- VBGF(list(Linf=80, K=0.6, t0=-0.1),t=t)
#' plot(t, Lt, t="l")
#'
#' # with L0
#' t <- seq(0,6,0.1)
#' Lt <- VBGF(list(Linf=80, K=0.6, L0=2),t=t)
#' plot(t, Lt, t="l")
#'
#' # with Winf
#' t <- seq(0,6,0.1)
#' Wt <- VBGF(list(Winf=4000, K=0.8), t=t)
#' plot(t, Wt, t="l")
#'
#' # seasonalised VBGF
#' t <- seq(0,6,0.1)
#' Lt <- VBGF(list(Linf=80, K=0.6, t0=-0.1, ts=0.5, C=0.75),t=t)
#' plot(t, Lt, t="l")
#'
#'
#' # calculation of ages
#' L <- seq(2,200,0.1)
#' t <- VBGF(L = L, list(Linf=210, K=0.8, C= 0.75))
#' plot(t, L, t="l")
#'
#' @details
#' Based upon which input parameters are given one of the following
#' VBGF types is applied: "special", "generalised", or "seasonalised" VBGF.
#'
#' @return A vector with estimated lengths corresponding to provided ages.
#'
#' @references
#' Somers, I. F. (1988). On a seasonally oscillating growth function.
#' Fishbyte, 6(1), 8-11
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

VBGF <- function(param, t = NA, L = NA){
  res <- param
  if(is.na(t[1]) & is.na(L[1])) stop("Either length L or age t has to be provided to calculate the corresponding.")
  Linf <- ifelse("Linf" %in% names(res),res$Linf, NA)
  Winf <- ifelse("Winf" %in% names(res),res$Winf, NA)
  if("K" %in% names(res)){
    K <- res$K
  }else stop("Please provide a K parameter in 'param'.")
  t0 <- ifelse("t0" %in% names(res),res$t0, 0)
  b <- ifelse("b" %in% names(res),res$b, 3)
  D <- ifelse("D" %in% names(res),res$D, 1)
  L0 <- ifelse("L0" %in% names(res),res$L0, NA)
  ts <- ifelse("ts" %in% names(res),res$ts, 0)
  C <- ifelse("C" %in% names(res),res$C, 0)
  if("ta" %in% names(res)) t0 <- res$ta


  if(is.na(Linf) & is.na(Winf)) stop("You have to provide either Linf or Winf.")
  if(is.na(L[1]) & is.na(t[1])) stop("Please provide a vector with either potential ages or lengths.")

  if(is.na(L0)){
    # generalised seasonalised VBGF for length
    if((is.na(Winf) & is.na(L[1])) |
       (!is.na(Winf) & !is.na(Linf) & is.na(L[1]))){
      res <- Linf * (1 - exp(-(K * D * (t - t0) + (((C*K*D)/(2*pi)) * sin(2*pi*(t-ts))) -
                                  (((C*K*D)/(2*pi)) * sin(2*pi*(t0-ts)))))) ^ (1/D)
    }

    # OLD:
    # res <- Linf * (1 - exp(-K * D * (t - t0) -
    #                                       ((C*K*D)/(2*pi)) *
    #                                       ((sin(2*pi*(t-ts))) +
    #                                          (sin(2*pi*(t0-ts))))))^ (1/D)


    if((is.na(Winf) & is.na(t[1])) |
       (!is.na(Winf) & !is.na(Linf) & is.na(t[1]))){
      if(D == 1 & C == 0) res <- t0 - (log(1-L/Linf)/K)
      # lookup table for soVBGF
      if(D != 1 | C != 0){
        tmax <- (t0 * K * D - log(1 - exp(D*log(L[length(L)]/Linf)))) / (K*D)
        if(is.na(tmax)) tmax <- 40
        lookup_age <- seq((t0 - 1),(tmax+10),0.001)
        lookup_length <-  Linf * (1 - exp(-(K * D * (lookup_age - t0) + (((C*K*D)/(2*pi)) *
                                                                          sin(2*pi*(lookup_age-ts))) -
                                            (((C*K*D)/(2*pi)) * sin(2*pi*(t0-ts)))))) ^ (1/D)

        # OLD:
        # lookup_length <-  Linf * (1 - exp(-K * D * (lookup_age - t0) -
        #                                     ((C*K*D)/(2*pi)) *
        #                                     ((sin(2*pi*(lookup_age-ts))) +
        #                                        (sin(2*pi*(t0-ts))))))^ (1/D)

        lookup_ind <- lapply(X = L, FUN = function(x) which.min((lookup_length - x)^2))
        res <- lookup_age[unlist(lookup_ind)]
      }
    }


    # generalised VBGF for weight
    if(is.na(Linf) & is.na(L[1])) res <- Winf * (1 - exp(-K * D * (t - t0)))^(b/D)
    if(is.na(Linf) & is.na(t[1])) res <- (D * K * t0 - log(1 - exp((D * log(L/Winf))/(b)))) / (K*D)
  }

  # special VBGF for length with L0
  if(!is.na(L0) & is.na(L[1])) res <- Linf - (Linf - L0) * exp(-K * t)
  if(!is.na(L0) & is.na(t[1])) res <- log((L - Linf) / -(Linf - L0)) / -K

  return(res)
}

#' @title Fitting VBGF growth curves through lfq data
#'
#' @description Thsi function estimates von Bertalanffy growth function (VBGF)
#'    curves for a set of growth parameters.
#'
#' @param lfq a list of the class "lfq" consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row)
#'      and sampling date (column),
#'   \item \strong{rcounts} restructured frequencies,
#'   \item \strong{peaks_mat} matrix with positive peaks with distinct values,
#'   \item \strong{ASP} available sum of peaks, sum of posititve peaks which
#'      could be potential be hit by growth curves;
#' }
#' @param par a list with following growth parameters:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default: 100),
#'   \item \strong{K} curving coefficient (default: 0.1),
#'   \item \strong{ta} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0.25),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param agemax maximum age of species; default NULL, then estimated from Linf
#' @param flagging.out logical; should positive peaks be flagged out?
#'    (Default : TRUE)
#' @param lty The line type. Line types can either be specified as an integer
#'    (0=blank, 1=solid, 2=dashed (default), 3=dotted, 4=dotdash, 5=longdash,
#'    6=twodash) or as one of the character strings "blank", "solid", "dashed",
#'    "dotted", "dotdash", "longdash", or "twodash", where "blank" uses 'invisible
#'    lines' (i.e., does not draw them).
#' @param lwd The line width, a positive number, defaulting to 2. The interpretation
#'    is device-specific, and some devices do not implement line widths less
#'    than one. (See the help on the device for details of the interpretation.)
#' @param col A specification for the default plotting color. See section
#'    'Color Specification'.
#' @param draw logical; indicating whether growth curves should be added to
#'    existing lfq plot
#' @param tincr step for plotting
#'
#' @examples
#' data(synLFQ5)
#' res <- lfqRestructure(synLFQ5, MA=11)
#' plot(res)
#' tmp <- lfqFitCurves(res, par=list(Linf=80,K=0.5,ta=0.25), draw=TRUE)
#'
#' @details \code{ta} subsitutes the starting point from known from Fisat 2.
#'    This parameter is necessary for anchoring the growth curves on the time axis.
#'    It does not subsitute \code{t0}. However, it corresponds to the peak spawning
#'    of the species (x intercept of growth curve) and has values between 0 and 1,
#'    where 0 corresponds to spawning at the 1st of January and 0.999 corresponds to the
#'    31st of December. The default value of 0.25 or 3/12 corresponds the third month
#'    of the year, March.
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{Lt}: dataframe with ages and lengths of the cohorts,
#'   \item \strong{agemax}: maximum age of species.
#'   \item \strong{ncohort}: number of cohorts,
#'   \item \strong{ASP}: available sum of peaks, sum of posititve peaks which
#'   could be potential be hit by growth curves. This is calculated as the sum of
#'   maximum values from each run of posive restructured scores.
#'   \item \strong{ESP}: available sum of peaks,
#'   \item \strong{fASP}: available sum of peaks,
#'   \item \strong{fESP}: available sum of peaks,
#' }
#'
#'
#' @references
#' Brey, T., Soriano, M., and Pauly, D. 1988. Electronic length frequency analysis: a revised and expanded
#' user's guide to ELEFAN 0, 1 and 2.
#'
#' Pauly, D. 1981. The relationship between gill surface area and growth performance in fish:
#' a generalization of von Bertalanffy's theory of growth. \emph{Meeresforsch}. 28:205-211
#'
#' Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the objective extraction of
#' growth parameters from length-frequency data. \emph{Meeresforschung}, 28(4):205-211
#'
#' Pauly, D., 1985. On improving operation and use of ELEFAN programs. Part I: Avoiding
#' "drift" of K towards low values. \emph{ICLARM Conf. Proc.}, 13-14
#'
#' Pauly, D., 1987. A review of the ELEFAN system for analysis of length-frequency data in
#' fish and aquatic invertebrates. \emph{ICLARM Conf. Proc.}, (13):7-34
#'
#' Pauly, D. and G. R. Morgan (Eds.), 1987. Length-based methods in fisheries research.
#' (No. 13). WorldFish
#'
#' Pauly, D. and G. Gaschuetz. 1979. A simple method for fitting oscillating length growth data, with a
#' program for pocket calculators. I.C.E.S. CM 1979/6:24. Demersal Fish Cttee, 26 p.
#'
#' Pauly, D. 1984. Fish population dynamics in tropical waters: a manual for use with programmable
#' calculators (Vol. 8). WorldFish.
#'
#' Quenouille, M. H., 1956. Notes on bias in estimation. \emph{Biometrika}, 43:353-360
#'
#' Somers, I. F., 1988. On a seasonally oscillating growth function. ICLARM Fishbyte 6(1): 8-11.
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2): 407 p.
#'
#' Tukey, J., 1958. Bias and confidence in not quite large samples.
#' \emph{Annals of Mathematical Statistics}, 29: 614
#'
#' Tukey, J., 1986. The future of processes of data analysis. In L. V. Jones (Eds.),
#' The Collected Works of John W. Tukey-philosophy and principles of data analysis:
#' 1965-1986 (Vol. 4, pp. 517-549). Monterey, CA, USA: Wadsworth & Brooks/Cole
#'
#' @export

lfqFitCurves <- function(lfq,
  par = list(Linf = 100, K = 0.1, ta = 0.25, C = 0, ts = 0),
  agemax = NULL, flagging.out = TRUE,
  lty = 2, lwd = 1, col = 1,
  draw = FALSE, tincr = 0.05
){

  if(is.null(par$Linf) | is.null(par$K) | is.null(par$ta)) stop("Linf, K and ta have to provided in par.")

  # ISSUE: if seasonalised max age can be different and 0.95 very coarse
  if(is.null(agemax)){
    agemax <- ceiling((1/-par$K)*log(1-((par$Linf*0.95)/par$Linf)))
  }

  yeardec <- date2yeardec(lfq$dates)

  tmax <- max(yeardec, na.rm = TRUE) # maximum sample date
  tmin <- min(yeardec, na.rm = TRUE) # minimum sample date
  ncohort <- agemax + ceiling(diff(range(yeardec))) # number of cohorts
  tA.use <- par$ta # VBGF anchor time in year
  tAs <- seq(floor(tmax-ncohort), floor(tmax), by=1) + (tA.use+1e6)%%1   # anchor times
  tAs <- tAs[which(tAs < tmax)]
  ncohort <- length(tAs)


  par2 <- par
  t <- c(yeardec,max(yeardec)+0.2)  # for plotting, otherwise growth curves are cut at last sampling time
  Lt <- vector(mode="list", ncohort)
  for(ct in seq(ncohort)){
    par2$ta <- tAs[ct]
    rel.age <- t-tAs[ct] # relative age to anchor time
    t.use <- which(rel.age <= agemax & rel.age > 0)
    t.ct <- t[t.use]
    rel.age <- rel.age[t.use]
    if(length(t.ct) > 0){
      Lt.ct <- VBGF(param = par2, t = t.ct) ## do.call(what = VBGF,  par2)    # Lt.ct <- VBGF(lfq = par2, t = yeardec)   #
      Lt[[ct]] <- data.frame(t=t.ct, Lt=Lt.ct, ct=ct, rel.age=rel.age)
      if(draw){
        tmp <- par2
        t_tmp <- seq(tAs[ct], max(t.ct)+tincr, by=tincr)
        tmp$L <- VBGF(tmp, t = t_tmp) ## do.call(what = VBGF,  tmp)
        tmp$t <- yeardec2date(t_tmp)
        lines(L ~ t, tmp, lty=lty, lwd=lwd, col=col)
      }
    }
  }
  Lt <- do.call(rbind, Lt)
  Lt <- Lt[which(Lt$Lt>=0),]# trim negative lengths


  # calc scores
  grd <- expand.grid(Lt=lfq$midLengths, t=date2yeardec(lfq$dates))
  grd$Fs <- c(lfq$rcounts)  # turn matrix into 1-D vector
  grd$cohort_peaks <- c(lfq$peaks_mat) # turn matrix into 1-D vector
  grd$hit <- 0 # empty vector to record bins that are "hit" by a growth curve
  bin.width <- diff(lfq$midLengths) # bin width (should allow for uneven bin sizes)
  grd$bin.lower <- lfq$midLengths - (c(bin.width[1], bin.width)/2) # upper bin limit
  grd$bin.upper <- lfq$midLengths + (c(bin.width, bin.width[length(bin.width)])/2) # lower bin limit

  # mark all length classes (in all sampling times) which are hit
  tmp <-
    outer(X = c(grd$t), Y = Lt$t, FUN = "==") &
    outer(X = c(grd$bin.lower), Y = Lt$Lt, FUN = "<=") &
    outer(X = c(grd$bin.upper), Y = Lt$Lt, FUN = ">"
  )

  tmp2 <- apply(X = tmp, MARGIN = 1, FUN = sum, na.rm = TRUE)
  if(flagging.out){
    grd$hit <- tmp2
  } else {
    grd$hit <- as.numeric(tmp2 > 0)
  }

  # Count only one crossing of positive peaks within a cohort peak
  if(flagging.out){
    ch <- unique(grd$cohort_peaks)
    if(0 %in% ch) ch <- ch[-which(0 %in% ch)]
    for(ci in seq(length(ch))){
      chi <- ch[ci]
      peaki <- which(grd$cohort_peaks == chi)
      dpch <- grd$hit[peaki]
      if(sum(dpch, na.rm = TRUE) > 1){
        grd$hit[peaki] <- 0
        maxi <- which.max(grd$Fs[peaki])
        grd$hit[peaki[maxi]] <- 1
      }
    }
  }


  ESP <- sum(grd$Fs * grd$hit, na.rm = TRUE)
  fASP <- ESP/lfq$ASP
  fESP <- (10^(ESP/lfq$ASP)) /10

  lfq$Lt <- Lt
  lfq$agemax <- agemax
  lfq$ncohort <- ncohort
  lfq$ESP <- ESP
  lfq$fASP <- fASP
  lfq$fESP <- fESP
  lfq$Rn_max <- fESP

  return(lfq)
}
#' @title ELEFAN
#'
#' @description Electronic LEngth Frequency ANalysis for estimating growth parameter.
#'
#' @param x a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row)
#'      and sampling date (column);
#' }
#' @param Linf_fix numeric; if used the K-Scan method is applied with a fixed
#'    Linf value (i.e. varying K only).
#' @param Linf_range numeric vector with potential Linf values. Default is the
#'    last length class plus/minus 5 cm
#' @param K_range K values for which the score of growth functions should be
#'    calculated
#'    (by default: exp(seq(log(0.1),log(10),length.out = 100)))
#' @param C growth oscillation amplitude (default: 0)
#' @param ts onset of the first oscillation relative to t0 (summer point, default: 0)
#' @param MA number indicating over how many length classes the moving average
#'    should be performed (default: 5, for
#'    more information see \link{lfqRestructure}).
#' @param addl.sqrt Passed to \link{lfqRestructure}. Applied an additional square-root transformation of positive values according to Brey et al. (1988).
#'    (default: FALSE, for more information see \link{lfqRestructure}).
#' @param agemax maximum age of species; default NULL, then estimated from Linf
#' @param flagging.out logical; should positive peaks be flagged out? (Default : TRUE)
#' @param method Choose between the old FiSAT option to force VBGF crossing of a pre-defined
#'    bin (method = "cross"), or the more sophisticated (but computationally expensive) option
#'    to solve for ta via a maximisation of reconstructed score
#'    (default: method = "optimise").
#' @param cross.date Date. For use with \code{method = "cross"}. In combination
#'    with \code{cross.midLength}, defines the date of the crossed bin.
#' @param cross.midLength Numeric. For use with \code{method = "cross"}. In combination
#'    with \code{cross.date}, defines the mid-length of the crossed bin.
#' @param cross.max logical. For use with \code{method = "cross"}. Forces growth function
#'    to cross the bin with maximum positive score.
#' @param hide.progressbar logical; should the progress bar be hidden? (default: FALSE)
#' @param plot logical; indicating if plot with restructured frequencies and growth curves should
#'    be displayed
#' @param contour if used in combination with response surface analysis, contour lines
#'    are displayed rather than the score as text in each field of the score plot. Usage
#'    can be logical (e.g. TRUE) or by providing a numeric which indicates the
#'    number of levels (\code{nlevels} in \code{\link{contour}}). By default FALSE.
#' @param add.values logical. Add values to Response Surface Analysis plot (default: TRUE).
#'    Overridden when \code{contour = TRUE}.
#' @param rsa.colors vector of colors to be used with the Response Surface Analysis plot.
#'    (default: terrain.colors(20))
#' @param plot_title logical; indicating whether title to score plots should be displayed
#'
#' @examples
#' \donttest{
#' data(alba)
#'
#' ### Response Surface Analysis (varies Linf and K) ###
#'
#' # 'cross' method (used in FiSAT)
#' fit1 <- ELEFAN(
#'    x = alba, method = "cross",
#'    Linf_range = seq(from = 10, to = 20, length.out = 10),
#'    K_range = exp(seq(from = log(0.1), to = log(2), length.out = 20)),
#'    cross.date = alba$dates[3], cross.midLength = alba$midLengths[4],
#'    contour = TRUE
#' )
#' fit1$Rn_max; unlist(fit1$par)
#' plot(fit1); points(alba$dates[3], alba$midLengths[4], col=2, cex=2, lwd=2)
#'
#' # 'cross' method (bin with maximum score crossed)
#' fit2 <- ELEFAN(
#'    x = alba, method = "cross",
#'    Linf_range = seq(from = 10, to = 20, length.out = 20),
#'    K_range = exp(seq(from = log(0.1), to = log(2), length.out = 20)),
#'    cross.max = TRUE,
#'    contour = TRUE
#' )
#' fit2$Rn_max; unlist(fit2$par)
#' plot(fit2); points(alba$dates[7], alba$midLengths[9], col=2, cex=2, lwd=2)
#'
#'
#' # 'optimise' method (default)
#' fit3 <- ELEFAN(
#'    x = alba, method = "optimise",
#'    Linf_range = seq(from = 10, to = 20, length.out = 10),
#'    K_range = exp(seq(from = log(0.1), to = log(2), length.out = 20)),
#'    contour = TRUE
#' )
#' fit3$Rn_max; unlist(fit3$par)
#' plot(fit3)
#'
#'
#' ### K-Scan (varies K, Linf is fixed) ###
#'
#' # 'cross' method
#' fit4 <- ELEFAN(
#'    x = alba, method = "cross",
#'    Linf_fix = 10,
#'    K_range = round(exp(seq(from = log(0.1), to = log(2), length.out = 50)),2),
#'    cross.date = alba$dates[3], cross.midLength = alba$midLengths[4]
#' )
#' fit4$Rn_max; unlist(fit4$par)
#' plot(fit4); points(alba$dates[3], alba$midLengths[4], col=2, cex=2, lwd=2)
#'
#' }
#'
#' @details This functions allows to perform the K-Scan and Response surface
#'    analysis to estimate growth parameters.
#'    It combines the step of restructuring length-frequency data
#'    (\link{lfqRestructure}) followed by the fitting of VBGF
#'    curves through the restructured data (\link{lfqFitCurves}). K-Scan is a
#'    method used to search for the K
#'    parameter with the best fit while keeping the Linf fixed. In contrast,
#'    with response surface analysis
#'    both parameters are estimated and the fits are displayed in a heatmap.
#'    Both methods use \code{\link[stats]{optimise}} to find the best \code{ta} value
#'    for each combination of \code{K} and \code{Linf}. To find out more about
#'    \code{ta}, please refer to the Details description of
#'    \code{\link{lfqFitCurves}}. The score value \code{Rn_max} is comparable with
#'    the score values of the other ELEFAN functions (\code{\link{ELEFAN_SA}} or
#'    \code{\link{ELEFAN_GA}}) when other settings are consistent
#'    (e.g. `MA`, `addl.sqrt`, `agemax`, `flagging.out`).
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{rcounts}: restructured frequencies,
#'   \item \strong{peaks_mat}: matrix with positive peaks with distinct values,
#'   \item \strong{ASP}: available sum of peaks, sum of posititve peaks which
#'      could be potential be hit by
#'      growth curves,
#'   \item \strong{score_mat}: matrix with scores for each Linf (only Linf_fix)
#'    and K combination,
#'   \item \strong{ncohort}: number of cohorts used for estimation,
#'   \item \strong{agemax}: maximum age of species,
#'   \item \strong{par}: a list with the parameters of the von Bertalanffy growth
#'      function:
#'      \itemize{
#'        \item \strong{Linf}: length infinity in cm,
#'        \item \strong{K}: curving coefficient;
#'        \item \strong{ta}: time point anchoring growth curves in year-length
#'          coordinate system, corrsponds to peak spawning month,
#'        \item \strong{C}: amplitude of growth oscillation
#'          (if \code{seasonalised} = TRUE),
#'        \item \strong{ts}: summer point of oscillation (ts = WP - 0.5)
#'          (if \code{seasonalised} = TRUE),
#'        \item \strong{phiL}: growth performance index defined as
#'          phiL = log10(K) + 2 * log10(Linf);
#'      }
#'   \item \strong{Rn_max}: highest score value
#' }
#'
#' @importFrom grDevices terrain.colors
#' @importFrom graphics abline axis grid image mtext par plot text title
#' @importFrom utils setTxtProgressBar txtProgressBar flush.console
#' @importFrom stats optimise
#'
#' @references
#' Brey, T., Soriano, M., and Pauly, D. 1988. Electronic length frequency analysis: a revised and expanded
#' user's guide to ELEFAN 0, 1 and 2.
#'
#' Pauly, D. 1981. The relationship between gill surface area and growth performance in fish:
#' a generalization of von Bertalanffy's theory of growth. \emph{Meeresforsch}. 28:205-211
#'
#' Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the objective extraction of
#' growth parameters from length-frequency data. \emph{Meeresforschung}, 28(4):205-211
#'
#' Pauly, D., 1985. On improving operation and use of ELEFAN programs. Part I: Avoiding
#' "drift" of K towards low values. \emph{ICLARM Conf. Proc.}, 13-14
#'
#' Pauly, D., 1987. A review of the ELEFAN system for analysis of length-frequency data in
#' fish and aquatic invertebrates. \emph{ICLARM Conf. Proc.}, (13):7-34
#'
#' Pauly, D. and G. R. Morgan (Eds.), 1987. Length-based methods in fisheries research.
#' (No. 13). WorldFish
#'
#' Pauly, D. and G. Gaschuetz. 1979. A simple method for fitting oscillating length growth data, with a
#' program for pocket calculators. I.C.E.S. CM 1979/6:24. Demersal Fish Cttee, 26 p.
#'
#' Pauly, D. 1984. Fish population dynamics in tropical waters: a manual for use with programmable
#' calculators (Vol. 8). WorldFish.
#'
#' Quenouille, M. H., 1956. Notes on bias in estimation. \emph{Biometrika}, 43:353-360
#'
#' Somers, I. F., 1988. On a seasonally oscillating growth function. ICLARM Fishbyte 6(1): 8-11.
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2): 407 p.
#'
#' Tukey, J., 1958. Bias and confidence in not quite large samples.
#' \emph{Annals of Mathematical Statistics}, 29: 614
#'
#' Tukey, J., 1986. The future of processes of data analysis. In L. V. Jones (Eds.),
#' The Collected Works of John W. Tukey-philosophy and principles of data analysis:
#' 1965-1986 (Vol. 4, pp. 517-549). Monterey, CA, USA: Wadsworth & Brooks/Cole
#'
#' @export

ELEFAN <- function(lfq, Linf_fix = NA, Linf_range = NA,
                   K_range = exp(seq(log(0.1), log(10), length.out=100)),
                   C = 0, ts = 0,
                   MA = 5, addl.sqrt = FALSE,
                   agemax = NULL, flagging.out = TRUE, method = "optimise",
                   cross.date = NULL, cross.midLength = NULL, cross.max = FALSE,
                   hide.progressbar = FALSE,
                   plot = FALSE, contour = FALSE,
                   add.values = TRUE, rsa.colors = terrain.colors(20),
                   plot_title = TRUE){

    res <- lfq
    classes <- res$midLengths
    catch <- res$catch
    dates <- res$dates
    interval <- (classes[2]-classes[1])/2

    n_samples <- dim(catch)[2]
    n_classes <- length(classes)

    ## ts <- WP - 0.5
    if(is.na(Linf_fix) & is.na(Linf_range[1])) Linf_range <- seq(classes[n_classes]-5,classes[n_classes]+5,1) ### OLD: c(classes[n_classes]-5,classes[n_classes]+5)

    ## ELEFAN 0
    res <- lfqRestructure(res, MA = MA, addl.sqrt = addl.sqrt)
    catch_aAF_F <- res$rcounts
    peaks_mat <- res$peaks_mat
    ASP <- res$ASP

    if(!is.na(Linf_fix)){
        Linfs <- Linf_fix
    }else Linfs <-  Linf_range
    Ks <- K_range




    ## optimisation function
    sofun <- function(tanch, lfq, par, agemax, flagging.out){
        Lt <- lfqFitCurves(lfq,
                           par=list(Linf=par[1], K=par[2], ta=tanch,
                                    C=par[4], ts=par[5]),
                           flagging.out = flagging.out, agemax = agemax)
        return(Lt$ESP)
    }

    ESP_tanch_L <- matrix(NA,nrow=length(Ks),ncol=length(Linfs))
    ESP_list_L <- matrix(NA,nrow=length(Ks),ncol=length(Linfs))
    writeLines(paste("Optimisation procuedure of ELEFAN is running. \nThis will take some time. \nThe process bar will inform you about the process of the calculations.",sep=" "))
    flush.console()


    ## check that both cross.date and cross.midLength are identified
    if(method == "cross"){
        if(cross.max){ ## overrides 'cross.date' and 'cross.midLength'
            max.rcount <- which.max(res$rcounts)
            COMB <- expand.grid(midLengths = rev(rev(res$midLengths)), dates = res$dates)
            cross.date <- COMB$dates[max.rcount]
            cross.midLength <- COMB$midLengths[max.rcount]
        } else {
            if(is.null(cross.date) | is.null(cross.midLength)){
                stop("Must define both 'cross.date' and 'cross.midLength' when 'cross.max' equals FALSE")
            }
        }
    }


    if(!hide.progressbar){
        nlk <- prod(dim(ESP_list_L))
        pb <- txtProgressBar(min=1, max=nlk, style=3)
        counter <- 1
    }
    for(li in 1:length(Linfs)){

        ESP_tanch_K <- rep(NA,length(Ks))
        ESP_list_K <- rep(NA,length(Ks))
        for(ki in 1:length(Ks)){

            ## determine agemax for Linf and K combination if not already defined
            if(is.null(agemax)){
                agemax.i <- ceiling((1/-Ks[ki])*log(1-((Linfs[li]*0.95)/Linfs[li])))
            } else {
                agemax.i <- agemax
            }


            if(method == "cross"){ ## Method for crossing center of a prescribed bin

                t0s <- seq(floor(min(date2yeardec(res$dates)) - agemax.i), ceiling(max(date2yeardec(res$dates))), by = 0.01)

                Ltlut <- Linfs[li] * (1-exp(-(
                    Ks[ki]*(date2yeardec(cross.date)-t0s)
                    + (((C*Ks[ki])/(2*pi))*sin(2*pi*(date2yeardec(cross.date)-ts)))
                    - (((C*Ks[ki])/(2*pi))*sin(2*pi*(t0s-ts)))
                )))

                ta <- t0s[which.min((Ltlut - cross.midLength)^2)] %% 1

                resis <- sofun(tanch = ta, lfq = res, par = c(Linfs[li], Ks[ki], NA, C, ts), agemax = agemax.i, flagging.out = flagging.out)
                ESP_list_K[ki] <- resis
                ESP_tanch_K[ki] <- ta

            }

            ## Optimised method for searching best scoring ta
            if(method == "optimise"){
                resis <- optimise(
                    f = sofun,
                    lower = 0,
                    upper = 1,
                    lfq = res,
                    par = c(Linfs[li], Ks[ki], NA, C, ts),
                    agemax = agemax.i,
                    flagging.out = flagging.out,
                    tol = 0.001,
                    maximum = TRUE
                )
                ESP_list_K[ki] <- resis$objective
                ESP_tanch_K[ki] <- resis$maximum
            }

            ## update counter and progress bar
            if(!hide.progressbar){
                setTxtProgressBar(pb, counter)
                counter <- counter + 1
            }
        }
        ESP_list_L[,li] <- ESP_list_K
        ESP_tanch_L[,li] <- ESP_tanch_K
    }

    dimnames(ESP_list_L) <- list(Ks,Linfs)
    score_mat <- round((10^(ESP_list_L/ASP)) /10,digits = 3)
    rownames(score_mat) <- round(as.numeric(rownames(score_mat)), digits = 2)
    dimnames(ESP_tanch_L) <- list(Ks,Linfs)


    ## Graphs
    if(is.na(Linf_fix)){
        plot_dat <- reshape2::melt(score_mat)

        image(
            x = Linfs,
            y = Ks,
            z = t(score_mat), col=rsa.colors,
            ylab = 'K', xlab='Linf'
        )

        if(plot_title)  title('Response surface analysis', line = 1)

        if(contour){
            contour(x = Linfs, y = Ks, z = t(score_mat), add = TRUE)
        } else {
            if(is.numeric(contour)){
                contour(x = Linfs, y = Ks, z = t(score_mat), add = TRUE, nlevels = contour)
            } else {
                if(add.values){
                    text(x=plot_dat$Var2,y=plot_dat$Var1,round(as.numeric(plot_dat$value),digits = 2),cex = 0.6)
                }
            }
        }
    } else {
        if(all(Ks %in% exp(seq(log(0.1),log(10),length.out=100)))){
            K_labels <- c(seq(0.1,1,0.1),2:10)
            K_plot <- log10(Ks)
            K_ats <- log10(K_labels)
        } else {
            K_labels <- Ks
            K_plot <- Ks
            K_ats <- Ks
        }
        phis <- round(log10(K_labels) + 2 * log10(Linfs), digits = 2)

        op <- par(mar = c(12,5,4,2))
        plot(K_plot,score_mat,type = 'l', ylim=c(0,max(score_mat, na.rm = TRUE)*1.4),
             ylab = "Score function", xlab = "Growth constant K (/year)", col = "red", lwd=2,xaxt='n')
        axis(1,at = K_ats,labels = K_labels)
        axis(1,at = K_ats,labels = phis,
             line = 5.5)
        mtext(text = expression(paste("Growth performance index (",phi,"')")),side = 1,line = 8.5)
        grid(nx = 0, NULL, lty = 6, col = "gray40")
        abline(v = K_ats, lty = 6, col = "gray40")
        if(plot_title) title("K-Scan", line = 1)
        par(op)
    }

    Rn_max <- max(score_mat, na.rm = TRUE)[1]
    idxs <- which(score_mat == Rn_max, arr.ind = TRUE)[1,]
    Linfest <- as.numeric(as.character(colnames(score_mat)[idxs[2]]))
    Kest <- as.numeric(as.character(rownames(score_mat)[idxs[1]]))
    tanchest <- as.numeric(as.character(ESP_tanch_L[idxs[1],idxs[2]]))

    ## growth performance index
    phiL <- log10(Kest) + 2 * log10(Linfest)

    pars <- list(Linf = Linfest,
                 K = Kest,
                 ta = tanchest,
                 C = C,
                 ts = ts,
                 phiL = phiL)

    final_res <- lfqFitCurves(lfq = res, par=pars,
                              flagging.out = flagging.out,
                              agemax = agemax)

    ## Results
    res$score_mat <- score_mat
    res$ncohort <- final_res$ncohort
    res$agemax <- final_res$agemax
    res$par <- pars
    res$fESP <- Rn_max
    res$Rn_max <- Rn_max


    class(res) <- "lfq"
    if(plot){
        plot(res, Fname = "rcounts")
        Lt <- lfqFitCurves(res, par = pars, draw=TRUE)
    }
    return(res)
}

#' @title ELEFAN_SA
#'
#' @description Electronic LEngth Frequency ANalysis with simulated annealing
#'    for estimating growth parameters.
#'
#' @param x a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row)
#'      and sampling date (column);
#' }
#' @param seasonalised logical; indicating if the seasonalised von Bertalanffy
#'    growth function should be applied (default: FALSE).
#' @param init_par a list providing the Initial values for the components to be
#' optimized. When set to NULL the following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is the maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.5),
#'   \item \strong{ta} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0.5),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param low_par a list providing the lower bounds for components. When set to
#' NULL the following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.01),
#'   \item \strong{ta} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param up_par a list providing the upper bounds for components. When set to
#' NULL the following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.01),
#'   \item \strong{ta} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param SA_time numeric; Maximum running time in seconds (default : 60 * 1).
#' @param maxit Integer. Maximum number of iterations of the
#'              algorithm. Default is NULL.
#' @param nb.stop.improvement Integer. The program will stop when
#'              there is no any improvement in 'nb.stop.improvement'
#'              steps. Default is NULL
#' @param SA_temp numeric; Initial value for temperature (default : 1e5).
#' @param verbose logical; TRUE means that messages from the algorithm
#'    are shown (default : TRUE).
#' @param MA number indicating over how many length classes the moving average
#'  should be performed (defalut: 5, for
#'    more information see \link{lfqRestructure}).
#' @param addl.sqrt Passed to \link{lfqRestructure}. Applied an additional
#'    square-root transformation of positive values according to Brey et al. (1988).
#'    (default: FALSE, for more information see \link{lfqRestructure}).
#' @param agemax maximum age of species; default NULL, then estimated from Linf
#' @param flagging.out logical; passed to \link{lfqFitCurves}. Default is TRUE
#' @param plot logical; Plot restructured counts with fitted lines using
#' \code{\link{plot.lfq}} and \code{\link{lfqFitCurves}} (default : FALSE).
#' @param plot.score logical; Plot simulated annealing score progression.
#'    (Default: plot.score=TRUE)
#'
#' @examples
#' \donttest{
#' ## synthetic lfq data example
#' data(synLFQ4)
#' plot(synLFQ4, Fname="catch")
#'
#' # ELEFAN_SA (takes approximately 2 minutes)
#' output <- ELEFAN_SA(synLFQ4, SA_time = 60*2, seasonalised = TRUE, MA = 11,
#'   init_par = list(Linf = 75, K = 0.5, ta = 0.5, C = 0.5, ts = 0.5),
#'   low_par = list(Linf = 70, K = 0.3, ta = 0, C = 0, ts = 0),
#'   up_par = list(Linf = 90, K = 0.7, ta = 1, C = 1, ts = 1))
#' output$par
#' output$Rn_max
#'
#' # view fit
#' plot(output)
#'
#' # or
#' plot(output, draw = FALSE)
#' lfqFitCurves(output, col=1, par=output$par, draw=TRUE)$ESP
#'
#' # compare to original parameters
#' tmp <- lfqFitCurves(output, col=4, lty=1,
#'    par=list(Linf=80, K=0.5, ta=0.25, C=0.75, ts=0.5), draw=TRUE)
#' tmp$fESP
#' output$Rn_max
#' }
#'
#' @details A more detailed description of the simulated annealing (SA) can be found in
#'    Xiang et al. (2013). The score value \code{cost_value} is not comparable with
#'    the score value of the other ELEFAN functions (\code{\link{ELEFAN}} or
#'    \code{\link{ELEFAN_GA}}).
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{rcounts}: restructured frequencies,
#'   \item \strong{peaks_mat}: matrix with positive peaks with distinct values,
#'   \item \strong{ASP}: available sum of peaks, sum of posititve peaks which
#'      could be potential be hit by
#'      growth curves,
#'   \item \strong{ncohort}: maximum age of species,
#'   \item \strong{agemax}: maximum age of species,
#'   \item \strong{par}: a list with the parameters of the von Bertalanffy growth
#'      function:
#'      \itemize{
#'        \item \strong{Linf}: length infinity in cm,
#'        \item \strong{K}: curving coefficient;
#'        \item \strong{ta}: time point anchoring growth curves in year-length
#'          coordinate system, corrsponds to peak spawning month,
#'        \item \strong{C}: amplitude of growth oscillation
#'          (if \code{seasonalised} = TRUE),
#'        \item \strong{ts}: summer point of oscillation (ts = WP - 0.5)
#'          (if \code{seasonalised} = TRUE),
#'        \item \strong{phiL}: growth performance index defined as
#'          phiL = log10(K) + 2 * log10(Linf);
#'      }
#'   \item \strong{Rn_max}:  highest score value (absolute value of cost function,
#'   comparable with ELEFAN and ELEFAN_GA).
#' }
#'
#' @importFrom graphics par plot title lines grid
#' @importFrom grDevices adjustcolor
#' @importFrom stats median
#' @importFrom GenSA GenSA
#' @importFrom utils flush.console
#'
#' @references
#' Brey, T., Soriano, M., and Pauly, D. 1988. Electronic length frequency analysis: a revised and
#' expanded
#' user's guide to ELEFAN 0, 1 and 2.
#'
#' Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the objective extraction of
#' growth parameters from length-frequency data. \emph{Meeresforschung}, 28(4):205-211
#'
#' Xiang, Y., Gubian, S., Suomela, B., & Hoeng, J. (2013). Generalized simulated
#' annealing for global optimization: the GenSA Package. R Journal, 5(1), 13-28.
#' @export

ELEFAN_SA <- function(lfq,
                      seasonalised = FALSE,
                      init_par = list(Linf = 50, K = 0.5, ta = 0.5, C = 0, ts = 0),
                      low_par = NULL,
                      up_par = NULL,
                      SA_time = 60 * 1,
                      maxit = NULL,
                      nb.stop.improvement = NULL,
                      SA_temp = 1e5,
                      verbose = TRUE,
                      MA = 5, addl.sqrt = FALSE,
                      agemax = NULL,
                      flagging.out = TRUE,
                      plot = FALSE,
                      plot.score = TRUE
                      ){

    res <- lfq
    classes <- res$midLengths
    n_classes <- length(classes)
    Linf_est <- classes[n_classes]


                                        # initial parameters
    init_par_ALL <- list(Linf = Linf_est,
                         K = 0.5, ta = 0.5, C = 0, ts = 0)
    init_Linf <- ifelse("Linf" %in% names(init_par),
                        get("Linf", init_par),
                        get("Linf", init_par_ALL))
    init_K <- ifelse("K" %in% names(init_par),
                     get("K", init_par),
                     get("K", init_par_ALL))
    init_tanc <- ifelse("ta" %in% names(init_par),
                        get("ta", init_par),
                        get("ta", init_par_ALL))
    init_C <- ifelse("C" %in% names(init_par),
                     get("C", init_par),
                     get("C", init_par_ALL))
    init_ts <- ifelse("ts" %in% names(init_par),
                      get("ts", init_par),
                      get("ts", init_par_ALL))

                                        # lower parameter bounds
    low_par_ALL <- list(Linf = Linf_est * 0.5,
                        K = 0.01,
                        ta = 0,
                        C = 0,
                        ts = 0)
    low_Linf <- ifelse("Linf" %in% names(low_par),
                       get("Linf", low_par),
                       get("Linf", low_par_ALL))
    low_K <- ifelse("K" %in% names(low_par),
                    get("K", low_par),
                    get("K", low_par_ALL))
    low_tanc <- ifelse("ta" %in% names(low_par),
                       get("ta", low_par),
                       get("ta", low_par_ALL))
    low_C <- ifelse("C" %in% names(low_par),
                    get("C", low_par),
                    get("C", low_par_ALL))
    low_ts <- ifelse("ts" %in% names(low_par),
                     get("ts", low_par),
                     get("ts", low_par_ALL))

                                        # upper parameter bounds
    up_par_ALL <- list(Linf = Linf_est * 1.5,
                       K = 1,
                       ta = 1,
                       C = 1,
                       ts = 1)
    up_Linf <- ifelse("Linf" %in% names(up_par),
                      get("Linf", up_par),
                      get("Linf", up_par_ALL))
    up_K <- ifelse("K" %in% names(up_par),
                   get("K", up_par),
                   get("K", up_par_ALL))
    up_tanc <- ifelse("ta" %in% names(up_par),
                      get("ta", up_par),
                      get("ta", up_par_ALL))
    up_C <- ifelse("C" %in% names(up_par),
                   get("C", up_par),
                   get("C", up_par_ALL))
    up_ts <- ifelse("ts" %in% names(up_par),
                    get("ts", up_par),
                    get("ts", up_par_ALL))


                                        # ELEFAN 0
    res <- lfqRestructure(res, MA = MA, addl.sqrt = addl.sqrt)
    catch_aAF_F <- res$rcounts
    peaks_mat <- res$peaks_mat
    ASP <- res$ASP

                                        # seasonalised cost function
    soSAfun <- function(lfq, par=c(init_Linf, init_K, init_tanc, init_C, init_ts),
                        agemax, flagging.out){
        Lt <- lfqFitCurves(lfq,
                           par=list(Linf=par[1], K=par[2], ta=par[3], C=par[4], ts=par[5]),
                           flagging.out = flagging.out, agemax = agemax)
        return(-Lt$fESP)
    }
                                        # cost function
    SAfun <- function(lfq, par=c(init_Linf, init_K, init_tanc),
                      agemax, flagging.out){
        Lt <- lfqFitCurves(lfq,
                           par=list(Linf=par[1], K=par[2], ta=par[3], C = 0, ts = 0),
                           flagging.out = flagging.out, agemax = agemax)
        return(-Lt$fESP)
    }


    ## control list
    control <- list(temperature = SA_temp,
                    verbose = verbose)
    if(!is.null(SA_time)) control$max.time = SA_time
    if(!is.null(maxit)) control$maxit = maxit
    if(!is.null(nb.stop.improvement)) control$nb.stop.improvement


    if(seasonalised){
                                        # Simulated annealing with seasonalised VBGF
        writeLines(paste(
            "Simulated annealing is running. \nThis will take approximately",
            round(SA_time/60,digits=2),
            "minutes.",
            sep=" "
        ))
        flush.console()
        SAfit <- GenSA::GenSA(
                            par = c(init_Linf, init_K, init_tanc, init_C, init_ts),
                            fn = soSAfun,
                            lower = c(low_Linf, low_K, low_tanc, low_C, low_ts),
                            upper = c(up_Linf, up_K, up_tanc, up_C, up_ts),
                            agemax = agemax,
                            flagging.out = flagging.out,
                            control = control,
                            lfq = res
                        )

        pars <- as.list(SAfit$par)
        names(pars) <- c("Linf","K","ta", "C", "ts")
    }else{
                                        # Simulated annealing
        writeLines(paste(
            "Simulated annealing is running. \nThis will take approximately",
            round(SA_time/60,digits=2),
            "minutes.",
            sep=" "
        ))
        flush.console()
        SAfit <- GenSA::GenSA(par = c(init_Linf, init_K, init_tanc),
                              fn = SAfun,
                              lower = c(low_Linf, low_K, low_tanc),
                              upper = c(up_Linf, up_K, up_tanc),
                              agemax = agemax,
                              flagging.out = flagging.out,
                              control = control,
                              lfq = res
                              )

        pars <- as.list(SAfit$par)
        names(pars) <- c("Linf","K","ta")
    }

                                        # Score graph in GA style
    tmp <- as.data.frame(SAfit$trace.mat)
    meani <- aggregate(tmp$function.value, list(step = tmp$nb.steps),mean, na.rm = TRUE)
    exe <- aggregate(tmp$current.minimum, list(step = tmp$nb.steps),mean, na.rm = TRUE)
    medi <- aggregate(tmp$function.value, list(step = tmp$nb.steps),median, na.rm = TRUE)
    ylim <- c(min(range(exe$x,na.rm = TRUE, finite = TRUE)),
              max(range(meani$x, na.rm = TRUE, finite = TRUE)))
    if(plot.score){
        op <- par(mar=c(5.1, 4.1, 1, 4.1))
        plot(tmp$nb.steps, tmp$function.value, type = "n", ylim = ylim, xlab = "Iteration",
             ylab = "Cost value")
        graphics::grid(equilogs = FALSE)
        points(tmp$nb.steps, tmp$current.minimum, type = "o", pch = 16, lty = 1,
               col = "green3", cex = 0.7)
        points(meani$step, meani$x, type = "o", pch = 1, lty = 2,
               col = "dodgerblue3", cex = 0.7)
        polygon(c(meani$step, rev(meani$step)),
                c(exe$x, rev(medi$x)),
                border = FALSE, col = adjustcolor("green3", alpha.f = 0.1))
        par(new=TRUE)
        plot(tmp$nb.steps, tmp$temperature, t="l", col=2, lty=2, log="y", axes = FALSE, xlab = "", ylab = "")
        axis(4, col=2, col.axis=2); mtext(text = "Temperature", side = 4, line = par()$mgp[1], col=2)
        legend("topright", legend = c("Best", "Mean", "Median", "Temperature"),
               col = c("green3", "dodgerblue3", adjustcolor("green3", alpha.f = 0.1), 2),
               pch = c(16, 1, NA, NA), lty = c(1,2,1,2),
               lwd = c(1, 1, 10, 1), pt.cex = c(rep(0.7,2), 2, NA),
               inset = 0.02)
        par(op)
    }

    final_res <- lfqFitCurves(lfq = res,par=pars,flagging.out = flagging.out,
                              agemax = agemax)

                                        # growth performance index
    phiL <- log10(pars$K) + 2 * log10(pars$Linf)
    pars$phiL <- phiL

                                        # Results
    res$ncohort <- final_res$ncohort
    res$agemax <- final_res$agemax
    res$par <- pars
    res$fESP <- abs(SAfit$value)
    res$Rn_max <- abs(SAfit$value)

    if(plot){
        plot(res, Fname = "rcounts")
        Lt <- lfqFitCurves(res, par = res$pars, draw=TRUE)
    }
    return(res)
}
#' @title ELEFAN_GA
#'
#' @description Electronic LEngth Frequency ANalysis with genetic algorithm
#' used for estimating growth parameters.
#'
#' @param x a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row)
#'      and sampling date (column);
#' }
#' @param seasonalised logical; indicating if the seasonalised von Bertalanffy
#'    growth function should be applied (default: FALSE).
#' @param low_par a list providing the minimum of the search space in case
#' of real-valued or permutation encoded optimizations. When set to NULL the
#' following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.01),
#'   \item \strong{ta} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param up_par a list providing the maximum of the search space in case of
#' real-valued or permutation encoded optimizations. When set to NULL the
#' following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 1),
#'   \item \strong{ta} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 1),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 1),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 1);
#' }
#' @param popSize the population size. Default: 50
#' @param maxiter the maximum number of iterations to run before the
#' GA search is halted. default:100
#' @param run the number of consecutive generations without any improvement
#' in the best fitness value before the GA is stopped. Default: equals maxiter
#' @param parallel a logical argument specifying if parallel computing
#' should be used (TRUE) or not (FALSE, default) for evaluating the
#' fitness function. See \code{\link[GA]{ga}} for details. Default:FALSE, but
#' setting to TRUE may substantially improve required calculation time. Use of
#' this functionality requires the following packages: parallel, doParallel.
#' @param pmutation the probability of mutation in a parent chromosome.
#' Usually mutation occurs with a small probability, and by default is set to 0.1.
#' @param pcrossover the probability of crossover between pairs of chromosomes.
#' Typically this is a large value and by default is set to 0.8.
#' @param elitism the number of best fitness individuals to survive at each generation.
#' By default the top 5\% individuals will survive at each iteration.
#' @param MA number indicating over how many length classes the moving average
#' should be performed (default: 5, for
#'    more information see \link{lfqRestructure})
#' @param addl.sqrt additional squareroot transformation of positive values
#' according to Brey et al. (1988) (default: FALSE, for
#'    more information see \link{lfqRestructure})
#' @param agemax maximum age of species; default NULL, then estimated from Linf
#' @param flagging.out logical; should positive peaks be flagged out? Original setting of
#' ELEFAN in TRUE. Default:TRUE
#' @param seed an integer value containing the random number generator state. This
#' argument can be used to replicate the results of a GA search. Note that
#' if parallel computing is required, the doRNG package must be installed.
#' (Default: 'seed = NULL')
#' @param monitor a logical or an R function which takes as input the current
#'                state of the 'ga-class' object and show the evolution of the
#'                search. By default, 'monitor = FALSE' so any
#'                output is suppressed. Possible also, the functions
#'                'gaMonitor' or 'gaMonitor2' (depending on whether or not is
#'                an RStudio session) which print the average and best fitness
#'                values at each iteration. If set to 'plot' these information
#'                are plotted on a graphical device. Other functions can be
#'                written by the user and supplied as argument.
#' @param plot logical; Plot restructured counts with fitted lines using
#' \code{\link{plot.lfq}} and \code{\link{lfqFitCurves}} (default : FALSE).
#' @param plot.score logical; Plot genetic algorithm fitness progression.
#'    (Default: plot.score=TRUE).
#' @param ... additional parameters to pass to \code{\link[GA]{ga}}
#'
#'
#' @details A more detailed description of the generic algorithm (GA) can be found in
#'    Scrucca (2013). The score value \code{fitnessValue} is not comparable with
#'    the score value of the other ELEFAN functions (\code{\link{ELEFAN}} or
#'    \code{\link{ELEFAN_SA}}).
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{samplingPeriod}: length of sampling period in years,
#'   \item \strong{samplingDays}: time of sampling times in relation to first sampling time,
#'   \item \strong{delta_t}: array with time differences between relative sampling time set to zero and
#'      other sampling times,
#'   \item \strong{rcounts}: restructured frequencies,
#'   \item \strong{peaks_mat}: matrix with positive peaks with distinct values,
#'   \item \strong{ASP}: available sum of peaks, sum of posititve peaks which could be
#'      potential be hit by
#'      growth curves,
#'   \item \strong{ncohort}: maximum age of species,
#'   \item \strong{agemax}: maximum age of species,
#'   \item \strong{par}: a list with the parameters of the von Bertalanffy growth
#'      function:
#'      \itemize{
#'        \item \strong{Linf}: length infinity in cm,
#'        \item \strong{K}: curving coefficient;
#'        \item \strong{ta}: time point anchoring growth curves in year-length
#'          coordinate system, corrsponds to peak spawning month,
#'        \item \strong{C}: amplitude of growth oscillation
#'          (if \code{seasonalised} = TRUE),
#'        \item \strong{ts}: summer point of oscillation (ts = WP - 0.5)
#'          (if \code{seasonalised} = TRUE),
#'        \item \strong{phiL}: growth performance index defined as
#'          phiL = log10(K) + 2 * log10(Linf);
#'      }
#'   \item \strong{Rn_max}: highest value of fitness function, (comparable with ELEFAN and ELEFAN_SA).
#' }
#'
#' @examples
#' \donttest{
#' # load data and view catch length frequencies
#' data(synLFQ4)
#' plot(synLFQ4, Fname="catch")
#'
#' # Genetic algorithm
#' # (if using a multicore processor,
#' #   consider adding the argument 'parallel=TRUE'
#' #   to reduce computation time)
#' output <- ELEFAN_GA(synLFQ4, seasonalised = TRUE,
#'    low_par = list(Linf = 70, K = 0.25, ta = 0, C = 0, ts= 0),
#'    up_par = list(Linf = 90, K = 0.7, ta = 1, C = 1, ts = 1),
#'    popSize = 40, maxiter = 50, run = 20,
#'    MA = 11, plot = TRUE, seed = 1111)
#' output$par
#' output$ASP
#' output$Rn_max
#'
#' # compare fitness score (fESP) to
#' # that calculated with "true" growth parameter values
#' plot(output, draw = FALSE)
#' lfqFitCurves(output, par=list(Linf=80, K=0.5, ta=0.25, C=0.75, ts=0.5),
#'        draw = TRUE, col=1, flagging.out = FALSE)$fESP
#' lfqFitCurves(output, par=output$par, draw = TRUE, col=2, flagging.out = FALSE)$fESP
#' legend("top", legend=c("orig.", "GA"), lty=2, col=1:2, ncol=2)
#'}
#'
#' @import parallel
#' @import doParallel
#' @importFrom GA ga
#' @importFrom utils flush.console
#'
#' @references
#' Brey, T., Soriano, M., and Pauly, D. 1988. Electronic length frequency analysis: a
#'    revised and expanded
#' user's guide to ELEFAN 0, 1 and 2.
#'
#' Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the objective extraction of
#' growth parameters from length-frequency data. \emph{Meeresforschung}, 28(4):205-211
#'
#' Scrucca, L. (2013). GA: a package for genetic algorithms in R. Journal of
#' Statistical Software, 53(4), 1-37.
#'
#' @export

ELEFAN_GA <- function(lfq,
                      seasonalised = FALSE,
                      low_par = NULL,
                      up_par = NULL,
                      popSize = 50,
                      maxiter = 100,
                      run = maxiter,
                      parallel = FALSE,
                      pmutation = 0.1,
                      pcrossover = 0.8,
                      elitism = base::max(1, round(popSize*0.05)),
                      MA = 5,
                      addl.sqrt = FALSE,
                      agemax = NULL,
                      flagging.out = TRUE,
                      seed = NULL,
                      monitor = FALSE,
                      plot = FALSE,
                      plot.score = TRUE,
                      ...
                      ){

    classes <- lfq$midLengths
    n_classes <- length(classes)
    Linf_est <- classes[n_classes]

    ## lower parameter bounds
    low_par_ALL <- list(Linf = Linf_est * 0.5,
                        K = 0.01,
                        ta = 0,
                        C = 0,
                        ts = 0)
    low_Linf <- ifelse("Linf" %in% names(low_par),
                       get("Linf", low_par),
                       get("Linf", low_par_ALL))
    low_K <- ifelse("K" %in% names(low_par),
                    get("K", low_par),
                    get("K", low_par_ALL))
    low_tanc <- ifelse("ta" %in% names(low_par),
                       get("ta", low_par),
                       get("ta", low_par_ALL))
    low_C <- ifelse("C" %in% names(low_par),
                    get("C", low_par),
                    get("C", low_par_ALL))
    low_ts <- ifelse("ts" %in% names(low_par),
                     get("ts", low_par),
                     get("ts", low_par_ALL))

    ## upper parameter bounds
    up_par_ALL <- list(Linf = Linf_est * 1.5,
                       K = 1,
                       ta = 1,
                       C = 1,
                       ts = 1)
    up_Linf <- ifelse("Linf" %in% names(up_par),
                      get("Linf", up_par),
                      get("Linf", up_par_ALL))
    up_K <- ifelse("K" %in% names(up_par),
                   get("K", up_par),
                   get("K", up_par_ALL))
    up_tanc <- ifelse("ta" %in% names(up_par),
                      get("ta", up_par),
                      get("ta", up_par_ALL))
    up_C <- ifelse("C" %in% names(up_par),
                   get("C", up_par),
                   get("C", up_par_ALL))
    up_ts <- ifelse("ts" %in% names(up_par),
                    get("ts", up_par),
                    get("ts", up_par_ALL))

    ## ELEFAN 0
    lfq <- lfqRestructure(lfq, MA = MA, addl.sqrt = addl.sqrt)

    ## seasonalised fitness function
    sofun <- function(lfq, par, agemax, flagging.out){
        Lt <- lfqFitCurves(lfq,
                           par=list(Linf=par[1], K=par[2], ta=par[3], C=par[4], ts=par[5]),
                           agemax = agemax, flagging.out = flagging.out)
        return(Lt$fESP)
    }
    ## non-seasonalised fitness function
    fun <- function(lfq, par, agemax, flagging.out){
        Lt <- lfqFitCurves(lfq,
                           par=list(Linf=par[1], K=par[2], ta=par[3], C = 0, ts = 0),
                           agemax = agemax, flagging.out = flagging.out)
        return(Lt$fESP)
    }


    ## Genetic algorithm
    if(seasonalised){
        min = c(low_Linf, low_K, low_tanc, low_C, low_ts)
        max = c(up_Linf, up_K, up_tanc, up_C, up_ts)
        writeLines("Genetic algorithm is running. This might take some time.")
        flush.console()
        fit <- GA::ga(
                       type = "real-valued",
                       fitness = sofun, lfq=lfq,
                       lower = min,
                       upper = max,
                       agemax = agemax,
                       flagging.out = flagging.out,
                       popSize = popSize, maxiter = maxiter, run = run, parallel = parallel,
                       pmutation = pmutation, pcrossover = pcrossover, elitism = elitism,
                       seed = seed, monitor = monitor,
                       ...
                   )
        pars <- as.list(fit@solution[1,])
        names(pars) <- c("Linf", "K", "ta", "C", "ts")
    }else{
        min = c(low_Linf, low_K, low_tanc)
        max = c(up_Linf, up_K, up_tanc)
        writeLines("Genetic algorithm is running. This might take some time.")
        flush.console()
        fit <- GA::ga(
                       type = "real-valued",
                       fitness = fun,
                       lfq=lfq,
                       lower = min,
                       upper = max,
                       agemax = agemax,
                       flagging.out = flagging.out,
                       popSize = popSize, maxiter = maxiter, run = run, parallel = parallel,
                       pmutation = pmutation, pcrossover = pcrossover, elitism = elitism,
                       seed = seed,
                       monitor = monitor,
                       ...
                   )
        pars <- as.list(fit@solution[1,])
        names(pars) <- c("Linf", "K", "ta")
    }

    ## Fitness graph
    if(plot.score){
        GA::plot(fit)
    }

    final_res <- lfqFitCurves(
        lfq = lfq, par=pars,
        flagging.out = flagging.out,
        agemax = agemax)

    ## growth performance index
    phiL <- log10(pars$K) + 2 * log10(pars$Linf)
    pars$phiL <- phiL

    ## Results
    lfq$ncohort <- final_res$ncohort
    lfq$agemax <- final_res$agemax
    lfq$par <- pars
    lfq$fESP <- fit@fitnessValue
    lfq$Rn_max <- fit@fitnessValue

    if(plot){
        plot(lfq, Fname = "rcounts")
        Lt <- lfqFitCurves(lfq, par = lfq$pars, draw=TRUE)
    }
    return(lfq)
}
#' @title Growth from tagging data
#'
#' @description  This function estimates growth parameters from tagging data. Munro plot
#'    is applied
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{L1}: length at tagging [cm],
#'   \item \strong{L2}: length at recapture [cm],
#'   \item \strong{delta_t}: time interval between tagging ang recapture
#'   (instead two vectors with \strong{t1} (age at tagging) and \strong{t2}
#'   (age at recapture) can be provided.
#' }
#' @param method indicating which of following methods should be applied: "GullandHolt"
#'    or "Munro".
#' @param Linf_range two values indicating the lower and upper limits of the range,
#'    in which the \link{optimise} searches for the Linf value with the best fit
#'    (lowest CV value ),
#' @param time_unit indicating the unit of the time interval, either "year", "month",
#'    "week", or "day"
#'
#' @examples
#' # from Wolff (1984)
#' dat <- list(L1 = c(40,46,29,30,18,31,48,49,59,58,61,65,57,55),
#'    L2 = c(85,53,55,56,25,43,70,59,62,80,72,83,65,56),
#'    delta_t = c(289,26,84,77,14,38,89,38,28,149,89,74,38,21))
#' growth_tagging(param = dat, "Munro", time_unit = "day", Linf_range=c(80,120))
#' growth_tagging(param = dat, "GullandHolt", time_unit = "day")
#'
#'
#' # from Sparre and Venema (1999)
#' dat <- list(L1 = c(9.7,10.5,10.9,11.1,12.4,12.8,14.0,16.1,16.3,17.0,17.7),
#'    L2 = c(10.2,10.9,11.8,12.0,15.5,13.6,14.3,16.4,16.5,17.2,18.0),
#'    delta_t = c(53,33,108,102,272,48,53,73,63,106,111))
#' growth_tagging(param = dat, "Munro", time_unit = "day", Linf_range = c(10,40))
#' growth_tagging(param = dat, "GullandHolt", time_unit = "day")
#'
#' @details
#' If Munro plot is applied the optimal Linf value is found by minimizing the coefficient
#' of variation (CV = sd(K)/mean(K)). For this iterative method the \link{optimise}
#' function is applied. The histogram of the individual K values allows to distinguish
#' potential differences in growth performance between individuals. t0 can not be
#' estimated by Munro plot, neither by the Gulland Holt method.
#'
#' @return A list with the input parameters and following parameters:
#' \itemize{
#'  \item \strong{x}: independent variable used for regression analysis,
#'  \item \strong{y}: dependent variable used for regression analysis,
#'  \item \strong{reg_coeffs}: regression coefficients,
#'  \item \strong{r2}: r squared of regression analysis,
#'  \item \strong{Linf}: infinite length for investigated species in cm [cm],
#'  \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'  \item \strong{conf_int_K}: confidence intervals of K (only if Gulland Holt method
#'      was applied).
#' }
#'
#' @importFrom graphics abline hist plot
#' @importFrom stats lm optimise qt sd
#'
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' Sparre, P., Venema, S.C., 1999. Introduction to tropical fish stock assessment.
#' Part 2. Excercises. FAO Fisheries Technical Paper, (306.2, Rev. 2). 94 p.
#'
#' Wolff, M., 1984. Early setback for scallop culture in Peru.
#'
#' @export

growth_tagging <- function(param, method,
                            Linf_range = c(5,600),time_unit = "year"){

  res <- param

  L1 <- res$L1
  L2 <- res$L2
  if("delta_t" %in% names(res)) delta_t <- res$delta_t
  if(!"delta_t" %in% names(res)) delta_t <- res$t2 - res$t1

  switch(time_unit,
         "year" ={
           time_u <- 1
         },
         "month"={
           time_u <- 12
         },
         "week"={
           time_u <- 52
         },
         "day"={
           time_u <- 365
         })

  # growth increment per time / growth rate
  incr_time <- ((L2-L1) / delta_t ) * time_u

  switch(method,
         "GullandHolt"={
           y <- incr_time

           #mean Lt
           x <- (L1 + L2) / 2

           sum_mod <- summary(lm(y ~ x))
           coeffs <- sum_mod$coefficients
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           K <- -b
           Linf <- a/K
           R2 <- sum_mod$r.squared
           # standard error and confidence limits
           sb2 <- (1/(length(L1)-2)) * ((sd(y,na.rm = TRUE)/sd(x,na.rm = TRUE))^2 - b^2)
           sb <- sqrt(sb2)
           SE_b <- abs(b) * qt(0.975,sum_mod$df[2])
           tg <- qt(0.975,sum_mod$df[2])
           conf_int <- c(K - tg * sb,K + tg * sb)

           plot(y ~ x, pch = 16,
                ylim = c(0,max(y,na.rm=TRUE)*1.1),
                xlim = c(min(x,na.rm=TRUE)*0.9,Linf*1.05),
                xlab = "mean L(t)",
                ylab=expression(paste(delta,"L / ",delta,"t")),
                main = "Gulland and Holt")
           abline(a,b)

           tmp <- list(y = y, x = x, coeffs = coeffs, R2 = R2, conf_int = conf_int)
         },

         "Munro"={
           func <- function(Linf){
             K <- ((log(Linf-L1) - log(Linf-L2)) / (delta_t)) * time_u
             CV <- sd(K, na.rm = TRUE) / mean(K, na.rm = TRUE)
           }

           opt_mod <- optimise(f = func, interval = Linf_range)
           Linf <- opt_mod$minimum
           K <- opt_mod$objective
           indi_Ks <- ((log(Linf-L1) - log(Linf-L2)) / (delta_t)) * time_u
           hist(indi_Ks, breaks = length(L1),
                xlab = "K",
                main = "Histogram  of individual Ks")

           tmp <- list(indi_Ks = indi_Ks)
         })

  res2 <- res
  if(exists("x", where = tmp)) res2$x <- x
  if(exists("y", where = tmp)) res2$y <- y
  if(exists("coeffs", where = tmp)) res2$reg_coeffs <- coeffs
  if(exists("R2", where = tmp)) res2$r2 <- R2
  if(exists("indi_Ks", where = tmp)) res2$indi_Ks <- indi_Ks
  ret <- c(res2,list(Linf = Linf, K = K))
  if(exists("conf_int", where = tmp)) ret$conf_int_K <- conf_int
  return(ret)
}
#' @title Estimation of growth parameter using length-at-age data
#'
#' @description  This function estimates growth parameters from
#'    length-at-age data. It
#'    allows to perform different methods: Gulland and Holt, Ford Walford plot,
#'    Chapman's method, Bertalanffy plot, or non linear least squares method (LSM).
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{age}: age measurements,
#'   \item \strong{length}: corresponding lengths in cm.
#' }
#' @param method indicating which of following methods should be applied:
#'    \code{"GullandHolt"},
#'    \code{"FordWalford"}, \code{"Chapman"}, \code{"BertalanffyPlot"},
#'    or \code{"LSM"}
#' @param Linf_est BertalanffyPlot requires an estimate for Linf to derive K and t0
#'    (for more information see Details).
#' @param Linf_init initital parameter of Linf for non-linear sqaures fitting (default 10)
#' @param K_init initital parameter of K for non-linear sqaures fitting (default 0.1)
#' @param t0_init initital parameter of t0 for non-linear sqaures fitting (default 0)
#' @param CI logical; Should confidence intervals be calculated? This option only
#'    works for the LSM method. Default is FALSE.
#' @param ci.level required confidence level (for LSM method only)
#' @param age_plot sequence with ages used for plotting (LSM method only). By default
#'      age_plot = seq(min(param$age),max(param$age),0.1)
#' @param do.sim logical. Should Monte Carlo simulation be applied? Default = FALSE
#' @param nsim the number of Monte Carlo simulations to be performed,
#'    minimum is 10000 (default).
#'
#' @examples
#' # synthetical length at age data
#' dat <- list(age = rep(1:7,each = 5),
#'    length = c(rnorm(5,25.7,0.9),rnorm(5,36,1.2),rnorm(5,42.9,1.5),rnorm(5,47.5,2),
#'    rnorm(5,50.7,0.4),rnorm(5,52.8,0.5),rnorm(5,54.2,0.7)))
#' growth_length_age(dat, method = "GullandHolt")
#'
#' # Bertalaffy plot
#' growth_length_age(dat, method = "BertalanffyPlot", Linf_est = 50)
#'
#' # non linear least squares method
#' output <- growth_length_age(param = dat, method = "LSM",
#'      Linf_init = 30, CI = TRUE, age_plot=NULL)
#' summary(output$mod)
#'
#' @details
#' Gulland and Holt plot assumes
#' infinitestimal delta t (only reasonable approximation of growth parameters if delta t
#' is small). Ford Walford plot and Chapman assume constant time intervals between ages
#' (delta t). The Bertalanffy plot is a robust method, however it requires an estimate of Linf. As
#' long as this estimate is reasonable the resulting estimate of K is reasonable. For
#' a first estimate of Linf the Powell Wetherall method \link{powell_wetherall} can
#' be used. Otherwise, the largest fish or the average of the ten largest fish can be
#' used for a small or large sample, respectively. All lengths have to be smaller than
#' Linf as otherwise the logarithm is not defined. Oldest fish (if larger than Linf) have
#' to be omitted. Non-linear least squares fitting is the preferred method to estimate
#' growth parameters according to Sparre and Venema (1998). If \code{CI = TRUE} the
#' confidence interval of parameters is calculated and plotted. For plotting the
#' confidence interval the \code{\link{predictNLS}} from the \link{propagate} package
#' is applied.
#'
#' @return A list with the input parameters and following parameters:
#' \itemize{
#'  \item \strong{x}: independent variable used for regression analysis,
#'  \item \strong{y}: dependent variable used for regression analysis,
#'  \item \strong{mod}: (non) linear model,
#'  \item \strong{Linf}: infinite length for investigated species in cm [cm],
#'  \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'  \item \strong{t0}: theoretical time zero, at which individuals of this species hatch
#'      (only for Bertalanffy plot and LSM method).
#'  \item \strong{estimates}: dataframe with growth parameters and confidence intervals
#'      (only if LSM method was applied).
#' }
#'
#' @importFrom propagate predictNLS
#' @importFrom graphics abline lines plot segments
#' @importFrom stats lm nls predict aggregate confint
#' @importFrom grDevices adjustcolor
#' @import MASS
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

growth_length_age <- function(param, method, Linf_est = NA,
                              Linf_init = 10, K_init = 0.1, t0_init = 0,
                              CI = FALSE, ci.level = 0.95,
                              age_plot = NULL,
                              do.sim = FALSE,
                              nsim = 10000){

  res <- param
  t <- res$age
  Lt <- res$length

  switch(method,
         "GullandHolt"={
           # Can not deal with multiple observations (lengths) per age
           if(length(unique(t)) == length(t)){
             delta_t <- diff(t)
           }else{
             delta_t <- diff(unique(t))
             Lt <- aggregate(Lt,list(t = t),mean,na.rm= TRUE)$x
           }

           # growth rate
           y <- rep(NA,length(delta_t))
           for(i in 1:length(delta_t)){
             y[i] <- (Lt[i+1] - Lt[i])/delta_t[i]
           }

           #mean Lt
           x <- rep(NA,(length(Lt)-1))
           for(i in 1:(length(Lt)-1)){
             x[i] <- (Lt[i] + Lt[i+1])/ 2
           }

           mod <- lm(y ~ x)
           sum_mod <- summary(lm(y ~ x))
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           R2 <- sum_mod$r.squared
           K <- -b
           Linf <- -a/b

           plot(y ~ x, pch = 16,
                ylim = c(0,max(y,na.rm=TRUE)*1.1),
                xlab = "mean L(t)",
                ylab=expression(paste(delta,"L / ",delta,"t")),
                main = "Gulland and Holt")
           abline(a,b)

           tmp <- list(delta_t=delta_t,x=x,y=y,R2=R2)

         },
         "FordWalford"={
           # Can not deal with multiple observations (lengths) per age
           if(length(unique(t)) == length(t)){
             delta_t <- diff(t)
           }else{
             delta_t <- diff(unique(t))
             Lt <- aggregate(Lt,list(t = t),mean,na.rm= TRUE)$x
           }

           if(!all(delta_t == 1)) stop("The Ford Walford method assumes constant time intervals!")

           delta_t <- delta_t[1]

           y <- Lt[2:length(Lt)]
           x <- Lt[-length(Lt)]

           mod <- lm(y ~ x)
           sum_mod <- summary(lm(y ~ x))
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           R2 <- sum_mod$r.squared
           K <- -1/delta_t * log(b)
           Linf <- a/(1-b)

           plot(y ~ x, pch = 16, ylim = c(0,max(y,na.rm=TRUE)*1.1),
                xlim = c(0,max(x,na.rm=TRUE)*1.1),
                xlab = "L(t)",
                ylab=expression(paste("L(t+",delta,"t)")),
                main = "Ford Walford plot")
           abline(a,b)
           abline(1,1, lty=2)
           segments(x0 = 0, y0 = Linf, x1 = Linf, y1 = Linf, lty=2)
           segments(x0 = Linf, y0 = 0, x1 = Linf, y1 = Linf, lty=2)

           tmp <- list(delta_t=delta_t,x=x,y=y,R2=R2)
         },
         "Chapman"={
           # Can not deal with multiple observations (lengths) per age
           if(length(unique(t)) == length(t)){
             delta_t <- diff(t)
           }else{
             delta_t <- diff(unique(t))
             Lt <- aggregate(Lt,list(t = t),mean,na.rm= TRUE)$x
           }

           if(!all(delta_t == 1)) stop("The Chapman method assumes constant time intervals!")

           delta_t <- delta_t[1]

           Lt_short <- Lt[-length(Lt)]
           y <- Lt[2:length(Lt)] - Lt_short
           x <- Lt_short


           mod <- lm(y ~ x)
           sum_mod <- summary(lm(y ~ x))
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           R2 <- sum_mod$r.squared
           c <- -b
           K <- -(1/delta_t) * log(1+b)
           Linf <- -a/b

           plot(y ~ x, pch = 16, ylim = c(0,max(y,na.rm=TRUE)*1.1),
                xlim = c(0,max(x,na.rm=TRUE)*1.1),
                xlab = "L(t)",
                ylab=expression(paste("L(t+",delta,"t)-L(t)")),
                main = "Chapman")
           abline(a,b)

           tmp <- list(delta_t=delta_t,x=x,y=y,R2=R2)
         },
         "BertalanffyPlot"={

           delta_t <- diff(t)

           #Linf should be given: otherwise using optim?
           if(is.na(Linf_est)) stop("For the Bertalanffy plot you have to provide an estimate for Linf_est!")
           y <- -log(1 - Lt / Linf_est)
           x <- t

           mod <- lm(y ~ x)
           sum_mod <- summary(lm(y ~ x))
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           R2 <- sum_mod$r.squared
           K <- b
           t0 <- a/-K
           Linf <- Linf_est

           plot(y ~ x, pch = 16, ylim = c(0,max(y,na.rm=TRUE)*1.1),
                xlim = c(0,max(x,na.rm=TRUE)*1.1),
                xlab = "t",
                ylab=expression(paste("-ln(1-L(t)/Linf)")),
                main = "Bertalanffy plot")
           abline(a,b)

           tmp <- list(x=x,y=y,R2=R2,t0=t0)


           if(max(Lt, na.rm=TRUE) > Linf_est) writeLines("Some lengths were larger as the Linf value which was \nprovided. These lengths are omitted.")
         },
         "LSM"={

           nls_mod <- nls(Lt ~ (Linf * (1 - exp(-K * (t - t0)))),
                      start = list(Linf = Linf_init, K = K_init, t0 = t0_init))

           if(!CI){
             plot(Lt ~ t,
                  ylab = "L(t)",
                  xlab = "t(age)",
                  main = "Non-linear least squares method")
             lines(t, predict(nls_mod))
           }
           sum_mod <- summary(nls_mod)
           Linf <- sum_mod$coefficients[1]
           K <- sum_mod$coefficients[2]
           t0 <- sum_mod$coefficients[3]

           mod <- nls_mod
           tmp <- list(t0 = t0)

           if(CI){
             suppressMessages(cis <- confint(nls_mod, level = ci.level))
             nls_res <- data.frame(Names=c("Linf","K","t0"),
                                   Value=c(round(Linf,2),
                                           round(K,2),
                                           round(t0,2)),
                                   Lower.CI = c(round(cis[1,1],2),
                                                round(cis[2,1],2),
                                                round(cis[3,1],2)),
                                   Upper.CI = c(round(cis[1,2],2),
                                                round(cis[2,2],2),
                                                round(cis[3,2],2)))
             tmp$nls_res = nls_res

             # plot with confidence interval
             if(is.null(age_plot)){
               age_plot <- seq(min(floor(t)),max(ceiling(t)),0.1)
             }else age_plot <- age_plot
             # Taylor error propagation and Monte Carlo simulation for confidence interval
             sink(tempfile())
             pred_L <- suppressMessages(propagate::predictNLS(nls_mod,
                                                              do.sim = do.sim,
                                                              nsim = nsim,
                                             newdata = data.frame(t = age_plot)))
             sink()
             # Taylor propagation
             if(!do.sim){
               predVals <- cbind(age_plot,as.data.frame(pred_L$summary[,c(1,5,6)]))
             }
             # Monte Carlo simulation
             if(do.sim){
               predVals <- cbind(age_plot,as.data.frame(pred_L$summary[,c(7,11,12)]))
             }
             names(predVals) <- c("age_plot","fit","lower","upper")
             plot(Lt ~ t, type = "n", ylab = "L(t)", xlab = "t(age)", xlim=c(min(age_plot),max(age_plot)),
                  main = "Non-linear least squares method")
             polygon(c(age_plot,rev(age_plot)),
                     c(predVals$lower,rev(predVals$upper)),
                     border = FALSE, col = adjustcolor("dodgerblue", alpha.f = 0.4))
             points(t, Lt)
             lines(age_plot, predict(nls_mod,newdata = data.frame(t=age_plot)))
           }
         },

         stop(paste("\n",method, "not recognised, possible options are \n",
                    "\"GullandHolt\", \"FordWalford\", \"Chapman\" \n",
                    "\"BertalanffyPlot\", and \"LSM\""))
         )

  res2 <- res
  if(exists("delta_t", where = tmp)) res2$delta_t <- delta_t
  if(exists("x", where = tmp)) res2$x <- x
  if(exists("y", where = tmp)) res2$y <- y
  if(exists("R2", where = tmp)) res2$r2 <- R2

  ret <- c(res2,list(mod = mod,
                    Linf = Linf,
                    K = K))

  if(exists("t0", where = tmp)) ret$t0 <- t0
  if(exists("nls_res", where = tmp)) ret$estimates <- nls_res

  return(ret)
}

