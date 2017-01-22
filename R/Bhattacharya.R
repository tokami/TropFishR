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
#'      Sparre and Venema (1998).
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
