#' @title Bhattacharya's method
#'
#' @description Find relative frequencies and frequency distribution of distinct cohorts
#'      in the observed length frequency distribution by resolving the distribution into
#'      Gaussian components.
#'
#' @param param a list consisting of following parameters:
#'   \code{$midLengths} midpoints of the length class as vector,
#'   \code{catch} a vector with the catch per length class or a matrix with
#'      catches per length class of subsequent years;
#'
#' @keywords Bhattacharya length-frequency
#'
#' @examples
#' \donttest{
#'  data(synLFQ1)
#'  Bhattacharya(synLFQ1)
#' }
#'
#' @details This method includes the \code{\link{identify}} function, which allows to
#'      choose points on a graphical device manually. To stop this process please press
#'      right mouse click on the graph, and in case of windows click on "Stop".
#'
#' @return A list with the input parameters and a data.frame (\strong{bhat_results})
#'      with the results of the Bhattacharya method. Most important is the column
#'      \strong{Nx}, where x is the number of the cohort, displaying the number of
#'      individuals per cohort. For an explanation of the other columns please refer to
#'      the FAO manual (Sparre & Venema, 1998).
#'
#' @references
#' Bhattacharya, C.G., 1967. A simple method of resolution of a distribution into Gaussian
#' components, \emph{Biometrics}, 23:115-135
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export



param <- synLFQ1

Bhattacharya <- function(param){

  res <- param
  if("midLengths" %in% names(res) == TRUE){
    midLengths <- as.character(res$midLengths)
  }else stop(
    noquote('The parameter list does not contain an object with name "midLengths"'))
  #if("catch_mat" %in% names(res) == TRUE) catch <- res$catch_mat
  if("catch" %in% names(res) == TRUE){
    catch <- res$catch
  }else stop(
    noquote('The parameter list does not contain an object with name "catch"'))

  #Transform all input variables to correct class
  midLengths <- as.numeric(midLengths)
  #catch <- as.numeric(as.character(catch))

  #Transform matrix into vector if provided
  if(class(catch) == 'matrix'){
    catch.vec <- rowSums(catch, na.rm = TRUE)
  }else catch.vec <- catch
  if(length(midLengths) != length(catch.vec)) stop(
    noquote("midLengths and catch do not have the same length!"))

  df.Bh <- data.frame(midLengths = midLengths,
                      catch = catch.vec)

  #calculate size class interval
  interval <- df.Bh$midLengths[2] - df.Bh$midLengths[1]

  #STEP 1: fills second column of bhat with counts per size class
  bhat.table <- data.frame('mean.length.classes' = df.Bh$midLengths,
                           'N1.plus' = df.Bh$catch,
                           'log.N1.plus' = NA,
                           'delta.log.N1.plus' = NA,
                           'L' = NA,
                           'delta.log.N' = NA,
                           'log.N1' = NA,
                           'N1' = NA,
                           'N2.plus' = NA)

  bhat.table.list <- vector("list", 31)
  a.b.list <- vector("list", 30)
  colour.xy <- c('blue','darkgreen','red','yellow','purple','orange',
                 'lightgreen','skyblue','brown','darkblue','darkorange','darkred')
  colour.vec <- rep('black', length(bhat.table$L))

  print(noquote("Starting on the left, please choose from black points which lie on a straight line! Do not include points which might be affected by the next distribution!"))
  print(noquote("To stop the process press right click (and choose 'Stop' if necessary)"))

  for(xy in 1:30){  #more than 30 cohorts?

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
#     for(i in 1:length(bhat.table$log.N1.plus)){
#       if(!is.na(bhat.table$delta.log.N1.plus[i]) &
#          bhat.table$delta.log.N1.plus[i] != 'Inf'){
#         min.size <- (bhat.table$mean.length.classes[i] - (interval/2))
#         bhat.table$L[i] <- min.size
#       }
#     }


    #STEP 5: plot of fifth against fourth column
    checki2 <- !is.na(bhat.table$delta.log.N1.plus) &
      bhat.table$delta.log.N1.plus != 'Inf'
    bhat.table$L[checki2] <- (bhat.table$mean.length.classes[checki2] - (interval/2))
#     for(i in 1:length(bhat.table$log.N1.plus)){
#       if(!is.na(bhat.table$delta.log.N1.plus[i]) & bhat.table$delta.log.N1.plus[i] != 'Inf'){
#         min.size <- (bhat.table$mean.length.classes[i] - (interval/2))
#         bhat.table$L[i] <- min.size
#       }
#     }
    #STEP 5: plot of fifth against fourth column
#     dev.new()
#     plot(bhat.table$delta.log.N1.plus ~ bhat.table$L, pch = 16)
#     abline(h=0)
#     p.bhat1 <- recordPlot()

    #STEP 6: select points for regression line
    if(xy == 1){bhat.table.list[[1]] <- bhat.table}

    # do not include in linear regression if number in N1.plus is under 3 or 4 if(bhat.table$N1.plus < 3){}

    #find new maximum for start of next cohort
#     if(xy > 1){
#       gradi <- NA
#       for(nm in 3:(length(bhat.table.list[[1]]$delta.log.N1.plus)-1)){
#         if((bhat.table.list[[1]]$delta.log.N1.plus[nm] >
#             bhat.table.list[[1]]$delta.log.N1.plus[nm-1]) == T &
#            (bhat.table.list[[1]]$delta.log.N1.plus[nm] >
#             bhat.table.list[[1]]$delta.log.N1.plus[nm+1]) == T
#         ){gradi[nm] <- T}else gradi[nm] <- F
#       }
#
#       cut.max <- max.rem + min(which(gradi[(max.rem+1):length(gradi)] == T))
#
#       bhat.table.list[[1]]$delta.log.N1.plus[
#         which(as.numeric(rownames(bhat.table.list[[1]])) >= cut.max)] <-
#         bhat.table$delta.log.N1.plus[
#           which(as.numeric(rownames(bhat.table.list[[1]])) >= cut.max)]
#       }


    par(mfrow=c(2,1))
    if(xy == 1){
      dev.new()
      plot(bhat.table.list[[1]]$L, bhat.table.list[[1]]$delta.log.N1.plus, pch = 16,
           col = colour.vec)
      abline(h = 0)
      text(bhat.table.list[[1]]$L, bhat.table.list[[1]]$delta.log.N1.plus,
           labels = row.names(bhat.table.list[[1]]), cex= 0.7, pos=3)
    }
    if(xy > 1){
      dev.new()
      plot(bhat.table.list[[1]]$L, bhat.table.list[[1]]$delta.log.N1.plus, pch = 4,
           col = colour.vec)
      abline(h = 0)
      text(bhat.table.list[[1]]$L[1:last_choices[2]],
           bhat.table.list[[1]]$delta.log.N1.plus[1:last_choices[2]],
           labels = row.names(bhat.table.list[[1]][(1:last_choices[2]),]),
           cex= 0.7, pos=3)

      abline(a = a.b.list[[xy-1]][1], b=a.b.list[[xy-1]][2], col=colour.xy[xy-1])

      points(bhat.table$L[(last_choices[2]+1):length(bhat.table$L)],
             bhat.table$delta.log.N1.plus[(last_choices[2]+1):length(bhat.table$L)],
             col = 'black', pch = 16)
      text(bhat.table$L[(last_choices[2]+1):length(bhat.table$L)],
           bhat.table$delta.log.N1.plus[(last_choices[2]+1):length(bhat.table$L)],
           labels=row.names(bhat.table.list[[1]][
             (last_choices[2]+1):length(bhat.table$delta.log.N1.plus),]),
           cex= 0.7, pos=3)
    }
    # Second plot with the histogramms
#     hist(bhat.table.list[[1]]$N1.plus)
#     plot(catch.vec ~ bhat.table.list[[1]]$mean.length.classes, type  = 'h')


    id.co1 <- identify(bhat.table$L,bhat.table$delta.log.N1.plus,
                       n = 2, pos = TRUE)
    if(length(id.co1$ind) == 0){
      break
    }
    colour.vec[id.co1$ind[1]:id.co1$ind[2]] <- colour.xy[xy]
    dev.off()

    #STEP 7: calculate mean length and standard deviation of regression line
    x.co1 <- bhat.table$L[id.co1$ind[1]:id.co1$ind[2]]
    y.co1 <- bhat.table$delta.log.N1.plus[id.co1$ind[1]:id.co1$ind[2]]
    m.co1 <- lm(y.co1 ~ x.co1)
    sum.m.co1 <- summary(m.co1)
    a.co1 <- sum.m.co1$coefficients[1]
    b.co1 <- sum.m.co1$coefficients[2]
    l.mean.co1 <- -a.co1/b.co1   #mean length: L(mean)(N1) = -a/b
    s.co1 <- sqrt(-1/b.co1)      #standard deviation: s(N1) = sqrt(-1/b)

    a.b.list[[xy]] <- c(a.co1,b.co1)

    #STEP 8: fill sixth column
    #set.seed ????
    normal.dis.co1 <- rnorm(n=1000, mean = l.mean.co1, sd=s.co1)
    max.class.ind.co1 <- which(round(max(normal.dis.co1), digits = 0) >= (
      bhat.table$mean.length.classes - (interval/2)) &
        round(max(normal.dis.co1),digits = 0) < (
          bhat.table$mean.length.classes + (
            interval/2)))
    max.class.co1 <- bhat.table$mean.length.classes[max.class.ind.co1]
    for(i in 1:max.class.ind.co1){
      delta.N <- round(a.co1 + b.co1 * bhat.table$L[i],digits=3)
      bhat.table$delta.log.N[i] <- delta.N
    }


    #STEP 9: fills one value in seventh column and one in eigth column, get clean starting value
#     one.left.max <- which(bhat.table$N1.plus ==
#                             max(bhat.table$N1.plus[which(!is.na(
#                               bhat.table$delta.log.N))],na.rm=T)) - 1   #OWN ASSUMPTION: one size class left of the maximum should be not influenced by the normal distribution right to this distribution
#     one.left.max <- one.left.max[1]
#     if(is.na(bhat.table$log.N1.plus[one.left.max])){
#       one.left.max <- one.left.max + 1
#     }
#     bhat.table$log.N1[one.left.max] <- log(bhat.table$N1.plus[
#       one.left.max])
#     bhat.table$N1[one.left.max] <- bhat.table$N1.plus[one.left.max]

    # NEW ASSUMPTION: (FIND OUT!)
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
    #for(i in (one.left.max+1):length(bhat.table$delta.log.N)){
    for(i in (mid_choice+1):length(bhat.table$delta.log.N)){
      log.N1 <- bhat.table$log.N1[i-1] + bhat.table$delta.log.N[i]
      bhat.table$log.N1[i] <- log.N1
    }

    #STEP 11: fills rest of eigth column
    #for(i in (one.left.max+1):length(bhat.table$delta.log.N)){
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
    #bhat.table.list[[xy+1]] <- bhat.table2 #or bhat.table1?? which one more important
    bhat.table.list[[xy+1]] <- bhat.table2 #or bhat.table1?? which one more important

    #necessary for stop function to work
    max.rem <- id.co1$ind[2]
    last_choices <- id.co1$ind
    rm(id.co1)

  }

  # Combine results of regression analyses for output
  bhat.table.list <- bhat.table.list[!sapply(bhat.table.list, is.null)]
  bhat.table.list2 <- bhat.table.list[-1]
  bhat.table.list_new <- lapply(bhat.table.list2, function(x)  return(x[,3:9]))
  bhat.results <- do.call(cbind, bhat.table.list_new)
  if(!is.null(bhat.results)){
    ret_pri <- cbind(bhat.table.list[[1]], bhat.results)
  }else ret_pri <- bhat.table.list[[1]]


  ret <- list(res, bhat_results = ret_pri)
  return(ret)
}



#put that into a loop for the number of cohorts one wants to analyse,
#in the plot give colour to all points which already were analysed or used for regression,
#save tables to a list
#save other important information to a list, such as total number in cohort (sum of N1 or N2), and  a and b ans stuffs
#question if satisfied with choosing
# have a splitted graphic device with plot to choose points and histogram where choosen
# Gaussian probability function is displayed. Ask user if "redo" choosing.
# First thing to show shoud be a histogram, then question if fixed number of cohort, if
# yes how many, or no, then user can choose freely



