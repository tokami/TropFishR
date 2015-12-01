#' @title Bhattacharya's method
#'
#' @description Separating normal distributions of several cohorts in distinct
#'   distributions representing cohorts. Second sentence.
#'
#' @param param A list consisting of following parameters:
#'   \code{$age} or \code{$midLengths} midpoints of the length class as vector (length frequency
#'   data) or ages as vector (age composition data),
#'   \code{catch} Catch as vector, or a matrix with catches of subsequent years if
#'   the catch curve with constat time intervals should be applied;
#'
#' @examples
#' \donttest{
#'  data(synLFQ1)
#'  Bhattacharya(synLFQ1)
#' }
#' @details Bhattacharya
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export
#'

Bhattacharya <- function(param){

  res <- param
  midLengths <- as.character(res$midLengths)
  if("catch_mat" %in% names(res) == TRUE) catch <- res$catch_mat
  if("catch" %in% names(res) == TRUE) catch <- res$catch

  #Transform all input variables to correct class
  midLengths <- as.numeric(as.character(midLengths))
  #catch <- as.numeric(as.character(catch))

  #Transform matrix into vector if provided
  if(class(catch) == 'matrix'){
    catch.vec <- rowSums(catch)
  }else catch.vec = catch
  if(length(midLengths) != length(catch.vec)) stop("midLengths and catch do not have the same length!")
  df.Bh <- data.frame(midLengths.Bh = midLengths,
                      catch.Bh = catch.vec)

  #calculate size class interval
  interval.Bh <- df.Bh$midLengths.Bh[2] - df.Bh$midLengths.Bh[1]

  #STEP 1: fills second column of bhat with counts per size class
  bhat.table <- data.frame('mean.length.classes' = df.Bh$midLengths.Bh,
                           'N1.plus'= df.Bh$catch.Bh,
                           'log.N1.plus'=NA,
                           'delta.log.N1.plus'=NA,
                           'L'=NA,
                           'delta.log.N'=NA,
                           'log.N1' =NA,
                           'N1' =NA,
                           'N2.plus'=NA)

  bhat.table.list <- list()
  a.b.list <- list()
  colour.xy <- c('blue','darkgreen','red','yellow','purple','orange','lightgreen','skyblue')
  colour.vec <- rep('black',length(bhat.table$L))


  for(xy in 1:30){  #more than 30 cohorts?
    #STEP 2: fills third column
    bhat.table$log.N1.plus <- round(log(bhat.table$N1.plus),digits=3)

    #STEP 3: fills fourth coulmn
    for(i in 2:length(bhat.table$log.N1.plus)){

      delta.value <- bhat.table$log.N1.plus[i] - bhat.table$log.N1.plus[i-1]
      bhat.table$delta.log.N1.plus[i] <- delta.value
    }

    #STEP 4: fills fifth coulmn
    for(i in 1:length(bhat.table$log.N1.plus)){
      if(!is.na(bhat.table$delta.log.N1.plus[i]) & bhat.table$delta.log.N1.plus[i] != 'Inf'){
        min.size <- (bhat.table$mean.length.classes[i] - (interval.Bh/2))
        bhat.table$L[i] <- min.size
      }
    }

    #STEP 5: plot of fifth against fourth column
    plot(bhat.table$delta.log.N1.plus ~ bhat.table$L, pch = 16)
    abline(h=0)
    p.bhat1 <- recordPlot()

    #STEP 6: select points for regression line
    if(xy == 1){bhat.table.list[[1]] <- bhat.table}

    # do not include in linear regression if number in N1.plus is under 3 or 4 if(bhat.table$N1.plus < 3){}

    #find new maximum for start of next cohort
    if(xy > 1){
      gradi <- NA
      for(nm in 3:(length(bhat.table.list[[1]]$delta.log.N1.plus)-1)){
        if((bhat.table.list[[1]]$delta.log.N1.plus[nm] >
            bhat.table.list[[1]]$delta.log.N1.plus[nm-1]) == T &
           (bhat.table.list[[1]]$delta.log.N1.plus[nm] >
            bhat.table.list[[1]]$delta.log.N1.plus[nm+1]) == T
        ){gradi[nm] <- T}else gradi[nm] <- F
      }

      cut.max <- max.rem + min(which(gradi[(max.rem+1):length(gradi)] == T))

      bhat.table.list[[1]]$delta.log.N1.plus[
        which(as.numeric(rownames(bhat.table.list[[1]])) >= cut.max)] <-
        bhat.table$delta.log.N1.plus[
          which(as.numeric(rownames(bhat.table.list[[1]])) >= cut.max)]
      }


    plot(bhat.table.list[[1]]$L,bhat.table.list[[1]]$delta.log.N1.plus, pch = 16,
         col = colour.vec)
    abline(h = 0)
    text(bhat.table.list[[1]]$L,bhat.table.list[[1]]$delta.log.N1.plus,
         labels=row.names(bhat.table.list[[1]]), cex= 0.7, pos=3)
    if(xy > 1){
      abline(a = a.b.list[[xy-1]][1],b=a.b.list[[xy-1]][2],col=colour.xy[xy-1])
    }
    print("Starting on the left, please choose from black points which lie on a straight line! Do not include points which might be affected by the next distribution!")
    id.co1 <- identify(bhat.table$L,bhat.table$delta.log.N1.plus,
                       n = 2, pos = TRUE)
    if(length(id.co1$ind) == 0){
      break
    }
    colour.vec[id.co1$ind[1]:id.co1$ind[2]] <- colour.xy[xy]


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
    normal.dis.co1 <- rnorm(n=1000,mean = l.mean.co1,sd=s.co1)
    max.class.ind.co1 <- which(round(max(normal.dis.co1),digits = 0) >= (
      bhat.table$mean.length.classes - (interval.Bh/2)) &
        round(max(normal.dis.co1),digits = 0) < (
          bhat.table$mean.length.classes + (
            interval.Bh/2)))
    max.class.co1 <- bhat.table$mean.length.classes[max.class.ind.co1]
    for(i in 1:max.class.ind.co1){
      delta.N <- round(a.co1 + b.co1 * bhat.table$L[i],digits=3)
      bhat.table$delta.log.N[i] <- delta.N
    }


    #STEP 9: fills one value in seventh column and one in eigth column, get clean starting value
    one.left.max <- which(bhat.table$N1.plus ==
                            max(bhat.table$N1.plus[which(!is.na(
                              bhat.table$delta.log.N))],na.rm=T)) - 1   #OWN ASSUMPTION: one size class left of the maximum should be not influenced by the normal distribution right to this distribution
    one.left.max <- one.left.max[1]
    if(is.na(bhat.table$log.N1.plus[one.left.max])){
      one.left.max <- one.left.max + 1
    }
    bhat.table$log.N1[one.left.max] <- log(bhat.table$N1.plus[
      one.left.max])
    bhat.table$N1[one.left.max] <- bhat.table$N1.plus[one.left.max]

    #STEP 10: fills rest of seventh column
    for(i in (one.left.max+1):length(bhat.table$delta.log.N)){
      log.N1 <- bhat.table$log.N1[i-1] + bhat.table$delta.log.N[i]
      bhat.table$log.N1[i] <- log.N1
    }

    #STEP 11: fills rest of eigth column
    for(i in (one.left.max+1):length(bhat.table$delta.log.N)){
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
                             'N1.plus'= bhat.table2$N2.plus,
                             'log.N1.plus'=NA,
                             'delta.log.N1.plus'=NA,
                             'L'=NA,
                             'delta.log.N'=NA,
                             'log.N1' =NA,
                             'N1' =NA,
                             'N2.plus'=NA)

    #arrange dataframe for saving
    colnames(bhat.table2) <- c('mean.length.classes',paste('N',xy,'.plus',sep=''),
                               paste('log.N',xy,'.plus',sep=''),
                               paste('delta.log.N',xy,'.plus',sep=''),'L',
                               'delta.log.N',paste('log.N',xy,sep=''),
                               paste('N',xy,sep=''),paste('N',xy+1,'.plus',sep=''))

    #save data
    bhat.table.list[[xy+1]] <- bhat.table2 #or bhat.table1?? which one more important

    #necessary for stop function to work
    max.rem <- id.co1$ind[2]
    rm(id.co1)
    print(bhat.table2)
  }
}



#put that into a loop for the number of cohorts one wants to analyse,
#in the plot give colour to all points which already were analysed or used for regression,
#save tables to a list
#save other important information to a list, such as total number in cohort (sum of N1 or N2), and  a and b ans stuffs
#question if satisfied with choosing



