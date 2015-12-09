#' @title ELEFAN
#'
#' @description Electronic length frequency analysis for estimating growth parameter of exploited fish populations.
#'
#' @param param A list consisting of following parameters:
#'   \code{$age} or \code{$midLengths} midpoints of the length class as vector (length frequency
#'   data) or ages as vector (age composition data),
#'   \code{catch} Catch as vector, or a matrix with catches of subsequent years if
#'   the catch curve with constat time intervals should be applied;
#' @param range.Linf lower and upper limit of range of Linf values
#' @param step.Linf step size between Linf values
#' @param range.K lower and upper limit of range of K values
#' @param step.K step size between K values
#' @param t0 theoretical age of fish at length zero
#' @param tmax maximal age of fish from literature only estimate needed
#'
#' @examples
#' \donttest{
#' data(trout)
#' param = trout
#' range.Linf <- c(11,15)
#' step.Linf <- 1
#' range.K <- c(0.2,1)
#' step.K <- 0.1
#' t0 <- 0
#' interval = 1
#' tmax = 5
#' }
#' @details ELEFAN
#'
#' @references
#' Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the objective extraction of
#' growth parameters from length-frequency data. \emph{Meeresforschung}, 28(4):205-211
#'
#' Pauly, D., 1987. A review of the ELEFAN system for analysis of length-frequency data in
#' fish and aquatic invertebrates. \emph{ICLARM Conf. Proc.}, (13):7-34
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export
#'

ELEFAN <- function(param, range.Linf, step.Linf,
                   range.K, step.K, t0 = 0,
                   tmax){

  res <- param
  classes <- res$midLengths
  catch <- res$catch
  interval <- classes[1]-classes[2]


  #________________________DATA ARRANGEMENT________________________#
  #--------------

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #         moving average           #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  catch.MA <- matrix(ncol=dim(catch)[2],nrow=dim(catch)[1],NA)
  vec.MA <- rep(NA,length(classes))
  for(i in 1:dim(catch)[2]){
    actual.month <- catch[,i]

    for(j in 1:length(actual.month)){

      if(j == 1){
        MA <- sum(actual.month[j:(j+2)]) / 5
      }else if(j == 2){
        MA <- sum(actual.month[(j-1):(j+2)]) / 5
      }else if(j > 2 & j < (length(actual.month)-1)){
        MA <- sum(actual.month[(j-2):(j+2)]) / 5
      }else if(j == (length(actual.month)-1)){
        MA <- sum(actual.month[(j-2):(j+1)]) / 5
      }else if(j == length(actual.month)){
        MA <- sum(actual.month[(j-2):j]) / 5
      }

      vec.MA[j] <- MA
    }

    #vec.MA[1:2] <- 0
    #vec.MA[(length(actual.month)-1):length(actual.month)] <- 0

    catch.MA[,i] <- vec.MA
  }
  #make sure: that in beginning and in the end is exactly two 0!!!

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #        adjusted frequency        #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  catch.AF <- matrix(ncol=dim(catch)[2],nrow=dim(catch)[1],NA)
  vec.AF <- rep(NA,length(classes))
  for(i in 1:dim(catch)[2]){

    actual.month <- catch[,i]
    actual.month.MA <- catch.MA[,i]

    for(j in 1:length(actual.month.MA)){

      if(actual.month[j] != 0 & actual.month.MA[j] != 0){
        AF <- (actual.month[j]/ actual.month.MA[j])
      }else AF <- 0

      vec.AF[j] <- AF
    }

    catch.AF[,i] <- vec.AF
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #    relative adjusted frequency   #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  catch.rAF <- matrix(ncol=dim(catch)[2],nrow=dim(catch)[1],NA)
  vec.rAF <- rep(NA,length(classes))
  for(i in 1:dim(catch)[2]){
    actual.month.AF <- catch.AF[,i]
    mean.AF <- mean(actual.month.AF,na.rm = T)

    for(j in 1:length(actual.month.AF)){

      if(mean.AF != 0){
        rAF <- (actual.month.AF[j]/ mean.AF) -1
      }else stop('The mean of adjusted frequencies must be unequal to zero')

      vec.rAF[j] <- rAF
    }

    catch.rAF[,i] <- vec.rAF
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #           Adjustements           #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #account for inflation of adjacent zero frequencies = nz
  #low values
  # -1 to 0
  catch.aAF <- matrix(ncol=dim(catch)[2],nrow=dim(catch)[1],NA)
  vec.aAF <- rep(NA,length(classes))
  vec.aAF.adj <- rep(NA,length(classes))

  for(i in 1:dim(catch)[2]){
    actual.month <- catch[,i]
    actual.month.rAF <- catch.rAF[,i]

    for(j in 1:length(actual.month.rAF)){
      #account for inflation of adjacent zero frequencies = nz
      if(actual.month.rAF[j] != -1 & j == 1){
        vec <- c(actual.month[j+1],actual.month[j+2])
        nu.zero <- length(vec[which(vec == 0)])
        if(nu.zero == 0){
          aAF <- actual.month.rAF[j] * 0.25
        }else if(nu.zero == 1){
          aAF <- actual.month.rAF[j] * 0.125
        }else if(nu.zero == 2){
          aAF <- actual.month.rAF[j] * 0.0625
        }
      }else if(actual.month.rAF[j] != -1 & j == 2){
        vec <- c(actual.month[j-1],actual.month[j+1],actual.month[j+2])
        nu.zero <- length(vec[which(vec == 0)])
        if(nu.zero == 0){
          aAF <- actual.month.rAF[j] * 0.5
        }else if(nu.zero == 1){
          aAF <- actual.month.rAF[j] * 0.25
        }else if(nu.zero == 2){
          aAF <- actual.month.rAF[j] * 0.125
        }else if(nu.zero == 3){
          aAF <- actual.month.rAF[j] * 0.0625
        }
      }else if(actual.month.rAF[j] != -1 & j > 2 & j < (length(actual.month.rAF)-1)){
        vec <- c(actual.month[j-2],actual.month[j-1],actual.month[j+1],actual.month[j+2])
        nu.zero <- length(vec[which(vec == 0)])
        if(nu.zero == 0){
          aAF <- actual.month.rAF[j]
        }else if(nu.zero == 1){
          aAF <- actual.month.rAF[j] * 0.5
        }else if(nu.zero == 2){
          aAF <- actual.month.rAF[j] * 0.25
        }else if(nu.zero == 3){
          aAF <- actual.month.rAF[j] * 0.125
        }else if(nu.zero == 4){
          aAF <- actual.month.rAF[j] * 0.0625
        }
      }else if(actual.month.rAF[j] != -1 & j == (length(actual.month.rAF)-1)){
        vec <- c(actual.month[j-2],actual.month[j-1],actual.month[j+1])
        nu.zero <- length(vec[which(vec == 0)])
        if(nu.zero == 0){
          aAF <- actual.month.rAF[j] * 0.5
        }else if(nu.zero == 1){
          aAF <- actual.month.rAF[j] * 0.25
        }else if(nu.zero == 2){
          aAF <- actual.month.rAF[j] * 0.125
        }else if(nu.zero == 3){
          aAF <- actual.month.rAF[j] * 0.0625
        }
      }else if(actual.month.rAF[j] != -1 & j == length(actual.month.rAF)){
        vec <- c(actual.month[j-2],actual.month[j-1])
        nu.zero <- length(vec[which(vec == 0)])
        if(nu.zero == 0){
          aAF <- actual.month.rAF[j] * 0.25
        }else if(nu.zero == 1){
          aAF <- actual.month.rAF[j] * 0.125
        }else if(nu.zero == 2){
          aAF <- actual.month.rAF[j] * 0.0625
        }
      }else if(actual.month.rAF[j] == -1){
        aAF <- -1
      }


      #low values
      if(actual.month[j] <= 10 & actual.month[j] > 0 & actual.month.rAF[j] > 0){
        aAF <- aAF / (sqrt(1+(2/actual.month[j])))
      }

      # -1 to 0
      if(aAF == -1){
        aAF <- 0
      }

      vec.aAF[j] <- aAF
    }

    #sum(+)/sum(-)
    vec.aAF.pos <- vec.aAF[which(vec.aAF > 0)]
    sum.pos <- sum(vec.aAF.pos)
    vec.aAF.neg <- vec.aAF[which(vec.aAF < 0)]
    sum.neg <- sum(vec.aAF.neg)

    ratio.SUM <- sum.pos/abs(sum.neg)

    vec.aAF[which(vec.aAF <= 0)] <- vec.aAF[which(vec.aAF <= 0)] * ratio.SUM

    vec.aAF.adj <- vec.aAF

    #ultimate negative to zero
    if((vec.aAF.adj[length(vec.aAF.adj)] < 0) == T){
      vec.aAF.adj[length(vec.aAF.adj)] <- 0
    }

    #penultimate negative dividing by 2
    if((vec.aAF.adj[(length(vec.aAF.adj)-1)] < 0) == T){
      vec.aAF.adj[(length(vec.aAF.adj)-1)] <- vec.aAF.adj[(length(vec.aAF.adj)-1)] / 2
    }


    catch.aAF[,i] <- vec.aAF.adj
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #        Calculation of ASP        #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  pos.NA.df <- catch.aAF
  for(i in 1:dim(catch)[2]){
    pos.NA <- catch.aAF[,i]
    pos.NA[which(catch.aAF[,i] <= 0)]  <- NA

    j = 0
    next.one = 1
    repeat{
      j = j + next.one

      if(!is.na(pos.NA[j]) & j < length(pos.NA)){
        k = 0
        pos <- NA
        pos[1] <- pos.NA[j]
        repeat{
          k = k + 1

          if(!is.na(pos.NA[j+k])){
            pos[k+1] <- pos.NA[j+k]
          }else if(is.na(pos.NA[j+k]) | (j+k) == length(pos.NA)) break
        }
        if(length(pos) > 1){
          highest <- (which(pos == max(pos)) - 1) + j
          to.delete <- seq(j,(j+(k-1)),1)
          to.delete <- to.delete[which(to.delete != highest)]
          pos.NA[to.delete] <- NA
          next.one <- k
        }else next.one <- 1
      }else next.one <- 1
      if(j > length(pos.NA)) break
    }
    pos.NA.df[,i] <- pos.NA
  }

  vec.ASP <- rep(NA,length(classes))
  for(i in 1:dim(pos.NA.df)[2]){
    vec.ASP[i-2] <- sum(pos.NA.df[,i],na.rm=T)
  }

  #assuming that there is one ASP value for a whole data set, going on:
  ASP <- sum(vec.ASP,na.rm=T)


  # create peak matrix
  prep_mat <- catch.aAF
  prep_mat <- ifelse(prep_mat > 0,1,prep_mat)
  prep_mat <- ifelse(prep_mat <= 0,0,prep_mat)
  peak.list <- list()
  for(coli in 1:dim(prep_mat)[2]){
    vec.peaki <- prep_mat[,coli]
    rle.val <- rle(vec.peaki)$values
    rle.val[which(rle.val == 1)] <- 1:length(rle.val[which(rle.val == 1)])
    peak.list[[coli]] <- rep(rle.val,rle(vec.peaki)$lengths)
  }
  peaks_mat <- do.call(cbind,peak.list)


  #     calculates the maximum sum of points "available" in a (set of) length-frequency sample(s) (see Fig. 1 C) ["available points" refers here to points which can possibly be "accumulated" by one single growth curve; see below]. This sum is termed "available sum of peaks" (ASP).
  #     "traces" through the (set of) length-frequency sample(s) sequentially arranged in time, for any arbitrary "seed" input of Lex> and K, a series of growth curves started from the base of each of the peaks, and projected backward and forward in time to meet all other samples of the sample set (Fig. 2) and/or the same sample repeated again and again (Fig. 3).
  #     accumulates the "points" obtained by each growth curve when passing through peaks (positive points) or through the troughs separating peaks (negative points) (see Fig. 1 B and C).
  #     selects the curve which, by passing through most peaks and avoiding most troughs best "explains" the peaks in the (set of) sample(s) and therefore accumulates the largest number of points. This new sum is called "explained sum of peaks" (ESP). decrements or increments the "seeded" values of Lex> and K until the ratio ESP/ASP reaches a maximum, and gives the growth parameters corresponding to this optimum ratio.

  #-------

  #________________________GROWTH CURVE________________________#
  #problem of arrangement of t, in literature: samples sequentially arranged in time

  #for t
  #retrieve dates from headers of columns
  names(catch)[1:length(names(catch))] <- unlist(strsplit(names(catch)[1:length(names(catch))],split = 'X'))[seq(2,(length(names(catch))*2),2)]
  #bring sample dates in right format
  dates <- as.POSIXlt(as.Date(names(catch)[1:length(names(catch))],format = "%d.%m.%Y"))
  dates.all <- as.Date(dates)
  #get continuous time in years
  #alternative: + - 100 days as full range
  days <- as.numeric(dates.all - as.Date((dates.all[1] - 100)))     #cut(dates.all[1], "year")))
  days.years <- days/365
  #sampling period
  sample.period.days <- days[length(days)] - days[1]
  sp.years <- sample.period.days/365

  t.days <- seq(1,days[length(days)],1)
  t.years <- (t.days)/365
  #  full range + / - days to sampling times
  t.years.fr <- seq(1,days[length(days)]+100,1)/365

  rownames(catch.aAF) <- classes
  colnames(catch.aAF) <- days.years   ## OLD: days

  Linfs <- seq(range.Linf[1],range.Linf[2],step.Linf)
  Ks <- seq(range.K[1],range.K[2],step.K)

  #t <- t.days
  t <- t.years   ## for whoel range of t max?
  #for whole live of fish:
  t <- 0:(tmax * 365) / 365


  # get other growth curves going through data which is dependent on tmax
  if(sp.years > 0 & sp.years < 1) repro.add <- 0
  if(sp.years > 1 & sp.years < 2) repro.add <- 1
  if(sp.years > 2 & sp.years < 3) repro.add <- 2
  if(sp.years > 3 & sp.years < 4) repro.add <- 3
  if(sp.years > 4 & sp.years < 5) repro.add <- 4

  samp.times <- as.numeric(colnames(catch.aAF))   # what about the next cohort one year later or one year ago if totap period is over one year???



sys_timeBF <- Sys.time()


  ESP.tshift.L <- list()
  ESP.list.L <- list()
  for(li in 1:length(Linfs)){

    ESP.tshift.k <- list()
    ESP.list.k <- list()
    for(ki in 1:length(Ks)){

      # VBGF
      Lt <-  Linfs[li] * ( 1 - exp( - Ks[ki] * (t - t0)))
      Lt.cohis <- tmax + repro.add

      # get lengths which fall into sampling times
      # improve these loops by defining when smaller than smallest length group or larger than largest length group it turns into NA and not the extreme groups
      options(warn = -1) # turn warning messages off globally
      lt.classes.list <- list()
      for(tshift in 1:(length(t.days)-1)){ # produces warnings because latest cohort will not be represented by data and thereby produce NA
        ttx <- t.years[tshift]

        lt.class.cohi <- list()
        for(cohi in -repro.add:tmax){
          diff <- samp.times - (ttx - cohi)
          diff[which(diff < 0)] <- NA
          lt.classes <- Lt[(diff*365)]

          # get corresponding length groups instead of exact lengths
          for(hitx in 1:length(lt.classes)){
            if(all(is.na(classes - lt.classes[hitx]))){lt.classes[hitx] <- NA
            }else lt.classes[hitx] <- which.min(abs(classes - lt.classes[hitx]))
          }

          lt.class.cohi[[cohi+repro.add+1]] <- lt.classes
        }
        lt.class.cohi.df <- do.call(cbind,lt.class.cohi)
        lt.classes.list[[tshift]] <- data.frame(tshift = rep(tshift,
                                                             length(lt.class.cohi.df[,1])),
                                                lt.class.cohi.df)
      }
      options(warn = 0)


      # get ESP values from catch.aAF matrix entries which are hit by growth curve
      ESP.list <- list()
      for(inds in 1:length(lt.classes.list)){

        df <- lt.classes.list[[inds]]
        loop_mat <- catch.aAF
        ESPs <- NA
        ESP.loop <- list()
        for(indi in 2:length(df)){
          res.ind <- df[,indi]
          for(indx in 1:length(res.ind)){

            ESPs[indx] <- loop_mat[res.ind[indx],indx]

            #flagging out of hit peaks (Pauly,1985)
            if(!is.na(res.ind[indx]) & peaks_mat[res.ind[indx],indx] != 0){
              peakX <- peaks_mat[res.ind[indx],indx]
              loop_mat[which(peaks_mat[,indx] != 0 &
                               peaks_mat[,indx] == peakX ),indx] <- NA
            }

          }
        ESP.loop[[indi-1]] <- ESPs
        }
        # get ESP by summing up
        ESP.list[[inds]] <- c(inds,sum(unlist(ESP.loop),na.rm=TRUE))
      }


      # choose best ESP value of all starting times (tshift)
      ESP.pre.df <- do.call(rbind,ESP.list)
      ESPx <- ESP.pre.df[which(ESP.pre.df[,2] == max(ESP.pre.df[,2],na.rm=TRUE)),]   ### here more than one combination of thsift and ESP can be resulting

      if(is.matrix(ESPx)) ESPx <- ESPx[1,]           ### with the one before the comma it isdefined that always the first of multiple best fits is taking (= the one with the smallest tshift)
      ESP.tshift.k[[ki]] <- c(Linfs[li],Ks[ki],ESPx[1])
      ESP.list.k[[ki]] <- c(Linfs[li],Ks[ki],ESPx[2])
    }
    ESP.tshift.L[[li]] <- do.call(rbind,ESP.tshift.k)
    ESP.list.L[[li]] <- do.call(rbind,ESP.list.k)
  }

  ESP.df <- do.call(rbind,ESP.list.L)
  tapply(ESP.df[,3],list(ESP.df[,1],ESP.df[,2]),round,digits = 2)
  ESP.time <- do.call(rbind,ESP.tshift.L)
  tapply(ESP.time[,3],list(ESP.time[,1],ESP.time[,2]),round,digits = 2)
  sys_timeAW <- Sys.time()
  sys_timeAW - sys_timeBF
  #PEAKS which are hit have to be flagged out, troughs can be hit several times

10^(8.78/ASP) /10

  #   tvar.list <- list()
  #   for(tvar in 0:length(t.days)){
  #     Lt <-  Linfs[li] * ( 1 - exp( - K[ki] * ((t+tvar) - t0)))
  #     tvar.list[[tvar+1]] <- Lt
  #   }


  #for plotting
  xplot <- c(min(t.years.fr),(max(t.years.fr)))
  yplot <- c(classes[1]-interval,classes[length(classes)]+interval)

  # Plot with rearranged histogramms
  par(mar = c(5, 5, 1, 1) + .1)
  plot(xplot, yplot, type = "n",axes = FALSE,
       ann = FALSE,  xaxs = "i", yaxs = "i")
  text(days.years, par("usr")[3] - 0.7,
       labels = dates.all, srt = 45, pos = 1, xpd = TRUE)
  axis(2, cex.axis = 1.2)
  mtext(side = 2, outer = F, line = 3.5, "Length [cm]", cex = 1.2)
  box( col = "gray40") #bty = "L"

  #Histograms
  # positive histograms or peaks
  y.b <- classes - interval / 2
  y.t <- classes + interval / 2
  for(histi in 1:dim(catch)[2]){
    x.l <- rep(days.years[histi],length(classes))
    x.x <- catch.aAF[,histi]
    x.x[which(x.x < 0)] <- 0
    x.r <- x.l - x.x * 0.05   ### make this number dependent on data!?!?
    rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t,
         col = "gray40", border = "black")
  }

  # negative histograms or troughs
  for(histi in 1:dim(catch)[2]){
    x.l <- rep(days.years[histi],length(classes))
    x.x <- catch.aAF[,histi]
    x.x[which(x.x > 0)] <- 0
    x.r <- x.l - x.x * 0.05   ### make this number dependent on data!?!?
    rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t,
         col = "gray80", border = "black")
  }



#   for(xyx in 1:length(Lt.list)){
#     lines(y = Lt.list[[xyx]], x = t, lty=2, col='darkgreen')
#   }


  Linf.p = 44.3
  K.p = 0.41
  tshift.p = 422/365
  # growth curves
  Lt <- Linf.p * (1 - exp(-K.p * (t.years.fr - tshift.p)))    ### choose correct t -> tshift
  lines(y = Lt, x = t.years.fr, lty=2, col='blue')
  #t.years.frA <- seq(xplot[1],xplot[2],1/365)
  Lt <- Linf.p * (1 - exp(-K.p * (t.years.fr - tshift.p - 1)))
  lines(y = Lt, x = t.years.fr, lty=2, col='blue')
  Lt <- Linf.p * (1 - exp(-K.p * (t.years.fr - tshift.p + 1)))
  lines(y = Lt, x = t.years.fr , lty=2, col='blue')
  Lt <- Linf.p * (1 - exp(-K.p * (t.years.fr - tshift.p + 2)))
  lines(y = Lt, x = t.years.fr, lty=2, col='blue')
  Lt <- Linf.p * (1 - exp(-K.p * (t.years.fr - tshift.p + 3)))
  lines(y = Lt, x = t.years.fr, lty=2, col='blue')




  #   #start with seed values for L8 and K and create growth functions from every peak
  #   peaks <- which(!is.na(pos.NA.df))
  #   #project backward and forward in time? wtf
}




#for plotting
xplot <- c(min(t.years.fr),(max(t.years.fr)))     # display months on x axis
yplot <- c(0,classes[length(classes)]+interval)

# Plot with rearranged histogramms
par(mar = c(5, 5, 1, 1) + .1)
plot(xplot, yplot, type = "n",axes = FALSE,
     ann = FALSE,  xaxs = "i", yaxs = "i")
text(days.years, par("usr")[3] - 0.7,
     labels = dates.all, srt = 45, pos = 1, xpd = TRUE)
axis(2, cex.axis = 1.2)
mtext(side = 2, outer = F, line = 3.5, "Length [cm]", cex = 1.2)
box( col = "gray40") #bty = "L"

#Histograms
# true distributions
for(histi in 1:dim(catch)[2]){
  x.r <- rep(days.years[histi],length(classes))
  x.x <- catch[,histi]
  x.l <- x.r - x.x * 0.002   ### make this number dependent on data!?!?
  rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t,
       col = "gray80", border = "black")
}


#   for(xyx in 1:length(Lt.list)){
#     lines(y = Lt.list[[xyx]], x = t, lty=2, col='darkgreen')
#   }







## OLD
# get size classes in which calculated Lt values fall into
for(hits in 1:length(lt.classes.list)){
  df <- lt.classes.list[[hits]]
  for(hiti in 2:length(df)){
    hi.res <- df[,hiti]
    for(hitx in 1:length(hi.res)){
      if(all(is.na(classes - hi.res[hitx]))){hi.res[hitx] <- NA
      }else hi.res[hitx] <- which.min(abs(classes - hi.res[hitx]))
    }
    df[,hiti] <- hi.res
  }
  lt.classes.list[[hits]] <- df
}



# get size classes in which calculated Lt values fall into
for(hits in 1:length(lt.classes.list)){
  df <- lt.classes.list[[hits]]
  for(hiti in 2:length(df)){
    for(hitx in 1:length(df[,2])){
      if(all(is.na(classes - df[hitx,hiti]))){df[hitx,hiti] <- NA
      }else df[hitx,hiti] <- which.min(abs(classes - df[hitx,hiti]))
    }
  }
  lt.classes.list[[hits]] <- df
}
