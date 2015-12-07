#' @title ELEFAN
#'
#' @description Electronic length frequency analysis for estimating growth parameter of exploited fish populations.
#'
#' @param param A list consisting of following parameters:
#'   \code{$age} or \code{$midLengths} midpoints of the length class as vector (length frequency
#'   data) or ages as vector (age composition data),
#'   \code{catch} Catch as vector, or a matrix with catches of subsequent years if
#'   the catch curve with constat time intervals should be applied;
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
#' }
#' @details ELEFAN
#'
#' @references
#' Pauly and David 1981
#'
#' Pauly 1987
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export
#'


ELEFAN <- function(param, range.Linf, step.Linf, range.K, step.K, t0 = NA, interval){

  res <- param
  classes <- res$midLengths
  catch <- res$catch
  interval <- interval


  #________________________DATA ARRANGEMENT________________________#

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


#     calculates the maximum sum of points "available" in a (set of) length-frequency sample(s) (see Fig. 1 C) ["available points" refers here to points which can possibly be "accumulated" by one single growth curve; see below]. This sum is termed "available sum of peaks" (ASP).
#     "traces" through the (set of) length-frequency sample(s) sequentially arranged in time, for any arbitrary "seed" input of Lex> and K, a series of growth curves started from the base of each of the peaks, and projected backward and forward in time to meet all other samples of the sample set (Fig. 2) and/or the same sample repeated again and again (Fig. 3).
#     accumulates the "points" obtained by each growth curve when passing through peaks (positive points) or through the troughs separating peaks (negative points) (see Fig. 1 B and C).
#     selects the curve which, by passing through most peaks and avoiding most troughs best "explains" the peaks in the (set of) sample(s) and therefore accumulates the largest number of points. This new sum is called "explained sum of peaks" (ESP). decrements or increments the "seeded" values of Lex> and K until the ratio ESP/ASP reaches a maximum, and gives the growth parameters corresponding to this optimum ratio.


  #________________________GROWTH CURVE________________________#
  #problem of arrangement of t, in literature: samples sequentially arranged in time

  #for t
  #retrieve dates from headers of columns
  names(catch)[1:length(names(catch))] <- unlist(strsplit(names(catch)[1:length(names(catch))],split = 'X'))[seq(2,(length(names(catch))*2),2)]
  #bring sample dates in right format
  dates <- as.POSIXlt(as.Date(names(catch)[1:length(names(catch))],format = "%d.%m.%Y"))
  dates.all <- as.Date(dates)
  #get continuous time in years
  days <- as.numeric(dates.all - as.Date(cut(dates.all[1], "year")))  # takes beginning of year, what if one sampled on 1st of January? better to take a certain time period before first sampling time
  #alternative: + - 100 days as full range
  days <- as.numeric(dates.all - as.Date((dates.all[1] - 100)))     #cut(dates.all[1], "year")))
  days.years <- days/365
  t.days <- seq(1,days[length(days)],1)
  t.years <- (t.days)/365
  #  full range + / - days to sampling times
  t.years.fr <- seq(1,days[length(days)]+100,1)/365

  rownames(catch.aAF) <- classes
  colnames(catch.aAF) <- days.years   ## OLD: days

  Linfs <- seq(range.Linf[1],range.Linf[2],step.Linf)
  Ks <- seq(range.K[1],range.K[2],step.K)

  #t <- t.days
  t <- t.years


  ESP.tshift.L <- list()
  ESP.list.L <- list()
  for(li in 1:length(Linfs)){

    ESP.tshift.k <- list()
    ESP.list.k <- list()
    for(ki in 1:length(Ks)){

      # VBGF
      Lt <-  Linfs[li] * ( 1 - exp( - Ks[ki] * (t - t0)))

      # get lengths which fall into sampling times
      samp.times <- as.numeric(colnames(catch.aAF))   # what about the next cohort one year later or one year ago if totap period is over one year???

      lt.classes.list <- list()
      for(tshift in 1:(length(t.days)-1)){
        ttx <- t.years[tshift]
        diff <- samp.times - ttx
        diff[which(diff < 0)] <- NA
        lt.classes <- Lt[(diff*365)]
        lt.classes.list[[tshift]] <- data.frame(tshift = rep(tshift,length(lt.classes)),
                                                lt.classes)
      }


      # get size classes in which calculated Lt values fall into
      for(hits in 1:length(lt.classes.list)){
        df <- lt.classes.list[[hits]]
        for(hiti in 1:length(df[,1])){
          if(all(is.na(classes - df[,2][hiti]))){df$cor.class[hiti] <- NA
          }else df$cor.class[hiti] <- which.min(abs(classes - df[,2][hiti]))
        }
        lt.classes.list[[hits]] <- df
      }


      # get ESP values from catch.aAF matrix entries which are hit by growth curve
      ESP.list <- list()
      for(inds in 1:length(lt.classes.list)){
        df <- lt.classes.list[[inds]]
        ESPs <- NA
        for(indi in 1:length(df[,1])){
          ESPs[indi] <- catch.aAF[df$cor.class[indi],indi]
        }
        # get ESP by summing up
        ESP.list[[inds]] <- c(inds,sum(ESPs,na.rm=TRUE))
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



#   tvar.list <- list()
#   for(tvar in 0:length(t.days)){
#     Lt <-  Linfs[li] * ( 1 - exp( - K[ki] * ((t+tvar) - t0)))
#     tvar.list[[tvar+1]] <- Lt
#   }


  #for plotting
  xplot <- c(min(t.years.fr),(max(t.years.fr)))
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
  y.b <- classes - interval / 2
  y.t <- classes + interval / 2
  for(histi in 1:dim(catch)[2]){
    x.l <- rep(days.years[histi],length(classes))
    x.r <- x.l + catch.aAF[,histi] * 0.03   ### make this number dependent on data!?!?
    rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t,
         col = "gray80", border = "gray40")
  }

  # growth curves
  Lt <- 11 * (1 - exp(-0.9 * (t.years.fr - 0)))    ### choose correct t -> tshift
  lines(y = Lt, x = t.years.fr, lty=2, col='blue')
  #t.years.frA <- seq(xplot[1],xplot[2],1/365)
  Lt <- 11 * (1 - exp(-0.9 * (t.years.fr + 1)))
  lines(y = Lt, x = t.years.fr, lty=2, col='blue')
  Lt <- 11 * (1 - exp(-0.9 * (t.years.fr - 1)))
  lines(y = Lt, x = t.years.fr , lty=2, col='blue')




#   #start with seed values for L8 and K and create growth functions from every peak
#   peaks <- which(!is.na(pos.NA.df))
#   #project backward and forward in time? wtf
}


# works:

lt.classes.list <- list()
for(tshift in 0:(length(t.days)-1)){
  diff <- samp.times - tshift
  diff[which(diff < 0)] <- NA
  lt.classes <- Lt[diff]
  lt.classes.list[[tshift+1]] <- data.frame(tshift = rep(tshift,length(lt.classes)),
                                            lt.classes)
}
#lt.classes.df <- do.call(rbind,lt.classes.list)



ESP.list.L <- list()
for(li in 1:length(Linfs)){
  ESP.list.k <- list()
  for(ki in 1:length(Ks)){

    # VBGF
    Lt <-  Linfs[li] * ( 1 - exp( - Ks[ki] * (t - t0)))

    # get lengths which fall into sampling times
    lt.classes <- Lt[which(colnames(catch.aAF) %in% order(Lt))]

    # get size classes in which calculated Lt values fall into
    hit.list <- NA
    for(hits in 1:length(lt.classes)){
      hit.list[hits] <- which.min(abs(classes - lt.classes[hits]))
    }

    # get ESP values from catch.aAF matrix entries which are hit by growth curve
    ESPs <- NA
    for(indi in 1:length(hit.list)){
      ESPs[indi] <- catch.aAF[hit.list[indi],indi]
    }
    # get ESP by summing up
    ESP <- sum(ESPs)

    ESP.list.k[[ki]] <- c(Linfs[li],Ks[ki],ESP)
  }
  ESP.list.L[[li]] <- do.call(rbind,ESP.list.k)
}

ESP.df <- do.call(rbind,ESP.list.L)
tapply(ESP.df[,3],list(ESP.df[,1],ESP.df[,2]),round,digits = 2)

#for plotting
xplot <- c(0,(max(days)+50))
yplot <- c(classes[1],classes[length(classes)])

# Plot with rearranged histogramms
par(mar = c(5, 5, 1, 1) + .1)
plot(xplot, yplot, type = "n",axes = FALSE,
     ann = FALSE,  xaxs = "i", yaxs = "i")
text(days, par("usr")[3] - 0.7,
     labels = dates.all, srt = 45, pos = 1, xpd = TRUE)
axis(2, cex.axis = 1.2)
mtext(side = 2, outer = F, line = 3.5, "Length [cm]", cex = 1.2)
box( col = "gray40") #bty = "L"

#Histograms
y.b <- classes - interval / 2
y.t <- classes + interval / 2
for(histi in 1:dim(catch)[2]){
  x.l <- rep(days[histi],length(classes))
  x.r <- x.l + catch.aAF[,histi] * 10      ### make this number dependent on data!?!?
  rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t,
       col = "gray80", border = "gray40")
}

Lt <- 10 * (1 - exp(-0.1 * (t - 80)))    ### choose correct t -> tshift
lines(y = Lt, x = t+320, lty=2, col='blue')


# Plot with real histogramms
par(mar = c(5, 5, 1, 1) + .1)
plot(xplot, yplot, type = "n",axes = FALSE,
     ann = FALSE,  xaxs = "i", yaxs = "i")
text(days, par("usr")[3] - 0.7,
     labels = dates.all, srt = 45, pos = 1, xpd = TRUE)
axis(2, cex.axis = 1.2)
mtext(side = 2, outer = F, line = 3.5, "Length [cm]", cex = 1.2)
box( col = "gray40") #bty = "L"

#Histograms
y.b <- classes - interval / 2
y.t <- classes + interval / 2
for(histi in 1:dim(catch)[2]){
  x.l <- rep(days[histi],length(classes))
  x.r <- x.l + catch[,histi] * 0.6     ### make this number dependent on data!?!?
  rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t,
       col = "gray80", border = "gray40")
}

