#' @title ELEFAN
#'
#' @description Electronic length frequency analysis
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

ELEFAN <- function(param, srt.Linf, srt.K){

  res <- param
  classes <- res$midLengths
  catch <- res$catch


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
  names(data.aAF)[3:length(names(data.aAF))] <- unlist(strsplit(names(data)[3:length(names(data))],split = 'X'))[seq(2,((length(names(data))-2)*2),2)]
  names(data)[3:length(names(data))] <- unlist(strsplit(names(data)[3:length(names(data))],split = 'X'))[seq(2,((length(names(data))-2)*2),2)]
  #bring sample dates in right format
  dates <- as.POSIXlt(as.Date(names(data)[3:length(names(data))],format = "%d.%m.%Y"))
  #add same samples one year before and after
  dates.before <- dates
  dates.after <- dates
  dates.before$year <- dates$year - 1
  dates.after$year <- dates$year + 1
  dates.all <- c(dates.before,dates,dates.after)
  dates.all <- as.Date(dates.all)
  #get continuous time in years
  days <- as.numeric(dates.all - as.Date(cut(dates.all[1], "year")))
  t.days <- seq(1,days[length(days)],1)
  t.years <- (t.days)/365

  #start with seed values for L8 and K and create growth functions from every peak
  peaks <- which(!is.na(pos.NA.df[,3:length(names(pos.NA.df))]))
  #project backward and forward in time? wtf


  L8 = 3:30
  K = seq(0.1,0.9,0.1)
  t0 = seq(0,1,0.1)
  vec.Lt <- list()
  vec.K <- list()
  vec.t0 <- list()
  for(i in 1:length(L8)){
    for(j in 1:length(K)){
      for(k in 1:length(t0)){
        Lt  =  L8[i] * (1 - exp(-K[j] * (t.years - t0[k])))
        vec.t0[[k]] <- Lt
      }

      vec.K[[j]] <- vec.t0
    }

    vec.Lt[[i]] <- vec.K
  }

  plot(vec.Lt[[length(vec.Lt)]][[length(vec.K)]][[length(vec.t0)]] ~ t.years,type='l',col='blue',lwd=2,
       ylim=c(0,20))
  abline(v=((days - 4)) / 365)
  for(i in 1:(length(vec.Lt)-1)){
    for(j in 2:length(vec.Lt[[1]])){
      for(k in 2:length(vec.Lt[[1]][[1]])){
        lines(vec.Lt[[i]][[j]][[k]] ~ t.years,col='blue')
      }
    }
  }





}




