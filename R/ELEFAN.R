#' @title ELEFAN
#'
#' @description Electronic length frequency analysis for estimating growth parameter of exploited fish populations.
#'
#' @param param A list consisting of following parameters:
#' \itemize{
#'   \item \strong{age} or \strong{midLengths} midpoints of the length class as vector (length frequency
#'   data) or ages as vector (age composition data),
#'   \item \strong{catch} Catch as vector, or a matrix with catches of subsequent years if
#'   the catch curve with constat time intervals should be applied;
#' }
#' @param range.Linf lower and upper limit of range of Linf values
#' @param step.Linf step size between Linf values
#' @param range.K lower and upper limit of range of K values
#' @param step.K step size between K values
#' @param t0 theoretical age of fish at length zero
#' @param tmax maximal age of fish from literature only estimate needed
#' @param C growth oscillation amplitude
#' @param WP winter point
#' @param conf.int logical; if \code{TRUE} calculates the 95% confidence intervals for
#'   calculated Linf and K parameters by means of the jackknife method (Quenouille, 1956; Tukey, 1958, 1986),
#'
#' @examples
#' \donttest{
#' data(trout)
#' ELEFAN(trout, range.Linf = c(11,15), step.Linf = 1, range.K = c(0.2,1),
#'    step.K = 0.1, tmax = 5)
#'
#' ELEFAN(param = trout,
#'  range.Linf = c(14,18),
#'  step.Linf = 1,
#'  range.K = c(0.4,1.2),
#'  step.K = 0.1,
#'  t0 = 0,
#'  growth.curve = "VBGF",
#'  tmax = 5,
#'  conf.int = FALSE)
#'
#'
#' # FOR CODING:
#' param = trout
#' range.Linf = c(11,15)
#' step.Linf = 1
#' range.K = c(0.2,1)
#' step.K = 0.1
#' t0 = 0
#' C = 0
#' WP = 0
#' growth.curve = "VBGF"
#' tmax = 5
#' conf.int = FALSE
#' }
#'
#' @details Winter point (WP) designates the period of the year when growth is slowest, which
#'    in the northern hemisphere is often found around 0.2 (February) and in the southern
#'    hemisphere aorund 0.7 (October). C express the amplitude of growth oscillations and
#'    should be between 0 (no oscillation) and 1 (no growth), values above 1 suggest
#'    a reduction of length and thus are unlikely (Pauly and Morgan, 1987).
#'
#' @references
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
#' Pauly, D. and G.R.Morgan (Eds.), 1987. Length-based methods in fisheries research.
#' (No. 13). WorldFish
#'
#' Quenouille, M. H., 1956. Notes on bias in estimation. \emph{Biometrika}, 43:353-360
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2): 407 p.
#'
#' Tukey, J., 1958. Bias and confidence in not quite large samples.
#' \emph{Annals of Mathematical Statistics}, 29: 614
#'
#' Tukey, J., 1986. The future of processes of data analysis. In L. V. Jones (Eds.),
#' The Collected Works of John W. Tukey- philosophy and principles of data analysis:
#' 1965-1986 (Vol. 4, pp. 517-549). Monterey, CA, USA: Wadsworth & Brooks/Cole
#'
#' @export
#'

ELEFAN <- function(param, range.Linf, step.Linf,
                   range.K, step.K, t0 = 0,
                   tmax, growth.curve = "VBGF", seed.C = 0, seed.WP = 0, #default of WP? as long as C is 0 it doesnt matter
                   conf.int = FALSE){

  res <- param
  classes <- res$midLengths
  catch <- res$catch
  interval <- classes[1]-classes[2]


  #________________________DATA ARRANGEMENT________________________#
  #  ELEFAN 0

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
#-------------

  #__________________________GROWTH CURVE__________________________#
  #  ELEFAN I

#-------------
  #problem of arrangement of t, in literature: samples sequentially arranged in time
  #display months on the x axis

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

  #for whole live of fish:
  t <- 0:(tmax * 365) / 365   #t <- t.years   ## for whoel range of t max?

  # get other growth curves going through data which is dependent on tmax
  if(sp.years > 0 & sp.years < 1) repro.add <- 0
  if(sp.years > 1 & sp.years < 2) repro.add <- 1
  if(sp.years > 2 & sp.years < 3) repro.add <- 2
  if(sp.years > 3 & sp.years < 4) repro.add <- 3
  if(sp.years > 4 & sp.years < 5) repro.add <- 4

  samp.times <- as.numeric(colnames(catch.aAF))   # what about the next cohort one year later or one year ago if totap period is over one year???



  #for CHECKING run time of big loop
  ##### sys_timeBF <- Sys.time()


  ESP.tshift.L <- list()
  ESP.list.L <- list()
  for(li in 1:length(Linfs)){

    ESP.tshift.k <- list()
    ESP.list.k <- list()
    for(ki in 1:length(Ks)){

      # VBGF
      switch(growth.curve,

             "VBGF.old" ={
               Lt <- Linfs[li] * (1 - exp(-Ks[ki] * (t - t0)))
             },

             "VBGF.seasonalized" ={
               Lt <- Linfs[li] * (1 - exp(-Ks[ki] * D * (t - t0) +
                                            ((C*Ks[ki]*D)/(2*pi)) *
                                            sin(2*pi) * (t - (WP-0.5))) ^ (1/D))
               # D not calculatable with length-frequency data
               # empirical formula for coarse estimate:
               # D = 3 * (1 - (0.6742 + 0.03574 * log10(Wmax)))
             },

             "VBGF" ={   # K & C & WP
               Lt <- Linfs[li] * (1 - exp(-Ks[ki] * (t - t0) +
                                            ((C * Ks[ki])/(2*pi)) *
                                            sin(2*pi) * (t - (WP-0.5))))
             },

             "Gompertz" ={  # Winsor (1932)   # B & C
               Lt <- Linfs[li] * exp(-B * exp(-C * t))
             },

             "logistic" ={  # Verhulst (1838)   # b & K
               Lt <- Linfs[li] / (1 + b * exp(-K * t))
             },

             "Krueger" ={   # N & r
               Lt <- Linfs[li] / (N^((t + r)^-1))
             })

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

          # get corresponding length groups instead of exact lengths, IMPROVE smaller than smallest group, larger than largest group
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
  #old:
  #ESP.df[,3] <- ESP.df[,3]/ASP
  #new: (Pauly, 1987)
  ESP.df[,3] <- (10^(ESP.df[,3]/ASP)) /10
  score_mat <- tapply(ESP.df[,3],list(ESP.df[,1],ESP.df[,2]),round,digits = 2)

  ESP.time <- do.call(rbind,ESP.tshift.L)
  time_mat <- tapply(ESP.time[,3],list(ESP.time[,1],ESP.time[,2]),round,digits = 2)
#-------------

  # calculation of confidence interval by use of jackknife method
  if(conf.int == TRUE){
    for(CIs in 1:dim(catch.aAF)[2]){
      catch.aAF.CI <- catch.aAF[,-CIs]
      ### HOW TO GET BEST FIT OUT OF SURFACE RESPONSE MATRIX?? IT WAS SUBJECTIVE CHOICE
    }


  }

  #for CHECKING run time of big loop
  #### sys_timeAW <- Sys.time()
  #### sys_timeAW - sys_timeBF

  #____________________RESPONSE SURFACE ANALYSIS___________________#

  image(x = Linfs,
      y = Ks,
      z = score_mat, col=colorRampPalette(c("yellow","red"), space="Lab")(5),
      main = 'Response surface analysis', ylab = 'K', xlab='Linf')
  #grid (NULL,NULL, lty = 6, col = "cornsilk2")
  text(x=ESP.df[,1],y=ESP.df[,2],round(ESP.df[,3],digits = 2),cex = 0.6)

  #identify(x = Linfs,y = Ks)


  res2 <- list(score_mat = score_mat,
               time_mat = time_mat)
  ret <- c(res,res2)
  return(ret)
  class(ret) <- 'ELEFAN'
}


