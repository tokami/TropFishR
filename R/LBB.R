#' @title Length-based Bayesian biomass estimator (LBB)
#'
#' @description A new approach for estimating stock status from length-frequency
#'    data 
#'
#' @param lfq A list of the class "lfq" consisting of following parameters:
#' \itemize{
#'   \item \strong{species} species name,
#'   \item \strong{stock} stock ID or name,
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row)
#'      and sampling date (column),
#'   \item \strong{comments} comments;
#' }
#' @param startYear Start year of assessment. If NA (default), the first year in the
#'    \code{lfq$dates} is used.
#' @param endYear Final year of assessment. If NA (default), the last year in the
#'    \code{lfq$dates} is used.
#' @param years Manual selection of years for assessment. If NA (default), all years
#'    in the \code{lfq$dates} are used.
#' @param binSize Optional; determines bin size (class width) for length-frequency data.
#'    If NA (default) bin size remains unchanged (as in \code{lfq$midLengths}).
#' @param LinfUser Optional; user defined asymptotic length. Any length observation larger than
#'    this value will be removed from the data. If NA (default), Linf is estimated by
#'    the model.
#' @param LcutUser Optional, user defined minimum length. Any length observation smaller than
#'    this value will be removed form the data. By default all length obervations are used
#'    (\code{LcutUser} = NA).
#' @param LcUser Optional, user defined selectivity parameter. If NA (default) L10 and L90 are
#'    used to estimate a proxy for Lc.
#' @param LstartUser Optional, user defined length at which selectivity is 0.95. If NA (default)
#'    Lstart (L95) is estimated by the model.
#' @param MKUser Optional; user defined MK ratio. If NA (default) MK ratio is set to 1.5.
#' @param mmUser Logical; indicating the unit of length measurements, where TRUE
#'    indicates that lengths are in mm and FALSE (default) indicate that lengths are in cm.
#' @param GausSel Logical; indicating the selectivity pattern. If FALSE (default) trawl-like,
#'    if TRUE gaussian selectivity is assumed.
#' @param MergeLF Logical; indicating if the data of subsequent years should be merged
#'    with data of preceeding years. (Default: FALSE).
#' @param plot Logical; should the individual year plot be displayed? (Default: FALSE).
#' @param mfrow A vector of the form ‘c(nr, nc)’.  Subsequent figures will be drawn in an
#'    ‘nr’-by-‘nc’ array on the device by _rows_ (‘mfrow’). If NA (default), a panel with
#'    3 columns and several rows (dependent on number of years) is used.
#'
#' @details Requires the Gibbs sampler JAGS to be installed on your computer,
#'   available for your Operating System from the following web
#'   site: \link{http://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/}.
#'   LBB is a new method for the analysis of length frequency data from the
#'   commercial fishery. It works for species that grow throughout their lives,
#'   such as most commercial fish and invertebrates, and requires no input in
#'   addition to length frequency data. It estimates asymptotic length (Linf),
#'   length at first capture (Lc), relative natural mortality (M/K) and relative
#'   fishing mortality (F/M) as means over the age range represented in the
#'   length-frequency sample. With these parameters as input, standard fisheries
#'   equations can be used to estimate depletion or current exploited biomass
#'   relative to unexploited biomass (B/B0). In addition, these parameters allow
#'   the estimation of the length at first capture that would maximize catch and
#'   biomass for the given fishing effort (Lc_opt), and estimation of a proxy
#'   for the relative biomass capable of producing maximum sustainable yields
#'   (Bmsy/B0). Relative biomass estimates of LBB were not significantly
#'   different from the “true” values in simulated data and similar to independent
#'   estimates from full stock assessments.
#' 
#' @author Rainer Froese, (\email{rfroese@geomar.de})
#' 
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{GausSel}: indicating if gaussian Selection was used,
#'   \item \strong{priors}: priors,
#'   \item \strong{refLev}: matrix with all reference levels for all years,
#'   \item \strong{medianRefLev}: median reference levels (plus 95\% confidence intervals),
#'   \item \strong{lastRefLev}: reference levels in the last year,
#'   \item \strong{LFall}: matrix with lengths and relative frequencies.
#' }
#' 
#' @examples
#' \donttest{
#' ## load data
#' data(synLFQ8)
#'
#' ## arrange lfq data
#' lfq <- lfqModify(synLFQ8, aggregate = "year")
#'
#' ## plot lfq data traditionally
#' plot(lfq)
#'
#' ## plot data in LBB manner
#' plotLBB.data(lfq)
#'
#' ## add length at maturity to lfq data
#' lfq$Lm50 <- 36
#'
#' ## run LBB model
#' res <- LBB(lfq, plot = TRUE)
#'
#' ## plot results
#' plotLBB(res)
#'
#' ## plot time series
#' plotLBB.ts(res)
#' }
#'
#' @import R2jags
#' @importFrom Hmisc wtd.quantile
#'
#' @references
#' R. Froese, H. Winker, G. Coro, N. Demirel, A.C. Tsikliras, D. Dimarchopoulou,
#' G. Scarcella, W.N. Probst, M. Dureuil, and D. Pauly (2018) A new approach
#' for estimating stock status from length frequency data. ICES Journal of Marine Science. DOI: 10.1093/icesjms/fsy078
#'
#' @export
#' 
LBB <- function(lfq, startYear=NA, endYear=NA, years=NA, binSize=NA, LinfUser=NA, LcutUser=NA,
                LcUser=NA, LstartUser=NA, MKUser=NA, mmUser=FALSE, GausSel=FALSE, MergeLF=FALSE,
                plot=FALSE, mfrow = NA){

    ## informative warning messages if input not correct
    ## length at maturity
    if(!"Lm50" %in% names(lfq)){
        stop("Lm50 missing! Please provide the length at maturity as element 'Lm50' in the argument 'lfq'.")
    }

    ## make selection
    startDate <- endDate <- NA
    if(!is.na(startYear)) startDate <- as.Date(paste0(startYear,"-01-01"))
    if(!is.na(endYear)) endDate <- as.Date(paste0(endYear,"-01-01"))
    res <- lfqModify(lfq, minDate = startDate, maxDate=endDate, years = years,
                     bin_size=binSize, Lmax=LinfUser, Lmin=LcutUser)

    years <- as.numeric(format(lfq$dates, "%Y"))

    ##--------------------------------------------------------------------------------------    
    ## Put data into vectors
    ##--------------------------------------------------------------------------------------
    if(class(lfq$catch) == "matrix"){
        nrowC <- nrow(lfq$catch)
        ncolC <- ncol(lfq$catch)
    }else if(class(lfq$catch) == "numeric"){
        nrowC <- length(lfq$catch)
        ncolC <- 1
    }else{
        stop("lfq$catch has to be either a matrix or a numeric vector!")
    }
    
    StartYear  <- min(years)
    EndYear    <- max(years)
    AllYear    <- rep(years, each = nrowC)
    AllLength  <- rep(lfq$midLengths, ncolC)
    AllFreq    <- as.numeric(lfq$catch)
    Years      <- sort(unique(AllYear))
    nYears     <- length(years)
    Stock      <- ifelse(!is.null(res$stock), res$stock, "unknown stock")
    Comment    <- ifelse(!is.null(res$comment), res$comment, "no comment")
    Species    <- ifelse(!is.null(res$species), res$species, "unknown species")    

    ## transformation needed from lfq object with zero length classes to NA length classes
    idx <- which(AllFreq != 0)
    AllLength <- AllLength[idx]
    AllYear <- AllYear[idx]
    AllFreq <- AllFreq[idx]

    ##--------------------------------------------------------------------------------------
    ## Create matrix to store annual estimates
    ##--------------------------------------------------------------------------------------    
    Ldat <- data.frame(Stock=rep(Stock,nYears),Year=rep(NA,nYears),
                       Linf=rep(NA,nYears),
                       Linf.lcl=rep(NA,nYears),
                       Linf.ucl=rep(NA,nYears),
                       Lc=rep(NA,nYears), # for trawl selection
                       Lc.lcl=rep(NA,nYears),
                       Lc.ucl=rep(NA,nYears),
                       Lmean=rep(NA,nYears),
                       r.alpha=rep(NA,nYears),
                       r.alpha.lcl=rep(NA,nYears),
                       r.alpha.ucl=rep(NA,nYears),
                       r.GLmean=rep(NA,nYears),r.SD=rep(NA,nYears), # for gill net selection
                       MK=rep(NA,nYears),
                       MK.lcl=rep(NA,nYears),
                       MK.ucl=rep(NA,nYears),
                       FK=rep(NA,nYears),
                       FK.lcl=rep(NA,nYears),
                       FK.ucl=rep(NA,nYears),
                       ZK=rep(NA,nYears),
                       ZK.lcl=rep(NA,nYears),
                       ZK.ucl=rep(NA,nYears),
                       FM=rep(NA,nYears),
                       FM.lcl=rep(NA,nYears),
                       FM.ucl=rep(NA,nYears),
                       r.Lopt=rep(NA,nYears),
                       BB0=rep(NA,nYears),
                       BB0.lcl=rep(NA,nYears),
                       BB0.ucl=rep(NA,nYears),
                       YR=rep(NA,nYears),
                       YR.lcl=rep(NA,nYears),
                       YR.ucl=rep(NA,nYears),
                       perc.mat=rep(NA,nYears),
                       L95=rep(NA,nYears))

    ##--------------------------------------------------------------------------------------
    ## Use aggregated LF data for estimation of Linf (and overall Z/K)
    ##--------------------------------------------------------------------------------------    
    df        <- data.frame(AllYear,AllLength,AllFreq)
    names(df) <- c("Year","Length","Freq")

    LF.all    <- AG(dat=df) # function to aggregate data by year and across years

    # standardize to max Freq
    LF.all$Freq = LF.all$Freq/max(LF.all$Freq) 
    # remove leading empty records
    LF.all     <- LF.all[which(LF.all$Freq>0)[1] : length(LF.all$Length),]
    # remove trailing empty records
    LF.all     <- LF.all[1 : which(LF.all$Length==max(LF.all$Length[LF.all$Freq>0])),]

    # get number of records in LF.all
    n.LF.all   <- length(LF.all$Length) 
    # use largest fish as Lmax
    Lmax       <- LF.all$Length[n.LF.all]
    # use median of largest fish per year as Lmax.med
    Lmax.med   <- median(as.numeric(by(rep(res$midLengths,ncolC)[as.numeric(res$catch)>0],
                                       rep(res$dates,each=nrowC)[as.numeric(res$catch)>0],max)))/10    

    # If no Linf is provided by the user (preferred), determine Linf from fully selected LF:
    # Freq=Nstart*exp(ZK*(log(1-L/Linf)-log(1-Lstart/Linf)))
    # Nstart is canceled out when dividing both sides by their sums
    # ---------------------------------------------------------
    # determine start values of selection ogive to find first fully selected length class Lstart
    L10 <- LF.all$Length[which(LF.all$Freq>0.1)[1]] # use length at 10% of peak frequency as proxy for L10
    L90 <- LF.all$Length[which(LF.all$Freq>0.9)[1]] # use length at 90% of peak frequency as proxy for L90
    Lc.st <- ifelse(is.na(LcUser),(L10 + L90)/2,LcUser)  # use mean of L10 and L90 as proxy for Lc, else user input
    alpha.st <- -log(1/LF.all$Freq[which(LF.all$Freq>0.1)[1]])/(L10-Lc.st) # use rearranged logistic curve to estimate slope alpha

    ## determine start values for Linf and Z/K 
    Linf.st     <- max(LF.all$Length) # use Lmax as proxy for Linf
    Lmean.st    <- sum(LF.all$Length[LF.all$Length>=Lc.st]*LF.all$Freq[LF.all$Length>=Lc.st])/
                         sum(LF.all$Freq[LF.all$Length>=Lc.st])
    MK.st       <- ifelse(is.na(MKUser), 1.5,MKUser) # default 1.5
    ZK.st       <- (Linf.st-Lmean.st)/(Lmean.st-Lc.st)       # the Holt equation
    FK.st       <- ifelse((ZK.st-MK.st)>0,ZK.st-MK.st,0.3)   # prevent M/K being larger than Z/K

    # get vectors with fully selected length classes for Linf estimation
    if(!is.na(LstartUser)) {
        Lstart <- LstartUser
    } else {
        Lstart     <- (alpha.st*Lc.st-log(1/0.95-1))/alpha.st   # Length where selection probability is 0.95  
        # test if there are enough (>=4) length classes for estimation of aggregated Linf and ZK 
        Lstart.i   <- which(LF.all>=Lstart)[1]
        Lmax.i     <- length(LF.all$Length)
        peak.i     <- which.max(LF.all$Freq)
        if(Lstart.i<(peak.i+1)) Lstart <- LF.all$Length[peak.i+1] # make sure fully selected length starts after peak 
        if((Lmax.i-Lstart.i)<4) Lstart <- LF.all$Length[Lstart.i-1] # make sure enough length classes are available
    }
    
    # do not include Lmax to allow Linf < Lmax and to avoid error in nls when Linf-L becomes negative
    L.L         <- LF.all$Length[LF.all$Length >= Lstart  & LF.all$Length < Linf.st]
    L.Freq      <- LF.all$Freq[LF.all$Length>=L.L[1]& LF.all$Length < Linf.st]

    ## plottting
    if(length(L.L)<4){
        plot(x=LF.all$Length,y=LF.all$Freq, bty="l",main=Stock)
        lines(x=c(Lstart,Lstart),y=c(0,0.9*max(LF.all$Freq)),lty="dashed")
        text(x=Lstart,y=max(LF.all$Freq),"Lstart")
        lines(x=c(Linf.st,Linf.st),y=c(0,0.9*max(LF.all$Freq)),lty="dashed")
        text(x=Linf.st,y=max(LF.all$Freq),"Lmax")
        stop("Too few fully selected data points: set Lstart.user\n")
    }


    ## standardize frequencies by dividing by sum of observed frequencies, needed to drop NLstart from equation
    sum.L.Freq  <- sum(L.Freq)
    L.Freq      <- L.Freq/sum.L.Freq

    # use nls() to find Linf-ZK combination with least residuals
    if(is.na(LinfUser)){
      Linf.mod    <- nls(L.Freq ~ ((Linf-L.L)/(Linf-Lstart))^ZK /
                        sum(((Linf-L.L)/(Linf-Lstart))^ZK),
                      start=list(ZK=ZK.st,Linf=Linf.st),
                      lower=c(0.5*ZK.st,0.999*Linf.st), 
                      upper=c(1.5*ZK.st,1.2*Linf.st), 
                      algorithm = "port")
      ZK.nls       <- as.numeric(coef(Linf.mod)[1])
      ZK.nls.sd    <- as.numeric(coef(summary(Linf.mod))[,2][1])
      ZK.nls.lcl   <- ZK.nls-1.96*ZK.nls.sd
      ZK.nls.ucl   <- ZK.nls+1.96*ZK.nls.sd
      Linf.nls     <- as.numeric(coef(Linf.mod)[2])
      Linf.nls.sd  <- as.numeric(coef(summary(Linf.mod))[,2][2])
      Linf.lcl     <- Linf.nls-1.96*Linf.nls.sd
      Linf.ucl     <- Linf.nls+1.96*Linf.nls.sd
    } else {  # end of loop to determine Linf and ZK.L
              # use given Linf and determine ZK.L
      # use Linf provided by user if given
       Linf.nls    <- LinfUser 
       Linf.nls.sd <- 0.01*LinfUser
       ZK.mod      <- nls(L.Freq ~ exp(ZK*(log(1-L.L/Linf.nls)-log(1-L.L[1]/Linf.nls)))/
                           sum(exp(ZK*(log(1-L.L/Linf.nls)-log(1-L.L[1]/Linf.nls)))),
                         start=list(ZK=ZK.st),
                         lower=c(0.7*ZK.st),
                         upper=c(1.3*ZK.st), 
                         algorithm = "port")
       ZK.nls       <- as.numeric(coef(ZK.mod)[1])
       ZK.nls.sd    <- as.numeric(coef(summary(ZK.mod))[,2][1])
       ZK.nls.lcl   <- ZK.nls-1.96*ZK.nls.sd
       ZK.nls.ucl   <- ZK.nls+1.96*ZK.nls.sd
   } # end of loop if Linf is given by user




    ## get vector of all lengths <= prior Linf to avoid error in equation
    AllFreq       <- AllFreq[AllLength <= Linf.nls]
    AllYear       <- AllYear[AllLength <= Linf.nls]
    AllLength     <- AllLength[AllLength <= Linf.nls]



    
    #-----------------------------------------
    # Start LF analysis by year
    #-----------------------------------------
    cat("Running Jags model to fit SL and N distributions for",Species,"\nin", Years,"....\n")
    i = 0 # start counter

    for(Year in Years){

        i = i+1
        # if MergeLF==TRUE and if this is the second or heigher year and no simulation, aggregate LF with previous year LF
        if(i>1 & MergeLF & substr(Stock,start=nchar(Stock)-2,stop=nchar(Stock))!="Sim"){
            AG.yr <- c(Years[i-1],Year)
        } else AG.yr <- Year

        # aggregate data within the year (sometimes there are more than one sample per year)
        df        <- data.frame(AllYear[AllYear%in%AG.yr],
                                AllLength[AllYear%in%AG.yr],
                                AllFreq[AllYear%in%AG.yr])
        names(df) <- c("Year","Length","Freq")
        LF.y      <- AG(dat=df) # function to aggregate data by year and across years
        LF.y$Freq <- LF.y$Freq/sum(LF.y$Freq) # standardize frequencies

        # remove empty leading and trailing records
        LF.y        <- LF.y[which(LF.y$Freq>0)[1] : length(LF.y$Length),]
        LF.y        <- LF.y[1 : which.max(LF.y$Length[LF.y$Freq>0]),]
        # get vectors
        L.y         <- LF.y$Length
        r.Freq.y    <- LF.y$Freq

        # fill remaining zero frequencies with very small number, to avoid error
        r.Freq.y[r.Freq.y==0] <- min(r.Freq.y[r.Freq.y>0],na.rm=T)/100
        # enter data for this year into data frame
        Ldat$Year[i]     <- Year

        ##-------------------------------------------------------------------------
        ## Estimate annual parameters Lc, alpha, M/K, F/K from LF curve with trawl-type selection
        ##-------------------------------------------------------------------------
        # determine priors 
        n.L         <- length(L.y)
        Linf.pr     <- Linf.nls
        Linf.sd.pr  <- ifelse(Linf.nls.sd/Linf.nls<0.01,Linf.nls.sd,0.01*Linf.nls) # restict prior CV of Linf to < 0.01
        MK.pr       <- MK.st
        MK.sd.pr    <- ifelse(is.na(MKUser)==TRUE,0.15,0.075)


        if(GausSel==FALSE){ # apply trawl-like selection 

            ## with 1.02 multiplier to account for systematic small underestimation
            Lc.pr <- ifelse(is.na(LcUser)==TRUE,1.02*Lc.st,LcUser) 
            ## assume narrower SD if Lc is given by user
            Lc.sd.pr <- ifelse(is.na(LcUser)==TRUE,0.1*Lc.pr,0.05*Lc.pr) 
            r.max.Freq <- max(r.Freq.y,na.rm=T) 
            r.alpha.pr <- -log(r.max.Freq/r.Freq.y[which(r.Freq.y>(0.1*r.max.Freq))[1]])/
                (L10/Linf.nls-Lc.st/Linf.nls) # relative alpha for standardized data
            r.alpha.sd.pr <- 0.025*r.alpha.pr 
            FK.pr <- ifelse((ZK.nls-MK.st) > 0,ZK.nls-MK.st,0.3) # if Z/K <= M/K assume low F/K = 0.3 

            # list of data to pass to JAGS plus list of parameters to estimate   
            jags.data <- list("r.Freq.y","L.y","n.L","Linf.pr","Linf.sd.pr",
                              "Lc.pr","Lc.sd.pr","r.alpha.pr","r.alpha.sd.pr",
                              "MK.pr","MK.sd.pr","FK.pr")
            jags.params <- c("r.alpha.d","Lc.d","SL","xN","FK.d","MK.d","Linf.d")

            ##---------------------------------
            ## LBB JAGS model
            ##---------------------------------
            sink("SLNMod.jags")
            cat("
             model {
             r.alpha.d_tau  <- pow(r.alpha.sd.pr, -2) 
             r.alpha.d      ~ dnorm(r.alpha.pr,r.alpha.d_tau) 

              Lc.d_tau  <- pow(Lc.sd.pr,-2)
              Lc.d      ~ dnorm(Lc.pr,Lc.d_tau) #       

              MK.d_tau  <-pow(MK.sd.pr, -2) # strong prior on M/K
              MK.d      ~ dnorm(MK.pr, MK.d_tau)

              Linf.tau  <- pow(Linf.sd.pr,-2) 
              Linf.d    ~ dnorm(Linf.pr,Linf.tau)

              FK.d       ~ dlnorm(log(FK.pr),4) # wide prior range for F/K

              SL[1]       ~ dlogis(0,1000)
              Freq.pred[1]<-0
              xN[1]       <-1

              for(j in 2:n.L) {
               SL[j]<- 1/(1+exp(-r.alpha.d*(L.y[j]/Linf.d-Lc.d/Linf.d))) # selection at length L[j]

               xN[j] <- xN[j-1]*((Linf.d-L.y[j])/(Linf.d-L.y[j-1]))^(MK.d+FK.d*SL[j])

                           Freq.pred[j]<-xN[j]*SL[j]

               # normalize frequencies by dividing by sum of frequencies; multiply with 10 to avoid small numbers and with 1000 for effective sample size
               r.Freq.pred[j]<- Freq.pred[j]/sum(Freq.pred)*10*1000
             }	

             #><> LIKELIHOOD FUNCTION
             #><> Fit observed to predicted LF data using a Dirichlet distribution (more robust in JAGS)
             r.Freq.y[2:n.L] ~ ddirch(r.Freq.pred[2:n.L])  

             } # END OF MODEL
               ",fill = TRUE)
           sink()

           MODEL = "SLNMod.jags"
           jagsfitSLN <- R2jags::jags.parallel(data=jags.data, working.directory=NULL, inits=NULL, 
                                       parameters.to.save=jags.params, 
                                       model.file=paste(MODEL), 
                                       n.burnin=300, n.thin=10, n.iter=600, n.chains=3)

           # use median and percentiles
           Ldat$Lc[i]      <- median(jagsfitSLN$BUGSoutput$sims.list$Lc.d)
           Ldat$Lc.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Lc.d,0.025)
           Ldat$Lc.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Lc.d,0.975)
           Ldat$Lmean[i]   <- sum(L.y[L.y>=Ldat$Lc[i]]*r.Freq.y[L.y>=Ldat$Lc[i]])/sum(r.Freq.y[L.y>=Ldat$Lc[i]])
           Ldat$r.alpha[i] <- median(jagsfitSLN$BUGSoutput$sims.list$r.alpha.d)
           Ldat$r.alpha.lcl[i]<- quantile(jagsfitSLN$BUGSoutput$sims.list$r.alpha.d,0.025)
           Ldat$r.alpha.ucl[i]<- quantile(jagsfitSLN$BUGSoutput$sims.list$r.alpha.d,0.975)
           Ldat$MK[i]      <- median(jagsfitSLN$BUGSoutput$sims.list$MK.d)
           Ldat$MK.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d,0.025)
           Ldat$MK.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d,0.975)
           Ldat$FK[i]      <- median(jagsfitSLN$BUGSoutput$sims.list$FK.d)
           Ldat$FK.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d,0.025)
           Ldat$FK.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d,0.975)
           FMi             <- jagsfitSLN$BUGSoutput$sims.list$FK.d/jagsfitSLN$BUGSoutput$sims.list$MK.d
           Ldat$FM[i]      <- median(FMi)
           Ldat$FM.lcl[i]  <- quantile(FMi,0.025)
           Ldat$FM.ucl[i]  <- quantile(FMi,0.975)
           ZKi             <- jagsfitSLN$BUGSoutput$sims.list$MK.d + jagsfitSLN$BUGSoutput$sims.list$FK.d
           Ldat$ZK[i]      <- median(ZKi)
           Ldat$ZK.lcl[i]  <- quantile(ZKi,0.025)
           Ldat$ZK.ucl[i]  <- quantile(ZKi,0.975)
           Ldat$r.Lopt[i]  <- 3/(3+Ldat$MK[i])
           Ldat$Linf[i]    <- median((jagsfitSLN$BUGSoutput$sims.list$Linf.d))
           Ldat$Linf.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d,0.025)
           Ldat$Linf.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d,0.975)
        } ## end of trawl-like selection
        


        ##----------------------------------------------------------------------
        ## Estimate parameters GLmean, SD, F/K, M/K if selection is gillnet-like
        ##---------------------------------------------------------------------- 
        if(GausSel==TRUE) {
          # determine priors
          # assume length at peak Freq as mean and distance to length at 80% of peak as SD of mean
          GLmean.st <- L.y[which.max(r.Freq.y)]
          # assume SD of Gaussian selection as distance between length at peak and length at 50% of peak
          Lc.pr     <- L.y[which(r.Freq.y >= (0.5*max(r.Freq.y)))][1]
          SD.st      <- max(GLmean.st-Lc.pr,0.25*GLmean.st)

          cat("Running Jags model to fit SL and N distributions\n")

          n.L <- length(L.y)

          jags.data <- list ("n.L","GLmean.st","L.y","SD.st","ZK.nls","r.Freq.y","Linf.pr","Linf.sd.pr","MK.pr")
          jags.params <- c("GLmean.d","SD.d","SL","xN","FK.d","MK.d","Linf.d")

          #---------------------------
          # JAGS model L-based with integral
          #---------------------------
          sink("SLNMod.jags")
          cat("
              model {
              GLmean.tau <- pow(0.1*GLmean.st,-2) 
              GLmean.d   ~ dnorm(GLmean.st,GLmean.tau)

              SD.tau    <- pow(0.2*SD.st,-2)
              SD.d      ~ dnorm(SD.st,SD.tau)

              MK.d_tau  <-pow(0.15,-2)
              MK.d      ~ dnorm(MK.pr,MK.d_tau)

              Linf.tau  <- pow(Linf.sd.pr,-2)
              Linf.d    ~ dnorm(Linf.pr,Linf.tau)

              FK        <- (ZK.nls-1.5) # ZK overestimated in gillnet selection, used as upper range
              FK.d      ~ dunif(0,FK)  

              SL[1]~ dlogis(0,1000)
              Freq.pred[1]<-0
              xN[1]<-1

              for(j in 2:n.L) {
                SL[j]<- exp(-((L.y[j]-GLmean.d)^2/(2*SD.d^2)))

                xN[j]<-xN[j-1]*exp((MK.d+FK.d*SL[j])*(log(1-L.y[j]/Linf.d)-log(1-L.y[j-1]/Linf.d)))

                Freq.pred[j]<-xN[j]*SL[j]

                #><> add effective sample size (try 100 typical for LF data)
                r.Freq.pred[j]<- Freq.pred[j]/sum(Freq.pred)*10000
              }	

              #><> LIKELIHOOD FUNCTION
              #><> Fit observed to predicted LF data using a Dirichlet distribution (more robust in JAGS)
              r.Freq.y[2:n.L]~ddirch(r.Freq.pred[2:n.L])  

           } # END OF MODEL
              ",fill = TRUE)
          sink()

          MODEL = "SLNMod.jags"
          #jagsfitSLN <- jags(jags.data, inits=NULL, jags.params, paste(MODEL), n.chains = Nchains , n.thin =Nthin , n.iter =Niter , n.burnin = Nburnin)

          jagsfitSLN <- R2jags::jags.parallel(data=jags.data, working.directory=NULL, inits=NULL, 
                                      parameters.to.save=jags.params, 
                                      model.file=paste(MODEL), 
                                      n.burnin=300, n.thin=10, n.iter=1000, n.chains=3)

          # use median and percentiles
          Ldat$GLmean[i]    <- median(jagsfitSLN$BUGSoutput$sims.list$GLmean.d)
          Ldat$GLmean.lcl[i]<- quantile(jagsfitSLN$BUGSoutput$sims.list$GLmean.d,0.025)
          Ldat$GLmean.ucl[i]<- quantile(jagsfitSLN$BUGSoutput$sims.list$GLmean.d,0.975)
          Ldat$SD[i]        <- median(jagsfitSLN$BUGSoutput$sims.list$SD.d)
          Ldat$SD.lcl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$SD.d,0.025)
          Ldat$SD.ucl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$SD.d,0.975)
          Ldat$MK[i]        <- median(jagsfitSLN$BUGSoutput$sims.list$MK.d)
          Ldat$MK.lcl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d,0.025)
          Ldat$MK.ucl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d,0.975)
          Ldat$FK[i]        <- median(jagsfitSLN$BUGSoutput$sims.list$FK.d)
          Ldat$FK.lcl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d,0.025)
          Ldat$FK.ucl[i]    <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d,0.975)
          FMi               <- jagsfitSLN$BUGSoutput$sims.list$FK.d/jagsfitSLN$BUGSoutput$sims.list$MK.d
          Ldat$FM[i]        <- median(FMi)
          Ldat$FM.lcl[i]    <- quantile(FMi,0.025)
          Ldat$FM.ucl[i]    <- quantile(FMi,0.975)
          ZKi               <- jagsfitSLN$BUGSoutput$sims.list$MK.d + jagsfitSLN$BUGSoutput$sims.list$FK.d
          Ldat$ZK[i]        <- median(ZKi)
          Ldat$ZK.lcl[i]    <- quantile(ZKi,0.025)
          Ldat$ZK.ucl[i]    <- quantile(ZKi,0.975)
          Ldat$r.Lopt[i]    <- 3/(3+Ldat$MK[i])
          Ldat$Linf[i]      <- median((jagsfitSLN$BUGSoutput$sims.list$Linf.d))
          Ldat$Linf.lcl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d,0.025)
          Ldat$Linf.ucl[i]  <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d,0.975)

        } # end of gillnet loop



        ## call BH function to estimate B/B0 and YR for the given year [i] 
        BH.list  <- BH(AllLength=unique(AllLength[AllYear==Year]),
                       Linf=Ldat$Linf[i],MK=Ldat$MK[i],FK=Ldat$FK[i],GausSel=GausSel,
                       selpar1=ifelse(GausSel==T,Ldat$GLmean[i]/Ldat$Linf[i],Ldat$Lc[i]/Ldat$Linf[i]),
                       selpar2=ifelse(GausSel==T,Ldat$SD[i]/Ldat$Linf[i],Ldat$r.alpha[i]))

        Ldat$BB0[i]  <- as.numeric(BH.list[1])
        Ldat$YR[i]   <- as.numeric(BH.list[2])


        ## Error propagation, assuming that fractional uncertainties add in quadrature  

        rel.lcl <- sqrt(((Ldat$FM[i]-Ldat$FM.lcl[i])/Ldat$FM[i])^2+
                        ((Ldat$MK[i]-Ldat$MK.lcl[i])/Ldat$MK[i])^2+
                        ((Ldat$FK[i]-Ldat$FK.lcl[i])/Ldat$FK[i])^2+
                        ((Ldat$Linf[i]-Ldat$Linf.lcl[i])/Ldat$Linf[i])^2)
        rel.ucl <- sqrt(((Ldat$FM.ucl[i]-Ldat$FM[i])/Ldat$FM[i])^2+
                        ((Ldat$MK.ucl[i]-Ldat$MK[i])/Ldat$MK[i])^2+
                        ((Ldat$FK.ucl[i]-Ldat$FK[i])/Ldat$FK[i])^2+
                        ((Ldat$Linf.ucl[i]-Ldat$Linf[i])/Ldat$Linf[i])^2)   
        Ldat$BB0.lcl[i] <- Ldat$BB0[i]-Ldat$BB0[i]*rel.lcl
        Ldat$BB0.ucl[i] <- Ldat$BB0[i]+Ldat$BB0[i]*rel.ucl
        Ldat$YR.lcl[i] <- Ldat$YR[i]-Ldat$YR[i]*rel.lcl
        Ldat$YR.ucl[i] <- Ldat$YR[i]+Ldat$YR[i]*rel.ucl

        ## get MSFD D3.3 indicators
        Ldat$L95[i]      <- Hmisc::wtd.quantile(x=L.y,weights=r.Freq.y,probs=c(0.95))
        Ldat$perc.mat[i] <- ifelse(is.na(res$Lm50)==F,sum(r.Freq.y[L.y>res$Lm50])/sum(r.Freq.y),NA)

        ## annual plots
        if(plot){
            ## plot all years
            if(which(Years==Year)==1){
                if(is.na(mfrow)){
                    ncols <- ifelse(length(Years)<=4,2,3)
                    if(length(Years) == 1) ncols <- 1
                    opar <- par(mfrow=c(ceiling(length(Years)/3),ncols))
                }else{
                    opar <- par(mfrow=mfrow)
                }
            }
            r.L.y    <- L.y[L.y < Ldat$Linf[i]] / Ldat$Linf[i] 
            r.Freq.y <- r.Freq.y[L.y < Ldat$Linf[i]]
            plotLBB.year(r.L.y=r.L.y, r.Freq.y=r.Freq.y,r.Lopt=Ldat$r.Lopt[i],
                      SL1=ifelse(GausSel==T,Ldat$GLmean[i],Ldat$Lc[i]),
                      SL2=ifelse(GausSel==T,Ldat$SD[i],Ldat$r.alpha[i]),
                      MK=Ldat$MK[i],FK=Ldat$FK[i],Linf=Ldat$Linf[i],
                      Year = Year, GausSel = GausSel)
        if(which(Years==Year)==length(Years)) par(opar)
        }
    } ## end of annual loop
    

    ## get some reference points as median of time series
    Linf.med     <- median(Ldat$Linf)
    Linf.lcl     <- median(Ldat$Linf.lcl)
    Linf.ucl     <- median(Ldat$Linf.ucl)
    if(GausSel==F) {
      Lc.med       <- median(Ldat$Lc)
      r.alpha.med  <- median(Ldat$r.alpha)
    } else {
      GLmean.med   <- median(Ldat$GLmean)
      SD.med       <- median(Ldat$SD)
    }
    MK.med       <- median(Ldat$MK)
    MK.lcl       <- median(Ldat$MK.lcl)
    MK.ucl       <- median(Ldat$MK.ucl)
    FK.med       <- median(Ldat$FK)
    FK.lcl       <- median(Ldat$FK.lcl)
    FK.ucl       <- median(Ldat$FK.ucl)
    FM.med       <- median(Ldat$FM)
    FM.lcl       <- median(Ldat$FM.lcl)
    FM.ucl       <- median(Ldat$FM.ucl)
    ZK.med       <- median(Ldat$ZK)
    ZK.lcl       <- median(Ldat$ZK.lcl)
    ZK.ucl       <- median(Ldat$ZK.ucl)
    r.Lopt.med   <- median(Ldat$r.Lopt)
    Lopt.med     <- r.Lopt.med*Linf.med
    Lc_opt.med   <- Linf.med*(2+3*FM.med)/((1+FM.med)*(3+MK.med)) 
    BB0.med      <- median(Ldat$BB0)
    BB0.lcl      <- median(Ldat$BB0.lcl)
    BB0.ucl      <- median(Ldat$BB0.ucl)
    YR.med       <- median(Ldat$YR)
    YR.lcl       <- median(Ldat$YR.lcl)
    YR.ucl       <- median(Ldat$YR.ucl)

    BFM1B0.list  <- BH(AllLength=unique(AllLength),Linf=Linf.med,MK=MK.med,FK=MK.med,GausSel=GausSel,
                         selpar1=ifelse(GausSel==T,r.Lopt.med,5/(2*(3+MK.med))),
                         selpar2=ifelse(GausSel==T,SD.med/Linf.med,r.alpha.med))

    BFM1B0       <- as.numeric(BFM1B0.list[1])
    YRFM1        <- as.numeric(BFM1B0.list[2])


    

    ## print results to console
    cat("\n----------------------------------------------------------------------\n")
    cat("Results for",Species,", stock",Stock,",",StartYear,"-",EndYear,ifelse(GausSel==T,", Gaussian selection",""), "\n")
##    cat("(95% confidence limits in parentheses) File:",File,"\n")
    cat("-----------------------------------------------------------------------\n")
    cat("Linf prior =",Linf.pr,", SD =",Linf.sd.pr,"(cm)",ifelse(is.na(LinfUser)==TRUE,"","(user-defined)"),"\n")
    cat("Z/K prior  =",ZK.nls,", SD =", ZK.nls.sd,", M/K prior  =", MK.pr, ", SD =",MK.sd.pr,ifelse(is.na(MKUser)==TRUE,"","(user-defined)"),"\n") 
    if(GausSel==F) { 
      cat("F/K prior  =", FK.pr, "(wide range with tau=4 in log-normal distribution)\n")
      cat("Lc prior   =",Lc.pr,", SD =",Lc.sd.pr,"(cm)",ifelse(is.na(LcUser)==TRUE,"","(user-defined)"),
      ", alpha prior=",r.alpha.pr,", SD =",0.1*r.alpha.pr,"\n\n")
    }

    cat("General reference points [median across years]: \n")
    cat("Linf               =",format(Linf.med,digits=3),
        paste("(",format(Linf.lcl,digits=3),"-",format(Linf.ucl,digits=3),
                            ifelse(mmUser==F,") cm",") mm"), sep=""), "\n")  
    cat("Lopt               =",format(Lopt.med,digits=2),
        paste(ifelse(mmUser==F,"cm,","mm,"),"Lopt/Linf ="),format(r.Lopt.med,digits=2),"\n")
    cat("Lc_opt             =",format(Lc_opt.med,digits=2),
        paste(ifelse(mmUser==F,"cm,","mm,"),"Lc_opt/Linf ="),
        format(Lc_opt.med/Linf.med,digits=2),"\n")
    cat("M/K                =",format(MK.med,digits=3),
        paste("(",format(MK.lcl,digits=3),"-",format(MK.ucl,digits=3),
                            ")",sep=""),"\n")
    cat("F/K                =",format(FK.med,digits=3),
        paste("(",format(FK.lcl,digits=3),"-",format(FK.ucl,digits=3),
                                         ")",sep=""),"\n")
    cat("Z/K                =",format(ZK.med,digits=3),
        paste("(",format(ZK.lcl,digits=3),"-",format(ZK.ucl,digits=3),
                                         ")",sep=""),"\n")
    cat("F/M                =",format(FM.med,digits=3),
        paste("(",format(FM.lcl,digits=3),"-",format(FM.ucl,digits=3),
                                         ")",sep=""),"\n")
    cat(ifelse(GausSel==F,"B/B0 F=M Lc=Lc_opt =","B/B0 F=M Lmean=Lopt="),
        format(BFM1B0,digits=3),"\n")
    cat("B/B0               =",format(BB0.med,digits=3),
        paste("(",format(BB0.lcl,digits=3),"-",format(BB0.ucl,digits=3),
                                         ")",sep=""),"\n")
    cat(ifelse(GausSel==F,"Y/R' F=M Lc=Lc_opt =","Y/R' F=M Lmean=Lopt="),
        format(YRFM1,digits=3),"\n")
    cat("Y/R'               =",format(YR.med,digits=3),
        paste("(",format(YR.lcl,digits=3),"-",format(YR.ucl,digits=3),
                                          ")",sep=""),"(linearly reduced if B/B0 < 0.25)\n\n")
    cat("Estimates for last year",EndYear,":\n")
    last            <- which(Ldat$Year==EndYear)
    if(GausSel==F){
        cat("Lc         =",format(Ldat$Lc[last],digits=3),
            paste("(",format(Ldat$Lc.lcl[last],digits=3),
                        "-",format(Ldat$Lc.ucl[last],digits=3),ifelse(mmUser==F,") cm, Lc/Linf = ",") mm, Lc/Linf = "),
                  format(Ldat$Lc[last]/Ldat$Linf[last],digits=2)," (",
                  format(Ldat$Lc.lcl[last]/Ldat$Linf[last],digits=3),"-",
                        format(Ldat$Lc.ucl[last]/Ldat$Linf[last],digits=3),")",sep=""),"\n")
        cat("alpha      =",format(Ldat$r.alpha[last],digits=3),"(",
            format(Ldat$r.alpha.lcl[last],digits=3),"-",
            format(Ldat$r.alpha.ucl[last],digits=3),") \n")
     cat("Lmean/Lopt =",format(Ldat$Lmean[last]/(Ldat$r.Lopt[last]*Ldat$Linf[last]),digits=2),
         ", Lc/Lc_opt =",format(Ldat$Lc[last]/Lc_opt.med,digits=2),
         ", L95th =", format(Ldat$L95[last],digits=3),ifelse(mmUser==F,"cm","mm"),
         ", L95th/Linf =",format(Ldat$L95[last]/Ldat$Linf[last],digits=2),
         ", \nLm50 =", format(res$Lm50,digits=3),ifelse(mmUser==F,"cm","mm"),
         ", Mature =",format(Ldat$perc.mat[last]*100,digits=2),"%\n")
     } else if(GausSel==T){
         cat("GLmean/Linf=",format(Ldat$GLmean[last]/Ldat$Linf[last],digits=2),",SD/Linf =",
             format(Ldat$SD[last]/Ldat$Linf[last],digits=3),"\n")
         cat("GLmean     =",format(Ldat$GLmean[last],digits=3),",SD =",
             format(Ldat$SD[last],digits=3),"\n")
     }
    cat("F/K        =",format(Ldat$FK[last],digits=2),"(",
        format(Ldat$FK.lcl[last],digits=3),"-",
        format(Ldat$FK.ucl[last],digits=3),")\n")
    cat("F/M        =",format(Ldat$FK[last]/Ldat$MK[last],digits=2),"(",
        format(Ldat$FM.lcl[last],digits=3),"-",
        format(Ldat$FM.ucl[last],digits=3),")\n")
    cat("Z/K        =",format(Ldat$ZK[last],digits=3),"(",
        format(Ldat$ZK.lcl[last],digits=3),"-",
        format(Ldat$ZK.ucl[last],digits=3),")\n")
    cat("Y/R'       =",format(Ldat$YR[last],digits=2),"(",
        format(Ldat$YR.lcl[last],digits=3),"-",
        format(Ldat$YR.ucl[last],digits=3),") (linearly reduced if B/B0 < 0.25)\n")
    cat("B/B0       =",format(Ldat$BB0[last],digits=2),"(",
        format(Ldat$BB0.lcl[last],digits=3),"-",format(Ldat$BB0.ucl[last],digits=3),")\n")
    cat("B/Bmsy     =",format(Ldat$BB0[last]/BFM1B0,digits=2),"(",
        format(Ldat$BB0.lcl[last]/BFM1B0,digits=3),"-",
        format(Ldat$BB0.ucl[last]/BFM1B0,digits=3),")\n")
    if(Comment != "" & !is.na(Comment)) cat("Comment:",Comment,"\n")

    # point out questionable or impossible results
    # negative rates
    if(Ldat$MK[last] < 0 | Ldat$FK[i] < 0) cat("Data unsuitable for LF analysis, negative mortality rates are impossible\n")
    # Biomass larger than unexploited
    if(Ldat$BB0[last] >1.1) cat("Data unsuitable for LF analysis, biomass exceeds carrying capacity\n")

    flush.console()





    ## return results
    
    priors <- t(data.frame(Linf = c(Linf.pr, Linf.sd.pr, LinfUser),
                                    ZK = c(ZK.nls, ZK.nls.sd, NA),
                                    MK = c(MK.pr, MK.sd.pr, MKUser),
                                    FK = c(FK.pr, NA, NA),
                                    Lc = c(Lc.pr, Lc.sd.pr, LcUser),
                                    alpha = c(r.alpha.pr, 0.1*r.alpha.pr, NA)))
    colnames(priors) <- c("prior", "SD", "user")

    
    medianRefLev <- t(data.frame(Linf = c(Linf.med, Linf.lcl, Linf.ucl,
                                                 mmUser), ## ifelse(mmUser==F,") cm",") mm")
                                        Lopt = c(Lopt.med, NA, NA,
                                                 mmUser), ## ifelse(mmUser==F,") cm",") mm")
                                        Lc_opt = c(Lc_opt.med, NA, NA,
                                                   mmUser), ## ifelse(mmUser==F,") cm",") mm")
                                        ZK = c(ZK.med, ZK.lcl, ZK.ucl, NA),
                                        MK = c(MK.med, MK.lcl, MK.ucl, NA),
                                        FK = c(FK.med, FK.lcl, FK.ucl, NA),
                                        FM = c(FM.med, FM.lcl, FM.ucl, NA),
                                        BB0 = c(BB0.med, BB0.lcl, BB0.ucl, NA),
                                        YR = c(YR.med, YR.lcl, YR.ucl, NA) ## "(linearly reduced if B/B0 < 0.25)\n\n"
                                        ))
    ## cat(ifelse(GausSel==F,"B/B0 F=M Lc=Lc_opt =","B/B0 F=M Lmean=Lopt="),BFM1B0,"\n")
    ## cat(ifelse(GausSel==F,"Y/R' F=M Lc=Lc_opt =","Y/R' F=M Lmean=Lopt="),YRFM1,"\n")

    colnames(medianRefLev) <- c("median", "lo", "up", "user")


    last <- which(Ldat$Year == EndYear)
    lastRefLev <- t(data.frame(ZK = c(Ldat$ZK[last], Ldat$ZK.lcl[last],
                                             Ldat$ZK.ucl[last], NA),
                                      MK = c(Ldat$MK[last], Ldat$MK.lcl[last],
                                             Ldat$MK.ucl[last], NA),
                                      FK = c(Ldat$FK[last], Ldat$FK.lcl[last],
                                             Ldat$FK.ucl[last], NA),
                                      FM = c(Ldat$FM[last], Ldat$FM.lcl[last],
                                             Ldat$FM.ucl[last], NA),
                                      BB0 = c(Ldat$BB0[last], Ldat$BB0.lcl[last],
                                             Ldat$BB0.ucl[last], NA),
                                      YR = c(Ldat$YR[last], Ldat$YR.lcl[last],
                                             Ldat$YR.ucl[last], NA),
                                      BBmsy = c(Ldat$BB0[last]/BFM1B0, Ldat$BB0.lcl[last]/BFM1B0,
                                             Ldat$BB0.ucl[last]/BFM1B0, NA)
                                      ))
    colnames(lastRefLev) <- c("median", "lo", "up", "mm")

    if(!GausSel){
        lastRefLev <- rbind(lastRefLev,
                            t(data.frame(Lc = c(Ldat$Lc[last], Ldat$Lc.lcl[last],
                                             Ldat$Lc.ucl[last], mmUser),
                                         LcLinf = c(Ldat$Lc[last]/Ldat$Linf[last],
                                                    Ldat$Lc.lcl[last]/Ldat$Linf[last],
                                                    Ldat$Lc.ucl[last]/Ldat$Linf[last], NA),
                                         alpha = c(Ldat$r.alpha[last], Ldat$r.alpha.lcl[last],
                                                   Ldat$r.alpha.ucl[last], NA),
                                         LmeanLopt = c(Ldat$Lmean[last]/(Ldat$r.Lopt[last]*Ldat$Linf[last]),
                                                        NA, NA, NA),
                                         LcLc_opt = c(Ldat$Lc[last]/Lc_opt.med,
                                                       NA, NA, NA),
                                         L95th = c(Ldat$L95[last],
                                                   NA, NA, mmUser),
                                         L95thLinf = c(Ldat$L95[last]/Ldat$Linf[last],
                                                        NA, NA, NA),
                                         Lm50 = c(res$Lm50, NA, NA, mmUser),
                                         Mature = c(Ldat$perc.mat[last] * 100, NA,NA,NA))))
    }else if(GausSel){
        lastRefLev <- rbind(lastRefLev,
                            t(data.frame(GLmean = c(Ldat$GLmean[last], Ldat$SD[last],
                                             NA, NA),
                                         GLmeanLinf = c(Ldat$GLmean[last]/Ldat$Linf[last],
                                                    Ldat$SD[last]/Ldat$Linf[last],
                                                    NA, NA))))
    }
    
    ret <- c(res, list(GausSel = GausSel,
                      priors = priors,
                      refLev = Ldat,
                      medianRefLev = medianRefLev,
                      lastRefLev = lastRefLev,
                      LFall = LF.all
                  ))

    return(ret)
}







#' @title B & H equations for LBB
#'
#' @description Exploited B/B0 ratio from B&H equations for variable F
#'
#' @param AllLength Vector with all lengths
#' @param Linf Asymptotic length
#' @param MK MK ratio
#' @param FK FK ratio
#' @param GausSel Logical; indicating the selectivity pattern. If FALSE (default) trawl-like,
#'    if TRUE gaussian selectivity is assumed.
#' @param selpar1 Selectivity parameter 1
#' @param selpar2 Selectivity parameter 2
#'
#' @details Assuming that observed lengths are the lower bounds of length
#'     classes get lowest exploited (>= 0.01 F) length class and class width.
#' 
#' @return A list with B/B0 and Y/R
#'
#' @references
#' R. Froese, H. Winker, G. Coro, N. Demirel, A.C. Tsikliras, D. Dimarchopoulou,
#' G. Scarcella, W.N. Probst, M. Dureuil, and D. Pauly (2018) A new approach
#' for estimating stock status from length frequency data. ICES Journal of Marine Science. DOI: 10.1093/icesjms/fsy078
#'
BH <- function(AllLength,Linf,MK,FK,GausSel,selpar1,selpar2) {

    if(GausSel==FALSE){
        r.Lc     <- selpar1
        r.alpha  <- selpar2 
        Lx       <- AllLength[AllLength >= Linf*(r.Lc-4.59/r.alpha)][1]
    }else if(GausSel==TRUE){
        r.GLmean <- selpar1
        r.SD     <- selpar2
        Lx       <- AllLength[AllLength >= Linf*(r.GLmean-3*r.SD)][1]
    }
    ##    class.width  <- AllLength[2]-AllLength[1]   ## only works if class widths are constant
    ## updated:
    class.width  <- median(diff(sort(unique(AllLength))))    
    FM <- FK/MK

    # Linf=120;Lx=22.5;r.Lc=0.2917;r.alpha=60;MK=1.5385;FK=0.7692;FM=0.5;ZK=2.3077 
    # uncomment above row for comparison of Y'R= 0.0332, B/B0=0.467 with CodLightSim
    r            <- vector() # auxilliary reduction factor
    G            <- vector() # product of reduction factors
    SL.bh        <- vector() # selection at length
    YR1.2        <- vector() # relative yield per recruit per length class
    CPUER1.2     <- vector() # relative CPUE per recruit per length class
    B1.2         <- vector() # relative unexploited biomass per recruit by length class
    L.bh         <- seq(from=Lx, to=Linf, by=class.width) # lengths to be considered
    r.L.bh       <-  L.bh / Linf # standardized lengths

    # calculate selection, Y'/R and CPUE'/R for every length class
    for(o in 1 : length(r.L.bh)) { 
      if(GausSel==F) {
        if(o<length(r.L.bh)) { SL.bh[o] <- mean(c(1/(1+exp(-r.alpha*(r.L.bh[o]-r.Lc))), # mean selection in length class
                                                  1/(1+exp(-r.alpha*(r.L.bh[o+1]-r.Lc)))))
        } else SL.bh[o] <- 1/(1+exp(-r.alpha*(r.L.bh[o]-r.Lc)))
      } else if(GausSel==T) { # gill net selection 
        if(o<length(r.L.bh)) { SL.bh[o] <- mean(c(exp(-((r.L.bh[o]-r.GLmean)^2/(2*r.SD^2))), # mean selection in length class
                                                  exp(-((r.L.bh[o+1]-r.GLmean)^2/(2*r.SD^2)))))
        } else SL.bh[o] <- exp(-((r.L.bh[o]-r.GLmean)^2/(2*r.SD^2)))
    } # end of calculation of selectivity loop

      if(o<length(r.L.bh)) {
        r[o]       <- (1-r.L.bh[o+1])^(FK*SL.bh[o])/(1-r.L.bh[o])^(FK*SL.bh[o]) 
        G[o]       <- prod(r[1:o]) }
      if(o==1) {
        YR1.2[o] <-(FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/
                    (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o])^2/(1+2/(MK+FK*SL.bh[o]))-
                     (1-r.L.bh[o])^3/(1+3/(MK+FK*SL.bh[o])))) -
                      (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o+1])^MK*(1-3*(1-r.L.bh[o+1])/(1+1/
                       (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o+1])^2/(1+2/(MK+FK*SL.bh[o]))-
                        (1-r.L.bh[o+1])^3/(1+3/(MK+FK*SL.bh[o]))))*G[o] 
      } else if(o==length(r.L.bh)) {
        YR1.2[o] <- (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/
                     (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o])^2/(1+2/(MK+FK*SL.bh[o]))-
                      (1-r.L.bh[o])^3/(1+3/(MK+FK*SL.bh[o])))) * G[o-1] 
      } else {
        YR1.2[o] <- (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/
                     (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o])^2/(1+2/(MK+FK*SL.bh[o]))-
                      (1-r.L.bh[o])^3/(1+3/(MK+FK*SL.bh[o])))) * G[o-1] -
                       (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o+1])^MK*(1-3*(1-r.L.bh[o+1])/(1+1/
                        (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o+1])^2/(1+2/(MK+FK*SL.bh[o]))-
                         (1-r.L.bh[o+1])^3/(1+3/(MK+FK*SL.bh[o]))))*G[o]              
      } # end of loop to calculate yield per length class

      CPUER1.2[o] <- YR1.2[o] / FM # CPUE/R = Y/R divided by F/M

      if(o<length(r.L.bh)) {
        B1.2[o] <- ((1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/MK)+3*(1-r.L.bh[o])^2/
                                               (1+2/MK)-(1-r.L.bh[o])^3/(1+3/MK)) -
          (1-r.L.bh[o+1])^MK*(1-3*(1-r.L.bh[o+1])/(1+1/MK)+3*(1-r.L.bh[o+1])^2/
                                        (1+2/MK)-(1-r.L.bh[o+1])^3/(1+3/MK)))*SL.bh[o]
      } else {
        B1.2[o] <- ((1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/MK)+3*(1-r.L.bh[o])^2/
                                               (1+2/MK)-(1-r.L.bh[o])^3/(1+3/MK)))*SL.bh[o]
      }
    } # end of B&H loop through length classes
    BB0   <- sum(CPUER1.2)/sum(B1.2)
    YR    <- sum(YR1.2)
    if(BB0 < 0.25) YR <- YR * BB0 / 0.25 # reduce YR if recruitment and thus productivity is reduced
    return(list(BB0,YR))
}






#' @title Aggregation function for LBB
#'
#' @description Function to aggregate data by year
#'
#' @param dat Data.frame with columns: Year, Length (in cm), and CatchNo
#' 
#' @return A data.frame with columns Length and Freq
#'
#' @references
#' R. Froese, H. Winker, G. Coro, N. Demirel, A.C. Tsikliras, D. Dimarchopoulou,
#' G. Scarcella, W.N. Probst, M. Dureuil, and D. Pauly (2018) A new approach
#' for estimating stock status from length frequency data. ICES Journal of Marine Science. DOI: 10.1093/icesjms/fsy078
#'
AG  <- function(dat){
    # aggregate normalized annual LFs by weighing with square root of sample size
    # get sum of frequencies per year
    sum.Ny  <- aggregate(Freq~Year,dat,sum)$Freq  
    # get the sqrt of the sum of frequencies for every year
    sqrt.Ny <- sqrt(sum.Ny) 
    # get highest frequency in each year
    max.Ny <- aggregate(Freq~Year,dat,max)$Freq
    # get Number of Length bins in each year
    binsN <- aggregate(Freq~Year,dat,length)$Freq    
    # create vectors for sqrt.Ni and sum.Ni to weigh LF data
    sqrt.Ni = rep(sqrt.Ny,binsN)
    sum.Ni = rep(sum.Ny,binsN)
    #Do weighing
    # Divide all years by sum.Ni and multiply by sqrt.Ni
    LF.w = dat$Freq/sum.Ni*sqrt.Ni  
    # Aggregate
    LF = aggregate(LF.w, by=list(dat$Length),FUN=sum)
    # Add correct column names
    colnames(LF) <- c("Length","Freq")         
    return(LF)
}






#' @title Plotting LBB fit (single year)
#'
#' @description Function to plot LBB fit for a single year
#'
#' @param r.L.y Length classes relative to asymptotic length.
#' @param r.Freq.y Relative frequencies.
#' @param r.Lopt Optimum length relative to asymptotic length.
#' @param SL1 Selectivity parameter 1
#' @param SL2 Selectivity parameter 2
#' @param MK MK ratio
#' @param FK FK ratio
#' @param Linf Asymptotic length
#' @param Year Assessment year
#' @param GausSel Logical; indicating the selectivity pattern. If FALSE (default) trawl-like,
#'    if TRUE gaussian selectivity is assumed.
#'
#' @details Expects lengths relative to Linf (L/Linf).
#' 
#' @references
#' R. Froese, H. Winker, G. Coro, N. Demirel, A.C. Tsikliras, D. Dimarchopoulou,
#' G. Scarcella, W.N. Probst, M. Dureuil, and D. Pauly (2018) A new approach
#' for estimating stock status from length frequency data. ICES Journal of Marine Science. DOI: 10.1093/icesjms/fsy078
#'
plotLBB.year <- function(r.L.y,
                         r.Freq.y,
                         r.Lopt,
                         SL1, SL2,
                         MK, FK, Linf,
                         Year,
                         GausSel = FALSE){
    SL       <- vector()
    SL[1]  <- 0
    xN       <- vector()
    xN[1]  <- 1  
    for(p in 2:length(r.L.y)) {
      if(GausSel){
        SL[p]  <- exp(-((r.L.y[p]-SL1/Linf)^2/(2*(SL2/Linf)^2))) # selection at length r.L[i]
      } else {
          SL[p]  <- 1/(1+exp(-SL2*(r.L.y[p]-SL1/Linf)))
      }
      xN[p]  <- xN[p-1]*exp((MK+FK*SL[p])*(log(1-r.L.y[p])-log(1-r.L.y[p-1])))
    } # end of loop for length classes
    # determine sum of squared residuals, store parameters if less then min of previous
    Freq.pred      <- xN*SL
    sum.Freq.pred  <- sum(Freq.pred)
    r.Freq.pred    <- Freq.pred / sum.Freq.pred          
    plot(x=r.L.y, y= r.Freq.pred,
         xlab="Length / Linf",ylab="relative Frequency",
         xlim=c(0,1),ylim = c(0,1.2*max(r.Freq.y)),
         col="red", type="l", bty="l",main=Year,las=1)
    points(x=r.L.y,y=r.Freq.y, cex=0.5)
    lines(x=c(1,1), y=c(0,1.07*max(r.Freq.y,na.rm=T)),col="darkgreen")
    text(x=1,y=1.15*max(r.Freq.y,na.rm=T),"Linf",col="darkgreen") 
    lines(x=c(r.Lopt,r.Lopt), y=c(0,1.07*max(r.Freq.y,na.rm=T)),col="darkgreen")
    text(x=r.Lopt,y=1.15*max(r.Freq.y,na.rm=T),"Lopt",col="darkgreen") 
    text(x=0.15,y=0.8*max(r.Freq.y,na.rm=T),paste("Linf=",format(Linf,digits=3),sep=""))
    text(x=0.15,y=0.6*max(r.Freq.y,na.rm=T),paste("Z/K=",format(MK+FK,digits=3),sep=""))
} 




#' @title Plotting LBB data
#'
#' @description Function to plot length-frequency data in the LBB format
#'
#' @param lfq A list of the class "lfq" consisting of following parameters:
#' \itemize{
#'   \item \strong{species} species name,
#'   \item \strong{stock} stock ID or name,
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row)
#'      and sampling date (column),
#'   \item \strong{comments} comments;
#' }
#' @param mfrow A vector of the form ‘c(nr, nc)’.  Subsequent figures will be drawn in an
#'    ‘nr’-by-‘nc’ array on the device by _rows_ (‘mfrow’). If NA (default), a panel with
#'    3 columns and several rows (dependent on number of years) is used.
#' 
#' @details expects lengths relative to Linf (L/Linf)
#' 
#' @references
#' R. Froese, H. Winker, G. Coro, N. Demirel, A.C. Tsikliras, D. Dimarchopoulou,
#' G. Scarcella, W.N. Probst, M. Dureuil, and D. Pauly (2018) A new approach
#' for estimating stock status from length frequency data. ICES Journal of Marine Science. DOI: 10.1093/icesjms/fsy078
#'
#' @export
#' 
plotLBB.data <- function(lfq, mfrow=NA){

    res <- lfq

    years <- as.numeric(format(res$dates, "%Y"))
    
    ##--------------------------------------------------------------------------------------    
    ## Put data into vectors
    ##--------------------------------------------------------------------------------------
    if(class(lfq$catch) == "matrix"){
        nrowC <- nrow(lfq$catch)
        ncolC <- ncol(lfq$catch)
    }else if(class(lfq$catch) == "numeric"){
        nrowC <- length(lfq$catch)
        ncolC <- 1
    }else{
        stop("lfq$catch has to be either a matrix or a numeric vector!")
    }
    
    
    StartYear  <- min(years)
    EndYear    <- max(years)
    AllYear    <- rep(years, each = nrowC)
    AllLength  <- rep(res$midLengths, ncolC)
    AllFreq    <- as.numeric(res$catch)
    Years      <- sort(unique(AllYear))
    nYears     <- length(years)

    for(z in 1:ceiling(nYears/6)) {
        if(is.na(mfrow)){
            ncols <- ifelse(length(Years)<=4,2,3)
            if(length(Years) == 1) ncols <- 1            
            opar <- par(mfrow=c(ceiling(length(Years)/3),ncols))
        }else{
            opar <- par(mfrow=mfrow)
        }        
        for(v in 1 : 6) {
            w <- v+(z-1)*6
            if(w > nYears) break()
            df.p <- data.frame(AllYear[AllYear==Years[w]&AllFreq>0],
                               AllLength[AllYear==Years[w]&AllFreq>0],AllFreq[AllYear==Years[w]&AllFreq>0])
            names(df.p) <- c("Year","Length","Freq")
            LF.p        <- AG(dat=df.p) # function to aggregate data in case bins are not unique
            plot(x=LF.p$Length,y=LF.p$Freq,xlab="",ylab="Freq",bty="l",main=Years[w],cex=0.5)
         }
    }
    on.exit(par(opar))
}





#' @title Plotting LBB fit
#'
#' @description Function to plot data and fit of a LBB assessment
#'
#' @param lfq A list of the class "lfq" consisting of following parameters:
#' \itemize{
#'   \item \strong{species} species name,
#'   \item \strong{stock} stock ID or name,
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row)
#'      and sampling date (column),
#'   \item \strong{comments} comments;
#' }
#' @param LcUser Optional, user defined selectivity parameter. If NA (default) L10 and L90 are
#'    used to estimate a proxy for Lc.
#' @param LstartUser Optional, user defined length at which selectivity is 0.95. If NA (default)
#'    Lstart (L95) is estimated by the model.
#' @param MKUser Optional; user defined MK ratio. If NA (default) MK ratio is set to 1.5.
#' @param mmUser Logical; indicating the unit of length measurements, where TRUE
#'    indicates that lengths are in mm and FALSE (default) indicate that lengths are in cm.
#' @param GausSel Logical; indicating the selectivity pattern. If FALSE (default) trawl-like,
#'    if TRUE gaussian selectivity is assumed.
#' 
#' @details Plot aggregated histogram with fit to fully selected part.
#' 
#' @references
#' R. Froese, H. Winker, G. Coro, N. Demirel, A.C. Tsikliras, D. Dimarchopoulou,
#' G. Scarcella, W.N. Probst, M. Dureuil, and D. Pauly (2018) A new approach
#' for estimating stock status from length frequency data. ICES Journal of Marine Science. DOI: 10.1093/icesjms/fsy078
#'
#' @export
#' 
plotLBB <- function(lfq, GausSel = FALSE, LcUser = NA, mmUser = FALSE, MKUser=NA, LstartUser=NA){
    LF.all = lfq$LFall
    # If no Linf is provided by the user (preferred), determine Linf from fully selected LF:
    # Freq=Nstart*exp(ZK*(log(1-L/Linf)-log(1-Lstart/Linf)))
    # Nstart is canceled out when dividing both sides by their sums
    # ---------------------------------------------------------
    # determine start values of selection ogive to find first fully selected length class Lstart
    L10 <- LF.all$Length[which(LF.all$Freq>0.1)[1]] # use length at 10% of peak frequency as proxy for L10
    L90 <- LF.all$Length[which(LF.all$Freq>0.9)[1]] # use length at 90% of peak frequency as proxy for L90
    Lc.st <- ifelse(is.na(LcUser)==TRUE,(L10 + L90)/2,LcUser)  # use mean of L10 and L90 as proxy for Lc, else user input
    Linf.st     <- max(LF.all$Length) # use Lmax as proxy for Linf
    Lmean.st    <- sum(LF.all$Length[LF.all$Length>=Lc.st]*LF.all$Freq[LF.all$Length>=Lc.st])/
                         sum(LF.all$Freq[LF.all$Length>=Lc.st])
    MK.st       <- ifelse(is.na(MKUser)==TRUE, 1.5,MKUser) # default 1.5
    ZK.st       <- (Linf.st-Lmean.st)/(Lmean.st-Lc.st)       # the Holt equation
    FK.st       <- ifelse((ZK.st-MK.st)>0,ZK.st-MK.st,0.3)   # prevent M/K being larger than Z/K
    alpha.st <- -log(1/LF.all$Freq[which(LF.all$Freq>0.1)[1]])/(L10-Lc.st) # use rearranged logistic curve to estimate slope alpha    

    # get vectors with fully selected length classes for Linf estimation
    if(!is.na(LstartUser)) {
        Lstart <- LstartUser
    } else {
        Lstart     <- (alpha.st*Lc.st-log(1/0.95-1))/alpha.st   # Length where selection probability is 0.95  
        # test if there are enough (>=4) length classes for estimation of aggregated Linf and ZK 
        Lstart.i   <- which(LF.all>=Lstart)[1]
        Lmax.i     <- length(LF.all$Length)
        peak.i     <- which.max(LF.all$Freq)
        if(Lstart.i<(peak.i+1)) Lstart <- LF.all$Length[peak.i+1] # make sure fully selected length starts after peak 
        if((Lmax.i-Lstart.i)<4) Lstart <- LF.all$Length[Lstart.i-1] # make sure enough length classes are available
    }
    # do not include Lmax to allow Linf < Lmax and to avoid error in nls when Linf-L becomes negative
    L.L         <- LF.all$Length[LF.all$Length >= Lstart  & LF.all$Length < Linf.st]
    L.Freq      <- LF.all$Freq[LF.all$Length>=L.L[1]& LF.all$Length < Linf.st]
    ##
    ZK.nls <- lfq$priors[rownames(lfq$priors) == "ZK",1]
    Linf.nls <- lfq$priors[rownames(lfq$priors) == "Linf",1]
    ##
    plot(x=LF.all$Length,y=LF.all$Freq, bty="l",xlim=c(0,max(max(LF.all$Length),Linf.nls)),
         ylim=c(0,1.1*max(LF.all$Freq)),
         main=paste(lfq$stock,", aggregated LF"),
         xlab=ifelse(mmUser==F,"Length (cm)","Length (mm)"),ylab="Frequency")
    Lstart.i    <- which(LF.all$Length>=Lstart)[1]
    Lstart.Freq <- mean(c(LF.all$Freq[(Lstart.i-1):(Lstart.i+1)]))
    if(!GausSel){
        lines(x=L.L,y=Lstart.Freq*exp(ZK.nls*(log(1-L.L/Linf.nls)-log(1-L.L[1]/Linf.nls))),
              col="blue", lwd=3)
        lines(x=c(Lc.st,Lc.st), y=c(0,1), col="darkgreen")
        text(x=Lc.st,y=1, "Lc", col="darkgreen", adj=c(0.5,-0.5)) 
    }
    lines(x=c(Linf.nls,Linf.nls), y=c(0,1), col="darkgreen")
    text(x=Linf.nls,y=1, "Linf", col="darkgreen", adj=c(0.5,-0.5))
    text(x=0.05*Linf.nls,y=1,"Priors:")
    text(x=0.1*Linf.nls,y=0.9,paste("Linf=",format(Linf.nls,digits=3),sep=""))
    if(!GausSel) text(x=0.1*Linf.nls,y=0.8,paste("Z/K=",format(ZK.nls,digits=2),sep=""))
    text(x=0.1*Linf.nls,y=0.7,paste("Lc=",format(Lc.st,digits=3),sep=""))
}




#' @title Plotting LBB results over time
#'
#' @description Function to plot main results of a LBB assessment as time series graphs
#'
#' @param lfq A list of the class "lfq" consisting of following parameters:
#' \itemize{
#'   \item \strong{species} species name,
#'   \item \strong{stock} stock ID or name,
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row)
#'      and sampling date (column),
#'   \item \strong{comments} comments;
#' }
#' @param mmUser Logical; indicating the unit of length measurements, where TRUE
#'    indicates that lengths are in mm and FALSE (default) indicate that lengths are in cm.
#' @param GausSel Logical; indicating the selectivity pattern. If FALSE (default) trawl-like,
#'    if TRUE gaussian selectivity is assumed.
#'
#' @details Expects lengths relative to Linf (L/Linf)
#' 
#' @references
#' R. Froese, H. Winker, G. Coro, N. Demirel, A.C. Tsikliras, D. Dimarchopoulou,
#' G. Scarcella, W.N. Probst, M. Dureuil, and D. Pauly (2018) A new approach
#' for estimating stock status from length frequency data. ICES Journal of Marine Science. DOI: 10.1093/icesjms/fsy078
#'
#' @export
#' 
plotLBB.ts <- function(lfq, mmUser = FALSE, GausSel = FALSE){

    res <- lfq
    Ldat <- res$refLev

    years <- as.numeric(format(res$dates, "%Y"))
    
    ##--------------------------------------------------------------------------------------    
    ## Put data into vectors
    ##--------------------------------------------------------------------------------------
    if(class(lfq$catch) == "matrix"){
        nrowC <- nrow(lfq$catch)
        ncolC <- ncol(lfq$catch)
    }else if(class(lfq$catch) == "numeric"){
        nrowC <- length(lfq$catch)
        ncolC <- 1
    }else{
        stop("lfq$catch has to be either a matrix or a numeric vector!")
    }
    
    
    StartYear  <- min(years)
    EndYear    <- max(years)
    AllYear    <- rep(years, each = nrowC)
    AllLength  <- rep(res$midLengths, ncolC)
    AllFreq    <- as.numeric(res$catch)
    Years      <- sort(unique(AllYear))
    nYears     <- length(years)


    ##--------------------------------------------------------------------------------------    
    ## get some reference points as median of time series
    ##--------------------------------------------------------------------------------------        
    Linf.med     <- median(Ldat$Linf)
    Linf.lcl     <- median(Ldat$Linf.lcl)
    Linf.ucl     <- median(Ldat$Linf.ucl)
    if(GausSel==F) {
      Lc.med       <- median(Ldat$Lc)
      r.alpha.med  <- median(Ldat$r.alpha)
    } else {
      GLmean.med   <- median(Ldat$GLmean)
      SD.med       <- median(Ldat$SD)
    }
    MK.med       <- median(Ldat$MK)
    MK.lcl       <- median(Ldat$MK.lcl)
    MK.ucl       <- median(Ldat$MK.ucl)
    FK.med       <- median(Ldat$FK)
    FK.lcl       <- median(Ldat$FK.lcl)
    FK.ucl       <- median(Ldat$FK.ucl)
    FM.med       <- median(Ldat$FM)
    FM.lcl       <- median(Ldat$FM.lcl)
    FM.ucl       <- median(Ldat$FM.ucl)
    ZK.med       <- median(Ldat$ZK)
    ZK.lcl       <- median(Ldat$ZK.lcl)
    ZK.ucl       <- median(Ldat$ZK.ucl)
    r.Lopt.med   <- median(Ldat$r.Lopt)
    Lopt.med     <- r.Lopt.med*Linf.med
    Lc_opt.med   <- Linf.med*(2+3*FM.med)/((1+FM.med)*(3+MK.med)) 
    BB0.med      <- median(Ldat$BB0)
    BB0.lcl      <- median(Ldat$BB0.lcl)
    BB0.ucl      <- median(Ldat$BB0.ucl)
    YR.med       <- median(Ldat$YR)
    YR.lcl       <- median(Ldat$YR.lcl)
    YR.ucl       <- median(Ldat$YR.ucl)

    BFM1B0.list  <- BH(AllLength=unique(AllLength),Linf=Linf.med,MK=MK.med,FK=MK.med,GausSel=GausSel,
                         selpar1=ifelse(GausSel==T,r.Lopt.med,5/(2*(3+MK.med))),
                         selpar2=ifelse(GausSel==T,SD.med/Linf.med,r.alpha.med))

    BFM1B0       <- as.numeric(BFM1B0.list[1])
    YRFM1        <- as.numeric(BFM1B0.list[2])


    opar <- par(mfrow=c(1,3))
    #----------------------------------------------
    #  Plot time series of Lc and Lmean
    #----------------------------------------------
    if(nYears > 1) {

        if(!GausSel){
            plot(x=Ldat$Year,y=Ldat$Lmean, bty="l",type="l", 
                 xlim=c(Ldat$Year[1],Ldat$Year[nYears]),
                 xaxt="n",
                 ylim=c(0,max(c(1.1*Lopt.med,max(Ldat$Lmean,na.rm=T),max(Ldat$Lc.ucl),na.rm=T))),lwd=2,
                 xlab="Year",ylab = paste("Length",ifelse(mmUser==F,"(cm)","(mm)")),main="Lmean vs Lopt & Lc vs Lc_opt")
            axis(1,at=Ldat$Year)
            lines(x=Ldat$Year,y=Ldat$Lc,lwd=1,lty="dashed")
            #lines(x=Ldat$Year,y=Ldat$Lc.lcl,lty="dotted")
            #lines(x=Ldat$Year,y=Ldat$Lc.ucl,lty="dotted")
            lines(x=Ldat$Year,y=rep(Lc_opt.med,nYears),col="darkgreen", lty="dashed") # line for Lc_opt
            text(x=Ldat$Year[nYears],y=Lc_opt.med,"Lc_opt", adj=c(1,-0.5), col="darkgreen")
            lines(x=Ldat$Year,y=rep(Lopt.med,nYears),col="darkgreen") # line for Lopt
            text(x=Ldat$Year[nYears],y=Lopt.med,"Lopt", adj=c(1,-0.5), col="darkgreen")
        }

        #----------------------------------------------
        #  Plot time series of GLmean relative to Lopt 
        #----------------------------------------------
          if(GausSel==T){
            plot(x=Ldat$Year,y=Ldat$GLmean, bty="l",type="l", 
                 xlim=c(Ldat$Year[1],Ldat$Year[nYears]),
                 xaxt="n",
                 ylim=c(0,max(1.1*median(Ldat$r.Lopt)*Linf.med,max(Ldat$GLmean),na.rm=T)),lwd=2,
                 xlab="Year",ylab = "Lenght (cm)",main="Lmean vs Lopt")
            axis(1,at=Ldat$Year)
            lines(x=Ldat$Year,y=(Ldat$GLmean.lcl),lty="dotted")
            lines(x=Ldat$Year,y=(Ldat$GLmean.ucl),lty="dotted")
            lines(x=Ldat$Year,y=rep(Lopt.med,nYears),col="darkgreen") # line for Lopt
            text(x=Ldat$Year[nYears],y=Lopt.med,"Lopt", adj=c(1,-0.5), col="darkgreen")
        }    
        #---------------------------------------------
        # Plot time series of F/M 
        #---------------------------------------------
        plot(x=Ldat$Year,y=Ldat$FM,
               ylim=c(0,max(max(Ldat$FM.ucl),1.05)), 
               bty="l",type = "l", lwd=1.5, xaxt="n",
               main="previous F/M",xlab="Year",ylab="F/M")
        axis(1,at=Ldat$Year)
        lines(x=Ldat$Year,y=Ldat$FM.lcl,lty="dotted")
        lines(x=Ldat$Year,y=Ldat$FM.ucl,lty="dotted")
        abline(h=1.0,col="darkgreen") 
        text(x=Ldat$Year[nYears],y=1,"F=M", adj=c(0.8,-0.5), col="darkgreen")

        #---------------------------------------------
        # Plot time series of B/B0 
        #---------------------------------------------
        plot(x=Ldat$Year,y=Ldat$BB0,ylim=c(0,min(c(1.1,max(c(0.6,Ldat$BB0.ucl,1.1*BFM1B0))))), 
             bty="l",type = "l", lwd=1.5, xaxt="n",
             main="exploited B / B0",xlab="Year",ylab="B / B0")
        axis(1,at=Ldat$Year)
        lines(x=Ldat$Year,y=Ldat$BB0.lcl,lty="dotted")
        lines(x=Ldat$Year,y=Ldat$BB0.ucl,lty="dotted")
        abline(h=1.0,col="darkgreen") # B0
        text(x=Ldat$Year[nYears],y=1,"B0", adj=c(0.8,-0.5), col="darkgreen")
        lines(x=Ldat$Year,y=rep(BFM1B0,nYears),lty="dashed", col="darkgreen")
        text(x=Ldat$Year[nYears-1],y=BFM1B0,"B F=M, Lc=opt", adj=c(0.8,-0.5),col="darkgreen")
        lines(x=Ldat$Year,y=rep(BFM1B0/2,nYears),lty="dotted", col="red")
        text(x=Ldat$Year[nYears-1],y=BFM1B0/2,"proxy 0.5 Bmsy", adj=c(0.8,-0.5),col="red")
    }else{
        writeLines("Number of years is not larger than 1.")
    }
    par(opar)
}
