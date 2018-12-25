#' @title Stock simulation
#' 
#' @description  This function estimates stock size, biomass and yield of a stock from
#'    fishing mortality per age class or length group. This function is embedded in the
#'    Thompson and Bell model (prediction model: \code{\link{predict_mod}}).
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \code{age} or \code{midLengths}: midpoints of length classes (length-frequency
#'   data) or ages (age composition data),
#'   \item \code{meanWeight}: mean weight in kg per age class or length group,
#'   \item \code{meanValue}: mean value per kg fish per age class or length group,
#'   \item \code{FM}: fishing mortality rates per age class or length group,
#'   \item \code{M} or \code{Z}: natural or total instantaneous mortality rate.
#' }
#' @param age_unit indicates if the age groups are per month (\code{"month"}) or
#'    per year (\code{"year"}). Default: \code{"year"}
#' @param stock_size_1 stock size of smallest age/length group
#' @param plus_group indicates age/length group, which should be turned into a plus
#'    group (i.e. all groups above are comprised in one group)
#'
#' @keywords function prediction ypr
#'
#' @examples
#' # age-based stock simulation
#' data(shrimps)
#'
#' # option 1: without plus group
#' stock_sim(shrimps, age_unit = "month")
#'
#' # option 2: with plus group
#' stock_sim(param = shrimps, age_unit = "month", plus_group = 11)
#'
#' # length-based stock simulation
#' data(hake)
#'
#' stock_sim(param = hake, stock_size_1 = 98919.3)
#'
#' @details better to treat last group always as a plus group...
#'      if stock size 1 not provided assumes 1000 as intital population size
#'      make sure that FM is also in same unit as the classes, e.g. when classes in
#'      months than also FM has to be provided in 1/months
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{dt}: delta t,
#'   \item \strong{N}: population numbers,
#'   \item \strong{dead}: number of deaths due to natural mortality,
#'   \item \strong{C}: catch,
#'   \item \strong{Y}: yield,
#'   \item \strong{B}: biomass,
#'   \item \strong{V}: value,
#'   \item \strong{totals}: summarised output:
#'   \itemize{
#'   \item \strong{totC} total catch,
#'   \item \strong{totY} total yield,
#'   \item \strong{totV} total value,
#'   \item \strong{meanB} mean biomass.
#'   },
#' }
#'
#' @references
#' Garcia, S. and N.P. van Zalinge, 1982. Shrimp fishing in Kuwait: methodology
#'    for a joint analysis of the artisanal and industrial fisheries. pp. 119-142 In:
#'    Report on the Workshop on assessment of the shrimp stocks of the west coast of
#'    the Gulf between Iran and the Arabian Peninsula. Fisheries development in the
#'    Gulf. Rome, FAO, FI:DP/RAB/80/015/1, 163 p.
#'
#' Millar, R. B., & Holst, R. (1997). Estimation of gillnet and hook selectivity using
#' log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#' 54(3):471-477
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

stock_sim <- function(param, age_unit = "year",
                      stock_size_1 = NA, plus_group = NA){

    res <- param

    meanValue <- res$meanValue

    ## mortalities
    FM <- res$FM
    if(!is.null(res$M)){
        nM <- res$M
        Z <- FM + nM
    }else{
        Z <- res$Z
        nM <- mean(Z - FM,na.rm=T)
    }

    if('midLengths' %in% names(res)) classes <- as.character(res$midLengths)
    if('age' %in% names(res)) classes <- as.character(res$age)
    ##  create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    ##  age based
    if('age' %in% names(res)){

        ## length at maturity
        Lmat <- ifelse(is.null(res$Lmat), NA, res$Lmat)
        wmat <- ifelse(is.null(res$wmat), NA, res$wmat)        
        ## mature individuals
        if(all(c("Linf","K") %in% names(param))){
            mids <- VBGF(param, t = classes.num)
            matP <- 1 / (1 + exp(-(mids - Lmat) /
                                   (wmat / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25))))))
        }else{
            ## writeLines("For estimation of SSB please provide Lmat, wmat and VBGF parameters.")
            matP <- NA
        }        

        ## mean weight
        meanWeight <- res$meanWeight
        
        ## delta t
        dt <- c(diff(classes.num),NA)
        if(age_unit == 'month'){
            dt <- dt * 1/12
        }

        ## population per age group
        N <- rep(NA,length(classes.num))
        N[1] <- ifelse(is.na(stock_size_1),1000,stock_size_1)

        for(x1 in 2:length(N)){
            N[x1] <- N[x1-1] * exp(-Z[x1-1] * dt[x1-1])
            ##  if(x1 == length(N)){
            ##    N[x1] <- N[x1-1] * exp(-Z[x1-1] * dt[x1-1])
            ##  }
        }

        ## number of deaths per time step month or year
        dead <- c(abs(diff(N)),NA)

        ## catch in numbers
        C <- dead * (FM / Z)

        ## yield in kg or g
        Y <- C * meanWeight

        ## mean biomass in kg or g
        B <- Y / (FM * dt)

        ## value expressed in money units
        V <- Y * meanValue

        ## SSB
        SSB <- B * matP

        ## total catch, yield, value and average biomass
        totals <- data.frame(totC = sum(C, na.rm=TRUE),
                             totY = sum(Y, na.rm=TRUE),
                             totV = sum(V, na.rm=TRUE),
                             meanB = sum((B * dt), na.rm=TRUE) / sum(dt, na.rm=TRUE),
                             meanSSB = sum((SSB * dt), na.rm=TRUE) /
                                 sum(dt, na.rm=TRUE))
        ##  more complicated biomass concept if dt is not constant, see Chapter 5

        res2 <- list(dt = dt, N = N, dead = dead, C = C,
                     Y = Y, B = B, SSB = SSB, V = V,
                     totals = totals)

        ## with plus group
        if(!is.na(plus_group)){

            df <- do.call(cbind,res2[1:(length(res2)-1)])
            df <- df[1:plus_group,]
            classes <- classes[1:plus_group]

            ## new class
            classes[length(classes)] <-
                paste(classes[length(classes)],"+",sep='')

            ## num deaths
            df[plus_group, "dead"] <- df[plus_group, "N"]

            ## catch
            new.C <- (FM[plus_group] / Z[plus_group]) * df[plus_group, "N"]
            catch.plus.dif <- new.C - df[plus_group, "C"]
            df[plus_group, "C"] <- new.C

            ## yield
            df[plus_group, "Y"] <- meanWeight[plus_group] * catch.plus.dif

            ## value
            df[plus_group, "V"] <-
                df[plus_group, "Y"] * meanValue[plus_group]

            ## biomass       ####not sure....omitted in manual
            df[plus_group, "B"] <-
                df[plus_group, "Y"] / (FM[plus_group] * df[plus_group, "dt"])

            df[plus_group, "SSB"] <- df[plus_group, "B"] * matP[plus_group]
            

            res2 <- as.list(as.data.frame(df))

            df2 <- do.call(cbind,res)
            df2 <- df2[1:plus_group,]
            res <- as.list(as.data.frame(df2))
            res$age <- classes

            ## total catch, yield, value and average biomass
            totals <- data.frame(totC = sum(res2$C, na.rm=TRUE),
                                 totY = sum(res2$Y, na.rm=TRUE),
                                 totV = sum(res2$V, na.rm=TRUE),
                                 meanB = sum((res2$B * res2$dt), na.rm=TRUE) / sum(dt, na.rm=TRUE),
                                 meanSSB = sum((res2$SSB * res2$dt), na.rm=TRUE) / sum(dt, na.rm=TRUE))
            ## more complicated biomass concept if dt is not constant, see Chapter 5

            res2 <- c(res2,totals = list(totals))
        }
    }

    ##  length based
    if('midLengths' %in% names(res)){

        if("par" %in% names(res)){
            Winf <- ifelse("Winf" %in% names(res$par), res$par$Winf, NA)
            Linf <- ifelse("Linf" %in% names(res$par), res$par$Linf, NA)
            K <- res$par$K
            t0 <- ifelse("t0" %in% names(res$par), res$par$t0, 0)
            C <- ifelse("C" %in% names(res$par), res$par$C, 0)
            ts <- ifelse("ts" %in% names(res$par), res$par$ts, 0)
            ## maturity parameters            
            Lmat <- ifelse("Lmat" %in% names(res$par), res$par$Lmat, NA)
            wmat <- ifelse("wmat" %in% names(res$par), res$par$wmat, NA)     
        }else{
            Winf <- ifelse("Winf" %in% names(res), res$Winf, NA)
            Linf <- ifelse("Linf" %in% names(res), res$Linf, NA)
            K <- res$K
            t0 <- ifelse("t0" %in% names(res), res$t0, 0)
            C <- ifelse("C" %in% names(res), res$C, 0)
            ts <- ifelse("ts" %in% names(res), res$ts, 0)
            ## maturity parameters
            Lmat <- ifelse("Lmat" %in% names(res), res$Lmat, NA)
            wmat <- ifelse("wmat" %in% names(res), res$wmat, NA)            
        }

        a <- res$a
        b <- res$b

        ## mature individuals
        mids <- classes.num
        matP <- 1 / (1 + exp(-(mids - Lmat) /
                             (wmat / ( log(0.75/(1-0.75)) - log(0.25/(1-0.25))))))

        ## calculate size class interval
        int <- classes.num[2] - classes.num[1]

        ## t of lower and upper length classes
        lowL <- classes.num - (int / 2)
        upL <- classes.num + (int/2)

        ## H
        H <- ((Linf - lowL)/(Linf - upL)) ^ (nM/(2*K))

        ## population
        N <- rep(NA,length(classes.num))
        N[1] <- ifelse(is.na(stock_size_1), 1000, stock_size_1)
        for(x1 in 2:length(classes.num)){
            N[x1] <- N[x1-1] * ((1/H[x1-1]) - (FM[x1-1]/Z[x1-1])) /
                (H[x1-1] - (FM[x1-1]/Z[x1-1]))
        }

        ## number of deaths per time step month or year
        dead <- c(abs(diff(N)),NA)

        ## catch
        C <- dead * (FM/Z)

        ## average weight
        W <- a * ((lowL + upL)/2 ) ^ b

        ## yield
        Y <- C * W

        ## Value
        V <- Y * meanValue

        ## biomass
        B <- (dead / Z ) * W

        ## SSB
        SSB <- B * matP

        ## last length group
        x2 <- length(classes.num)
        C[x2] <- N[x2] * FM[x2]/Z[x2]
        W[x2] <- a * ((lowL[x2] + Linf)/2 ) ^ b
        Y[x2] <- C[x2] * W[x2]
        B[x2] <- N[x2] / Z[x2] * W[x2]
        V[x2] <- Y[x2] * meanValue[x2]
        SSB[x2] <- B[x2] * matP[x2]


        ## total catch, yield, value and average biomass
        totals <- data.frame(totC = sum(C, na.rm=TRUE),
                             totY = sum(Y, na.rm=TRUE),
                             totV = sum(V, na.rm=TRUE),
                             meanB = sum(B, na.rm=TRUE),
                             meanSSB = sum(SSB, na.rm=TRUE))

        res2 <- list(N = N,
                     dead = dead,
                     C = C,
                     Y = Y,
                     B = B,
                     SSB = SSB,
                     V = V,
                     totals = totals)
    }

    ret <- c(res,res2)
    return(ret)
}

