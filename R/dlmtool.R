#' @title Harvest control rules based on TropFishR
#'
#' @description This function creates an harvest control rule for use with the
#'   DLMtool package.
#' 
#' @param uncertaintyCap Logical; If true TAC is bound between two values set in lower and upper. Default: FALSE.
#' @param lower lower bound of the uncertainty cap. Default is 0.8, used if uncertaintyCap = TRUE.
#' @param upper upper bound of the uncertainty cap. Default is 1.2, used if uncertaintyCap = TRUE.
#'
#' @return A function which can estimate TAC recommendations based on SPiCT assessment,
#'   taking assessment uncertainty into account.
#' 
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(DLMtool)
#' 
#' ## Put together an operating model from the available DLM toolkit examples
#' StockEx <- Red_snapper
#' FleetEx <- Generic_IncE
#' ObsEx <- Precise_Unbiased
#' 
#' ## Remove changes in life history parameters
#' StockEx@Mgrad <- c(0,0)
#' StockEx@Kgrad <- c(0,0)
#' StockEx@Linfgrad <- c(0,0)
#' StockEx@Prob_staying <- c(1,1)
#' 
#' ## Set the depletion level 
#' StockEx@D <- c(0.3, 0.4)
#'
#' ## create Operation Model
#' OMex <- new("OM", Stock = StockEx, Fleet = FleetEx, 
#'                   Obs = ObsEx)
#' 
#' ## Set simulation options
#' OMex@nsim <- 10
#' OMex@nyears <- 25
#' OMex@proyears <- 5
#' 
#' ## Get TFR HCR
#' MPname <- TFR2DLMtool(uncertaintyCap = c(FALSE,TRUE))
#' MPname
#' 
#' ## run MSE
#' MSEex <- runMSE(OMex, MPs = MPname, CheckMPs = FALSE)
#' Pplot2(MSEex, traj="quant", quants=c(0.2, 0.8))
#' }
#'
#'
TFR2DLMtool <- function(uncertaintyCap = FALSE,
                        lower=0.8,
                        upper=1.2,
                        env=globalenv()){

    ## allowing for multiple generation of MPs
    argList <- list(uncertaintyCap, lower, upper)
    argLengths <- sapply(argList, length)
    maxi <- max(argLengths)
    maxl  <- which(argLengths == maxi)
    if(maxi>1){
        if(max(argLengths[(1:3)[-maxl]]) > 1)
            stop("Specified arguments have different lengths, they should have the same length or length = 1.")
    }
    argListCor <- argList
    argListCor[(1:3)[-maxl]] <- lapply(argList[(1:3)[-maxl]], function(x) rep(unlist(x), maxi))
    uncertaintyCapPrint <- argListCor[[1]]
    uncertaintyCapPrint[which(uncertaintyCapPrint == TRUE)] <- "T"
    uncertaintyCapPrint[which(uncertaintyCapPrint == FALSE)] <- "F"                        

    ## dimensions Data@CAL: nsim x nyears x length(CAL_bins)
    ## MP as function
    template <- expression(paste0(
    'structure(function(x, Data, reps=100,
          uncertaintyCap=',a,',
          lower=',b,',
          upper=',c,'){

              ## input from OM
              dependencies = "Data@Cat, Data@vbK, Data@CV_vbK, Data@Mort, Data@CV_Mort, Data@wla, Data@CV_wla, Data@wlb, Data@CV_wlb, Data@CAL, Data@CAL_bins, Data@vbLinf, Data@CV_vbLinf, Data@FMSY_M, Data@vbt0, Data@CV_vbt0"

              ## midlengths
              interval <- Data@CAL_bins[3] - Data@CAL_bins[2]  ## assuming constant interval (removing first lower bound == 0)
              midLengths <- Data@CAL_bins[-1] + interval/2
              ## sample LH parameters
              linfs <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
              ks <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
              t0s <- trlnorm(reps, Data@vbt0[x], Data@CV_vbt0[x])
              ms <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])
              wlas <- trlnorm(reps, Data@wla[x], Data@CV_wla[x])
              wlbs <- trlnorm(reps, Data@wlb[x], Data@CV_wlb[x])
              fmsys <- rep(Data@FMSY_M[x] * Data@Mort[x], reps)
              if(Data@vbt0[x] != 0 & Data@CV_vbt0[x] != 1e-9){
                   t0s <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
              }else{
                   t0s <- rep(Data@vbt0[x], reps)
              }
              t0s[!is.finite(t0s)] <- 0
              ## run TFR FSA
              res <- mapply(function(linf, k, t0, m, wla, wlb, fmsy){
                            ## create lfq data
                            lfq <- list(midLengths = midLengths,
                                        dates = as.Date("2000-01-01"),  ## not relevant
                                        catch = as.matrix(Data@CAL[x, dim(Data@CAL)[2], ]))
                            class(lfq) <- "lfq"
                            lfq <- lfqModify(lfq, par=list(Linf=linf,K=k,t0=t0), vectorise_catch = TRUE)
                            ## create selectivity
                            ##   slist <- list(selecType = trawl_ogive, L50 = l50s, L75 = l75s)
                            ffmsy <- TropFishR::TFSA(lfq, m, wla, wlb, fmsy)
                            if(length(ffmsy) == 0) ffmsy <- NA
                            return(ffmsy)
              }, linfs, ks, t0s, ms, wlas, wlbs, fmsys, SIMPLIFY = TRUE)
              ## apply HCR
              Cc <- trlnorm(reps, Data@Cat[x, length(Data@Cat[x, ])], Data@CV_Cat[x])
              r <- 1         ## account for trend in stock biomass
              f <- 1 / res   ## proxy for fmsy/current exploitation
              b <- 1         ## min(1, proxy for B/Btrigger)
              if(uncertaintyCap){
                  TAC <- Cc * sapply(r*f*b, function(x) min(upper, max(lower, x)))
              }else{
                  TAC <- Cc * r * f * b
              }
              ## create output in required format
              Rec <- new("Rec")
              Rec@TAC <- DLMtool::TACfilter(TAC)
              return(Rec)
        },
        class="MP")'))
        
    nami <- rep(NA,maxi)
    for(I in 1:maxi){

        ## create MPs as functions
        subList <- lapply(argListCor, "[[", I)
        names(subList) <- letters[1:3]
        templati <- eval(parse(text=paste(parse(text = eval(template, subList)),collapse=" ")))

        ## save names of MPs
        if(uncertaintyCap){
            nami[I] <- paste0("TFR_uC_",uncertaintyCapPrint[I],
                              "_lo_",argListCor[[2]][I],
                              "_up_",argListCor[[3]][I])            
        }else{
            nami[I] <- paste0("TFR_uC_",uncertaintyCapPrint[I])
        }

        assign(value=templati, x=nami[I], envir=env)
    }

    ## allow for assigning names
    invisible(nami)
}








#' @title Traditional Fish Stock Assessment with TropFishR
#' 
#' @description Estimate F/F01 based on the length converted catch curve, length-based VPA and
#'   the Thompson and Bell yield per recruit analysis
#' 
#' @param lfq lfq data
#' @param m natural mortality
#' @param wla coefficient from length-weight relationship
#' @param wlb exponent from length-weight relationship
#' 
#' @return F/F01 estimated by LCCC + length-based VPA + YPR
#' 
#' @export
#'
#'
TFSA <- function(lfq, m, wla, wlb, fmsy){   ## import slist
    ## catch curve
    rescc <- try(TropFishR::catchCurve(lfq,
                                       auto = TRUE,
                                       plot = FALSE),
                 silent = TRUE)
    if(is(rescc, "try-error")){
        return(NA)
    }else{
        z <- rescc$Z
        fm <- z - m
##        print(paste0("F:  ", fm, " -- "))
        if(is.na(fm) | length(fm) == 0 | fm < 0){
            return(NA)
        }else{
            if(lfq$par$Linf < max(lfq$midLengths)){
                interval <- lfq$midLengths[2] - lfq$midLengths[1]
                upperLength <- lfq$midLengths + (interval / 2)
                lfq <- lfqModify(lfq,
                                 plus_group = lfq$midLengths[which.min(abs(upperLength - floor(lfq$par$Linf)))])
                plusgroup <- TRUE
            }else{
                plusgroup <- FALSE
            }
            lfq$a <- wla
            lfq$b <- wlb
            lfq$M <- m
            ## vpa
            resvpa <- try(TropFishR::VPA(lfq, terminalF=fm,
                                         plus_group = plusgroup,
                                         plot = FALSE),
                          silent = TRUE)
            if(is(resvpa, "try-error")){
                return(NA)
            }else{
                fvec <- resvpa$FM_calc
                lfq$FM <- fvec
                ## ypr
                resypr <- try(predict_mod(lfq,
                                          type="ThompBell",
##                                          s_list = slist,
                                          plus_group = plusgroup,
                                          FM_change=seq(0.001,3,0.01),
                                          plot = FALSE),
                              silent = TRUE)
                if(is(resypr, "try-error")){
                    return(NA)
                }else{
                    f01 <- as.numeric(resypr$df_Es[which(names(resypr$df_Es) == "F01")])  ## F01 as proxy for Fmsy
                    ff01 <- max(fvec)/f01
##                    print(paste0("Ref:  ",f01," -- ", fmsy))
                    return(ff01)                
                }
            }            
        }

    }
}
    
