#' @title Catch curve
#'
#' @description  This function applies the (length-converted) linearised catch
#'    curve to age composition and length-frequency data,
#'    respectively. It allows to estimate the instantaneous total mortality rate (Z).
#'    Optionally, the gear selectivity can be estimated and the cumulative catch
#'    curve cna be applied.
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \code{midLengths} or \code{age}: midpoints of the length classes (length-frequency
#'   data) or ages (age composition data),
#'   \item \code{Linf}: infinite length for investigated species in cm [cm],
#'   \item \code{K}: growth coefficent for investigated species per year [1/year],
#'   \item \code{t0}: theoretical time zero, at which individuals of this species hatch,
#'   \item \code{catch}: catches, vector or matrix with catches of subsequent years if
#'   the catch curve with constat time intervals should be applied;
#' }
#' @param catch_columns numerical; indicating the column of the catch matrix which should be
#'   used for the analysis.
#' @param cumulative logical; if TRUE the cumulative
#'   catch curve is applied (Jones and van Zalinge method)
#' @param calc_ogive logical; if TRUE the selection ogive is additionally
#'   calculated from the catch curve (only if \code{cumulative = FALSE})
#' @param reg_int instead of using the identity method a range can be determined,
#'    which is to be used for the regression analysis. If equal to NULL identity method
#'    is applied (default). For multiple regression lines provide list with the two points
#'    for the regression line in each element of the list.
#' @param reg_num integer indicating how many separate regression lines should be applied to the
#'    data. Default 1.
#' @param auto logical; no interactive functions used instead regression line is chosen
#'    automatically. Default = FALSE
#' @param boot an object of class 'lfqBoot'
#' @param natMort optional argument; name of column with natural mortality estimates
#'    for application to 'lfqBoot' object
#' @param robustReg logical; indicating whether the robust automatic fitting of the regression line should be used
#'    for application to 'lfqBoot' object
#' @param plot logical; should a plot be displayed? Default = TRUE
#' 
#'
#' @keywords function mortality Z catchCurve
#'
#' @examples
#' \donttest{
#' #_______________________________________________
#' # Variable paramter system (with catch vector)
#' # based on length frequency data
#' data(goatfish)
#' output <- catchCurve(goatfish)
#' summary(output$linear_mod)
#'
#' # based on age composition data
#' data(whiting)
#' catchCurve(whiting, catch_columns = 1)
#'
#' #_______________________________________________
#' # Constant parameter system based on age composition data (with catch matrix)
#' catchCurve(whiting)
#'
#' #_______________________________________________
#' # Cumulative Catch Curve
#' # based on length frequency data
#' data(goatfish)
#' catchCurve(goatfish, cumulative = TRUE)
#'
#' # based on age composition data
#' data(synCAA2)
#' catchCurve(synCAA2, cumulative = TRUE)
#'
#' #_______________________________________________
#' # Catch Curve with estimation of selection ogive
#' data(synLFQ3)
#' output <- catchCurve(synLFQ3, calc_ogive = TRUE)
#' summary(output$linear_mod_sel)
#'
#' #_______________________________________________
#' # Catch curve with results from bootstrapped ELEFAN
#' # coming soon
#'  }
#'
#' # the same with predefined selection for regression line:
#' output <- catchCurve(synLFQ3, calc_ogive = TRUE, reg_int = c(9,21))
#' plot(output, plot_selec = TRUE)
#'
#' @details This function includes the \link{identify} function, which asks you to
#'   choose two points from a graph manually. The two points which you choose by clicking
#'   on the plot in the graphical device represent the start and end of the data points,
#'   which should be used for the analysis. Based on these points the regression line
#'   is calculated.
#'   When the selection ogive
#'   is calculated by means of the catch curve the assumption is made, that Z is constant
#'   for all year classes or length groups, respectively. Accoring to Sparre and Venema
#'   (1998) this assumption might be true, because F is smaller for young fish
#'   (Selectivity) while M is higher for young fish (high natural mortality). The selectivity
#'   for not fully exploited old fish (e.g. due to gillnet fishery) can not be calculated yet
#'   by use of the catch curve.
#'   Based on the format of the list argument \code{catch} and whether the argument
#'   \code{catch_columns} is defined, the function automatically
#'   distinguishes between the catch curve with variable parameter system (if catch is a
#'   vector) and the one with constant parameter system (if catch is a matrix or a
#'   data.frame and \code{catch_columns = NA}). In the case of the variable parameter
#'   system the catches of one year are
#'   assumed to represent the catches during the entire life span of a so called
#'   pseudo-cohort.
#'   The cumulative catch curve does not allow for the estimation of the selectivity
#'   ogive.
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{classes.num}, \strong{tplusdt_2}, \strong{t_midL}, or
#'      \strong{ln_Linf_L}: age, relative age or subsitute depending on input and method,
#'   \item \strong{lnC} or \strong{lnC_dt}: logarithm of (rearranged) catches,
#'   \item \strong{reg_int}: the interval used for the regression analysis,
#'   \item \strong{linear_mod}: linear model used for the regression analysis,
#'   \item \strong{Z}: instantaneous total mortality rate, confidenceInt
#'   \item \strong{se}: standard error of the total mortality;
#'   \item \strong{confidenceInt}: confidence interval of the total mortality;}
#' in case calc_ogive == TRUE, additionally:
#' \itemize{
#'   \item \strong{intercept}: intercept of regression analysis,
#'   \item \strong{linear_mod_sel}: linear model used for the selectivity analysis,
#'   \item \strong{Sobs}: observed selection ogive,
#'   \item \strong{ln_1_S_1}: dependent variable of regression analysis for
#'   selectivity parameters,
#'   \item \strong{Sest}: estimated selection ogive,
#'   \item \strong{t50}: age at first capture (age at which fish have a 50%
#'   probability to be caught),
#'   \item \strong{t75}: age at which fish have a 75% probability to be caught,
#'   \item \strong{L50}: length at first capture (length at which fish have a 50%
#'   probability to be caught),
#'   \item \strong{L75}: length at which fish have a 75% probability to be caught;
#' }
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics identify par plot
#' @importFrom stats lm na.omit
#' @importFrom utils flush.console
#'
#' @references
#' Baranov, F.I., 1926. On the question of the dynamics of the fishing industry.
#' \emph{Nauchn. Byull. Rybn. Khoz}, 8 (1925), 7-11
#'
#' Beverton, R.J.H. and S.J. Holt, 1956. A review of methods for estimating mortality
#' rates in exploited fish populations, with special reference to sources of bias in
#' catch sampling. \emph{Rapports et Proces verbaux des Reunions}, Conseil Table3
#'
#' Chapman, D., and D.S Robson, 1960. The analysis of a catch curve.
#' \emph{Biometrics}, 354-368
#'
#' Edser, T., 1908. Note on the number of plaice at each length, in certain samples
#' from the southern part of the North Sea, 1906. \emph{Journal of the Royal
#' Statistical Society}, 686-690
#'
#' Heincke, F., 1913. Investigations on the plaice. General report. 1. The plaice fishery
#' and protective regulations. Part I. \emph{Rapp.P.-v.Reun.CIEM}, 17A:1-153 + Annexes
#'
#' ICES, 1981. Report of the \emph{Ad hoc} working group on the use of effort data in
#' assessment, Copenhagen, 2-6 March 1981. \emph{ICES C.M.} 1981/G:5 (mimeo)
#'
#' Jones, R., and N.P. Van Zalinge, 1981. Estimates of mortality rate and population size
#' for shrimp in Kuwait waters. \emph{Kuwait Bull. Mar. Sci}, 2, 273-288
#'
#' Pauly, D., 1983. Length-converted catch curves: a powerful tool for fisheries research
#' in the tropics (part I). \emph{ICLARM Fishbyte}, 1(2), 9-13
#'
#' Pauly, D., 1984. Length-converted catch curves: a powerful tool for fisheries
#' research in the tropics (part II). \emph{ICLARM Fishbyte}, 2(1), 17-19
#'
#' Pauly, D., 1984. Length-converted catch curves: a powerful tool for fisheries
#' research in the tropics (III: Conclusion). \emph{ICLARM Fishbyte}, 2(3), 9-10
#'
#' Ricker, W.E., 1987. Computation and interpretation of biological statistics of fish
#' populations. \emph{Bull.Fish.Res.Board Can.}, (191):382 p.
#'
#' Robson, D.S., and D.G. Chapman, 1961. Catch curves and mortality rates.
#' \emph{Trans.Am.Fish.Soc.}, 90(2):181-189
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' Van Sickle, J. 1977. Mortality rates from size distributions: the application of a
#' conservation law. \emph{Oecologia, Berl.}, 27(4):311-318
#'
#' @export

catchCurve <- function(param,
                       catch_columns = NA,
                       cumulative = FALSE,
                       calc_ogive = FALSE,
                       reg_int = NULL,
                       reg_num = 1,
                       auto = FALSE,
                       boot = NULL,
                       natMort = NULL,
                       robustReg = FALSE,
                       plot = TRUE
                       ){

    ## catch curve with bootstrapping ELEFAN results
    if(!is.null(boot) & class(boot) == "lfqBoot"){

        bootOut <- boot
        bootRaw <- boot$bootRaw


        if (any(is.null(bootRaw$Linf), is.null(bootRaw$K)))
          stop("LCCC with boot requires a boot object with columns Linf and K")
        

        Zs <- vector("numeric",nrow(bootRaw))
        ## seZ <- vector("numeric",nrow(bootRaw))
        ## confZ <- vector("numeric",nrow(bootRaw))
        t50s <- vector("numeric", nrow(bootRaw))
        t75s <- vector("numeric", nrow(bootRaw))
        L50s <- vector("numeric", nrow(bootRaw))
        L75s <- vector("numeric", nrow(bootRaw))        
        
        for(bi in 1:nrow(bootRaw)){

            set.seed(boot$seed[bi])
            lfqTemp <- lfqPermutate(param)
            lfqLoop <- lfqModify(lfqTemp, vectorise_catch = TRUE)

            ## error if lfq data spans several years!
            if(class(lfqLoop$catch) == "matrix") stop("The lfq data spans several years, please subset for one year at a time!")

            catch <- lfqLoop$catch

            classes <- as.character(lfqLoop$midLengths)
            classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
            classes.num <- as.numeric(classes.num[,1])            

            ##calculate size class interval
            midLengths <- classes.num
            interval <- midLengths[2] - midLengths[1]            

            ## L and t of lower length classes
            lowerLengths <- midLengths - (interval / 2)
            if("C" %in% names(bootRaw) & "ts" %in% names(bootRaw)){
                t_L1 <- VBGF(param = list(Linf = bootRaw$Linf[bi], K = bootRaw$K[bi],
                                          t0 = 0, C=bootRaw$C[bi], ts=bootRaw$ts[bi]), L = lowerLengths)
            }else{
                t_L1 <- VBGF(param = list(Linf = bootRaw$Linf[bi], K = bootRaw$K[bi],
                                          t0 = 0), L = lowerLengths)
            }
            ## t0 - (1/K) * log(1 - (lowerLengths / Linf))

            ## delta t
            dt <- rep(NA,length(midLengths))
            for(x1 in 1:(length(dt)-1)){
              dt[x1] <- t_L1[x1+1] - t_L1[x1]
            }

            # x varaible
            #ln (Linf - L)
            ln_Linf_L <- log(bootRaw$Linf[bi] - lowerLengths)
            ## t of midlengths
            if("C" %in% names(bootRaw) & "ts" %in% names(bootRaw)){
                t_midL <- VBGF(param = list(Linf = bootRaw$Linf[bi], K = bootRaw$K[bi],
                                            t0 = 0, C=bootRaw$C[bi], ts=bootRaw$ts[bi]), L = midLengths)
            }else{
                t_midL <- VBGF(param = list(Linf = bootRaw$Linf[bi], K = bootRaw$K[bi],
                                            t0 = 0), L = midLengths)
            }
            ## t0 - (1/K) * log(1 - (midLengths / Linf))

            # y variable
            #ln C(L1,Linf)
            lnC <- log(catch)
            # ln( Catch / delta t)
            lnC_dt <- log(catch / dt)
            lnC_dt[which(lnC_dt == -Inf)] <- NA   ### OR zero???

            if(!cumulative){
              xvar = t_midL
              yvar = lnC_dt
            }
            if(cumulative){
              xvar = ln_Linf_L
              yvar = lnC
            }

            
            ## remove all NAs and Infs
            temp <- cbind(xvar,yvar)
            temp <- as.matrix(na.exclude(temp))
            temp <- temp[(!(temp[,1] == Inf | temp[,1] == -Inf)),]
            temp <- temp[(!(temp[,2] == Inf | temp[,2] == -Inf)),]
            xvar <- temp[,1]
            yvar <- temp[,2]


            ## improving automatic fitting of regression line:
            if(robustReg){
                if(FALSE){
                ## 1: if last point is outlier (e.g. due to gear selectivity) -> remove
                if(xvar[which.max(xvar)] > (xvar[which.max(xvar)-1]+5)){
                    yvar[which.max(xvar)] <- NA
                    xvar[which.max(xvar)] <- NA
                }
                }
                
                ## 2: fitting regression line through last 3 points and compare slopes
                ##      -> if very different selectivity might be responsible
                ##      -> then remove last point
                yvar2 <- as.numeric(yvar)
                xvar2 <- xvar[which(yvar2 > 0.2)]
                
                indX <- (which.max(xvar2)-2):which.max(xvar2)
                yvarX <- yvar2[indX]
                xvarX <- xvar2[indX]
                sx <- abs(coefficients(lm(yvarX ~ xvarX))[2])

                indXX <- (which.max(yvar2)+1) : which.max(xvar2)
                xvarXX <- xvar2[indXX]
                yvarXX <- yvar2[indXX]
                sxx <- abs(coefficients(lm(yvarXX ~ xvarXX))[2])

                if(sx/sxx < 0.8){
                    yvar[which.max(xvar)] <- NA
                    xvar[which.max(xvar)] <- NA                    
                }

                if(FALSE){
                    plot(xvar2,yvar2)
                    abline(lm(yvarXX ~ xvarXX))
                    abline(lm(yvarX ~ xvarX),col=4)
                }

                ## remove all NAs and Infs
                temp <- cbind(xvar,yvar)
                temp <- as.matrix(na.exclude(temp))
                temp <- temp[(!(temp[,1] == Inf | temp[,1] == -Inf)),]
                temp <- temp[(!(temp[,2] == Inf | temp[,2] == -Inf)),]
                xvar <- temp[,1]
                yvar <- temp[,2]
            }

            

            ## cut
            yvar2 <- as.numeric(yvar)
            xvar2 <- xvar[which(yvar2 > 0.2)]
            cutter <- c(which(yvar2 == max(as.numeric(yvar2),na.rm=TRUE))+1, which(xvar2 == max(xvar2,na.rm=TRUE)))

            ## calculations + model
            df.CC <- as.data.frame(cbind(xvar,yvar))
            df.CC.cut <- df.CC[cutter[1]:cutter[2],]
            lm1 <- try(lm(yvar ~ xvar, data = df.CC.cut), silent = TRUE)
            if(class(lm1) != "try-error"){
                sum_lm1 <- summary(lm1)
                r_lm1 <- sum_lm1$r.squared
                intercept_lm1 <- sum_lm1$coefficients[1]
                slope_lm1 <- sum_lm1$coefficients[2]
                se_slope_lm1 <- sum_lm1$coefficients[4]

                ## fit of regression line
                lm1.fit <- sum_lm1$r.squared
                Z_lm1 <- abs(slope_lm1)
                SE_Z_lm1 <- abs(se_slope_lm1)
                confi <-  abs(se_slope_lm1) * qt(0.975,sum_lm1$df[2])
                conf_Z_lm1 <- Z_lm1 + c(-confi,confi)

                ## special case when cumulative and length-frequency data
                if(cumulative & "midLengths" %in% names(param) == TRUE){
                  Z_lm1 <- Z_lm1 * K
                  SE_Z_lm1 <- SE_Z_lm1 * K
                }

                Zs[bi] <- Z_lm1
                ## seZ[bi] <- SE_Z_lm1
                ## confZ[bi] <- conf_Z_lm1
                if(calc_ogive){
                  ## Assumption that Z of smallest selected individuals is most appropriate
                  mini <- min(cutter)

                  ## only use part of catch and t which is not fully exploited by the gear
                  t_ogive <- xvar[1:(cutter[1]-1)]
                  dt_ogive <- dt[1:(cutter[1]-1)]
                  catch_ogive <- catch[1:(cutter[1]-1)]

                  ## it could be that the smallest length is already fully exploited
                    if(length(catch_ogive) < 2){
                      t50s[bi] <- t75s[bi] <- t_ogive ## assuming knife edge
                      t0 <- 0
                      L50s[bi] <- L75s[bi] <- bootRaw$Linf[bi]*(1-exp(-bootRaw$K[bi]*(t50s[bi]-t0)))
                  }else{
                      ## calculate observed selection ogive
                      Sobs <- catch_ogive/(dt_ogive *
                                           exp(unlist(intercept_lm1) -
                                               unlist(Z_lm1) * t_ogive))

                      # dependent vairable in following regression analysis
                      ln_1_S_1 <- log((1/Sobs) - 1)

                      # get rid of Inf
                      ln_1_S_1[which(ln_1_S_1 == Inf)] <- NA
                      t_ogive[which(t_ogive == Inf)] <- NA

                      #regression analysis to caluclate T1 and T2
                      mod_ogive <- lm(ln_1_S_1 ~ t_ogive, na.action = na.omit)
                      sum_lm_ogive <- summary(mod_ogive)
                      T1 <- sum_lm_ogive$coefficients[1]
                      T2 <- abs(sum_lm_ogive$coefficients[2])

                      # calculate estimated selection ogive
                      Sest <- 1/(1+exp(T1 - T2*xvar))

                      # selection parameters
                      t50s[bi] <- T1/T2
                      t75s[bi] <- (T1 + log(3))/T2
    ##                  t95 <-  (T1 - log((1 / 0.95) - 1)) / T2
                      t0 <- 0                  
                      L50s[bi] <- bootRaw$Linf[bi]*(1-exp(-bootRaw$K[bi]*(t50s[bi]-t0)))
                      L75s[bi] <- bootRaw$Linf[bi]*(1-exp(-bootRaw$K[bi]*(t75s[bi]-t0)))
    ##                  L95 <- bootRaw$Linf[bi]*(1-exp(-bootRaw$K[bi]*(t95-t0)))
                  }
                }
            }else{
                Zs[bi] <- NA
                if(calc_ogive){
                      L50s[bi] <- NA
                      L75s[bi] <- NA
                }
            }

        }
        bootRaw[,(ncol(bootRaw)+1)] <- Zs
        colnames(bootRaw) <- c(colnames(bootRaw)[-ncol(bootRaw)],"Z")        
        
        ## estimate FM if M in data frame
        if(!is.null(natMort) & natMort %in% colnames(bootRaw)){
            bootRaw[,(ncol(bootRaw)+1)] <- bootRaw$Z - bootRaw[,(colnames(bootRaw) == natMort)]
            colnames(bootRaw) <- c(colnames(bootRaw)[-ncol(bootRaw)],"FM")            
            tmp <- as.data.frame(bootRaw[,(ncol(boot$bootRaw)+(1:2))])
            nx <- 2
        }else{
            tmp <- as.data.frame(bootRaw[,(ncol(boot$bootRaw)+1)])
            nx <- 1
        }
        ## estimate L50 and L75 if calc_ogive
        if(calc_ogive){
            bootRaw[,(ncol(bootRaw)+1)] <- L50s
            bootRaw[,(ncol(bootRaw)+1)] <- L75s            
            colnames(bootRaw) <- c(colnames(bootRaw)[-(((ncol(bootRaw)-1):ncol(bootRaw)))],"L50","L75")
            tmp <- cbind(tmp,as.data.frame(bootRaw[,((ncol(bootRaw)-1):ncol(bootRaw))]))
            nx <- nx + 2
        }

        ## max density and CIS
        resMaxDen <- vector("numeric", ncol(tmp))
        ciList <- vector("list", ncol(tmp))
        for(i in seq(nx)){
            ## max densities
            x <- try(ks::kde(as.numeric(na.omit(tmp[,i]))), TRUE)
            if(class(x) != "try-error"){
                ind <- which(x$estimate > x$cont["99%"])
                resMaxDen[i] <- mean(x$eval.points[ind])
                ## confidence intervals
                CItxt <- paste0(round(5), "%")
                inCI <- rle( x$estimate > x$cont[CItxt] )
                start.idx <- c(1, cumsum(inCI$lengths[-length(inCI$lengths)])+1)
                end.idx <- cumsum(inCI$lengths)
                limCI <- range(x$eval.points[start.idx[min(which(inCI$values))]:end.idx[max(which(inCI$values))]])
                ciList[[i]] <- limCI                
            }else{
                if(length(unique(as.character(tmp[,i]))) == 1 && all(!is.na(tmp[,i]))){
                    resMaxDen[i] <- unique(tmp[,i])
                    ciList[[i]] <- NA
                }else{
                    resMaxDen[i] <- NA
                    ciList[[i]] <- NA
                }
            }
        }
        resCIs <- cbind(boot$bootCIs,t(do.call(rbind,ciList)))
        colnames(resCIs) <- colnames(bootRaw)
        rownames(resCIs) <- c("lo","up")
        resMaxDen <- c(boot$bootMaxDen, resMaxDen)
        names(resMaxDen) <- names(bootRaw)

        ret <- list()
        ret$bootRaw <- bootRaw
        ret$bootMaxDen <- resMaxDen
        ret$bootCIs <- resCIs
        ret$seed <- boot$seed

        class(ret) <- "lfqBoot"
        return(ret)

    }else if(!is.null(boot) & class(boot) != "lfqBoot"){
        stop("You provided an object for boot, but it does not have class 'lfqBoot'. Please check.")
    }else{

        res <- param

        if("midLengths" %in% names(res)) classes <- as.character(res$midLengths)
        if("age" %in% names(res)) classes <- as.character(res$age)
        # create column without plus group (sign) if present
        classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
        classes.num <- as.numeric(classes.num[,1])

        constant_dt <- FALSE
        if(is.na(catch_columns[1]) & (class(res$catch) == 'matrix' |
           class(res$catch) == 'data.frame')){
          writeLines("Please be aware that you provided the catch as a matrix without specifiying any columns for \n the analysis. In this case the methods applies by default the catch curve with constant \n parameter system (refer to the help file for more information).")
          flush.console()
          constant_dt <- TRUE
          catch <- res$catch
        }

        # non cumulative catch curve
        if(cumulative == FALSE){
          if(is.na(catch_columns[1]) & constant_dt == FALSE) catch <- res$catch
          if(!is.na(catch_columns[1])){
            catchmat <- res$catch[,(catch_columns)]
            if(length(catch_columns) > 1){
              catch <- rowSums(catchmat, na.rm = TRUE)
            }else catch <- catchmat
          }
        }
        #  cumulative catch curve
        if(cumulative){
          if(is.na(catch_columns[1])) catch <- rev(cumsum(rev(res$catch)))
          if(!is.na(catch_columns[1])){
            catchmat <- res$catch[,(catch_columns)]
            if(length(catch_columns) > 1){
              catchpre <- rowSums(catchmat, na.rm = TRUE)
            }else catchpre <- catchmat
            catch <- rev(cumsum(rev(catchpre)))
            #catch <- rev(cumsum(rev(res$catch[,(catch_columns)])))
          }
        }

        # Error message if catch and age do not have same length
        #   Linearised catch curve with constant time intervals
        if(constant_dt){
          if("midLengths" %in% names(res) == TRUE) stop(noquote(
            "The catch curve with constant time interval is not applicable to length-frequency data. Please provide a catch vector."))

          #if(length(classes) != length(catch[,1])) stop(noquote(
          #  "Age/length classes and catch matrix do not have the same length!"))

          if(length(classes) != length(diag(as.matrix(catch)))) writeLines("Age/length classes and the real cohort in the catch matrix \ndo not have the same length. The missing age/length \nclasses will be omitted.")

          # Aged based Catch curve
          if("age" %in% names(res) == TRUE){
            #find cohort to analyse
            real.cohort <- diag(as.matrix(catch))
            catch <- c(real.cohort, rep(NA,length(classes.num) - length(real.cohort)))
          }

        }else if(class(catch) == 'numeric'){
          if(length(classes) != length(catch)) stop(noquote(
            "Age/length classes and catch vector do not have the same length!"))
        }

        ## Length converted catch curve
        if("midLengths" %in% names(res)){
            if("par" %in% names(res)){
                  Linf <- res$par$Linf
                  K <- res$par$K
                  t0 <- ifelse("t0" %in% names(res$par), res$par$t0, 0)
                  C <- ifelse("C" %in% names(res$par), res$par$C, 0)
                  ts <- ifelse("ts" %in% names(res$par), res$par$ts, 0)            
              }else{
                  Linf <- res$Linf
                  K <- res$K
                  t0 <- ifelse("t0" %in% names(res), res$t0, 0)
                  C <- ifelse("C" %in% names(res), res$C, 0)
                  ts <- ifelse("ts" %in% names(res), res$ts, 0)
              }
            

          if((is.null(Linf) | is.null(K))) stop(noquote(
            "You need to assign values to Linf and K for the catch curve based on length-frequency data!"))

          #calculate size class interval
          midLengths <- classes.num
          interval <- midLengths[2] - midLengths[1]

            ## L and t of lower length classes

            lowerLengths <- midLengths - (interval / 2)
            t_L1 <- VBGF(param = list(Linf = Linf, K = K, t0 = t0, C=C, ts=ts), L = lowerLengths)
          ## t0 - (1/K) * log(1 - (lowerLengths / Linf))

          # delta t
          dt <- rep(NA,length(midLengths))
          for(x1 in 1:(length(dt)-1)){
            dt[x1] <- t_L1[x1+1] - t_L1[x1]
          }

          # x varaible
          #ln (Linf - L)
          ln_Linf_L <- log(Linf - lowerLengths)
          ## t of midlengths
          if("C" %in% names(res) & "ts" %in% names(res)){
              t_midL <- VBGF(param = list(Linf = Linf, K = K, t0 = t0, C=C, ts=ts), L = midLengths)
          }else{
              t_midL <- VBGF(param = list(Linf = Linf, K = K, t0 = t0), L = midLengths)
          }
          ## t0 - (1/K) * log(1 - (midLengths / Linf))

          # y variable
          #ln C(L1,Linf)
          lnC <- log(catch)
          # ln( Catch / delta t)
          lnC_dt <- log(catch / dt)
          lnC_dt[which(lnC_dt == -Inf)] <- NA   ### OR zero???

          if(cumulative == FALSE){
            xvar = t_midL
            yvar = lnC_dt
            xname = "t_midL"
            yname = "lnC_dt"
            xlabel = "Relative age [yrs]"
            ylabel = "ln(C/dt)"
          }

          if(cumulative){
            xvar = ln_Linf_L
            yvar = lnC
            xname = "ln_Linf_L"
            yname = "lnC"
            xlabel = "ln (Linf - L)"
            ylabel = "ln C(L, Linf)"
          }
        }

        # Aged based Catch curve
        if("age" %in% names(res) == TRUE){
          # delta t
          if(constant_dt == FALSE){
            dt <- rep(NA,length(classes.num))
            for(x1 in 1:(length(dt)-1)){
              dt[x1] <- classes.num[x1+1] - classes.num[x1]
            }
          }
          if(constant_dt) dt <- rep(classes.num[2] - classes.num[1], length(classes.num))

          # x variable
          # (t + dt) / 2   ==   x
          if(cumulative == FALSE) tplusdt_2 <- classes.num + (dt / 2)
          if(cumulative) tplusdt_2 <- classes.num

          # y variable
          #ln C(L1,Linf)
          lnC <- log(catch)
          # ln( Catch / delta t)     ==    y
          lnC_dt <- log(catch / dt)


          if(constant_dt){
            xvar = classes.num
            xname = "classes.num"
            xlabel = "Age [yrs]"
            yvar = lnC
            yname = "lnC"
            ylabel = "ln C(t, inf)"
          }

          if(cumulative == FALSE & constant_dt == FALSE){
            xvar = tplusdt_2
            xname = "tplusdt_2"
            xlabel = "Age [yrs]"
            yvar = lnC_dt
            yname = "lnC_dt"
            ylabel = "ln(C/dt)"
          }

          if(cumulative){
            xvar = tplusdt_2
            xname = "tplusdt_2"
            xlabel = "Age [yrs]"
            yvar = lnC
            yname = "lnC"
            ylabel = "ln C(t, inf)"
          }
        }



        ## remove all NAs and Infs
        temp <- cbind(xvar,yvar)
        temp <- as.matrix(na.exclude(temp))
        temp <- temp[(!(temp[,1] == Inf | temp[,1] == -Inf)),]
        temp <- temp[(!(temp[,2] == Inf | temp[,2] == -Inf)),]
        xvar <- temp[,1]
        yvar <- temp[,2]

        
        #for plot
        #minY <- ifelse(min(yvar,na.rm=TRUE) < 0, min(yvar,na.rm=TRUE),0)
        minY <- min(yvar,na.rm=TRUE)
        maxY <- max(yvar,na.rm=TRUE) + 1
        xlims <- c(0, max(xvar,na.rm=TRUE))

          cutterList <- vector("list", reg_num)
        #identify plot
        if(is.null(reg_int) & !auto){
          writeLines("Please choose the minimum and maximum point in the graph \nto include for the regression line!")
          flush.console()
          for(I in 1:reg_num){

              dev.new()#noRStudioGD = TRUE)
              op <- par(mfrow = c(1,1),
                        c(5, 4, 4, 2) + 0.1,
                        oma = c(2, 1, 0, 1) + 0.1)
              plot(x = xvar,y = yvar, ylim = c(minY,maxY), xlim = xlims,
                   xlab = xlabel, ylab = ylabel, type = "n")

              ## plot previous regression lines when using multiple regression lines
              if(I > 1){
                  for(II in 1:(I-1)){
                      points(xvar[cutterList[[II]]],yvar[cutterList[[II]]], col="darkgreen", pch=16)
                      lines(xvar[cutterList[[II]]],yvar[cutterList[[II]]], col="darkgreen", lwd=2)
                  }
              }
              mtext(side = 3, "Click on two numbers. Escape to Quit.",
                    xpd = NA, cex = 1.25)
              text(xvar, yvar, labels=as.character(order(xvar)), cex= 0.7)
              cutter <- identify(x = xvar, y = yvar,
                                 labels = order(xvar), n=2, col = "red")
              par(op)

              if(is.na(cutter[1]) | is.nan(cutter[1]) |
                 is.na(cutter[2]) | is.nan(cutter[2]) ) stop(noquote("You did not choose any points in the graph. Please re-run the function and choose points in the graph!"))

              dev.off()

              ## save results to list
              cutterList[[I]] <- cutter
          }


        }


          if(!is.null(reg_int)){
              cutterList <- reg_int
              if(class(cutterList) != "list" && length(cutterList) != 2) stop("You have to provide 2 numbers in reg_int.")
              if(class(cutterList) == "list" && any(unlist(lapply(cutterList,length)) != 2)) stop("You have to provide 2 numbers in reg_int.")
          }

          if(auto){
              yvar2 <- as.numeric(yvar)
              xvar2 <- xvar[which(yvar2 > 0.2)]
              cutter <- c(which(yvar2 == max(as.numeric(yvar2),na.rm=TRUE))+1, which(xvar2 == max(xvar2,na.rm=TRUE)))
              cutterList <- list()
              cutterList[[1]] <- cutter
          }

          ## define result lists
          lm1List <- vector("list",reg_num)
          Z_lm1List <- vector("list",reg_num)
          SE_Z_lm1List <- vector("list",reg_num)
          conf_Z_lm1List <- vector("list",reg_num)
          intercept_lm1List <- vector("list",reg_num)


          for(I in 1:reg_num){

              if(class(cutterList) == "list"){
                  cutter <- cutterList[[I]]
              }else{
                  cutter <- cutterList
              }

              ## calculations + model
              df.CC <- as.data.frame(cbind(xvar,yvar))
              df.CC.cut <- df.CC[cutter[1]:cutter[2],]
              lm1 <- lm(yvar ~ xvar, data = df.CC.cut)
              sum_lm1 <- summary(lm1)
              r_lm1 <- sum_lm1$r.squared
              intercept_lm1 <- sum_lm1$coefficients[1]
              slope_lm1 <- sum_lm1$coefficients[2]
              se_slope_lm1 <- sum_lm1$coefficients[4]

              ## fit of regression line
              lm1.fit <- sum_lm1$r.squared
              Z_lm1 <- abs(slope_lm1)
              SE_Z_lm1 <- abs(se_slope_lm1)
              confi <-  abs(se_slope_lm1) * qt(0.975,sum_lm1$df[2])
              conf_Z_lm1 <- Z_lm1 + c(-confi,confi)

              ## special case when cumulative and length-frequency data
              if(cumulative & "midLengths" %in% names(res) == TRUE){
                Z_lm1 <- Z_lm1 * K
                SE_Z_lm1 <- SE_Z_lm1 * K
              }

              ## save results to lists
              lm1List[[I]] <- lm1
              Z_lm1List[[I]] <- Z_lm1
              SE_Z_lm1List[[I]] <- SE_Z_lm1
              conf_Z_lm1List[[I]] <- conf_Z_lm1
              intercept_lm1List[[I]] <- intercept_lm1
          }

          ##save all in list
          if(reg_num > 1){
              ret <- c(res,list(
                               xvar = xvar,
                               yvar = yvar,
                               reg_int = cutterList,
                               linear_mod = lm1List,
                               Z =  Z_lm1List,
                               se = SE_Z_lm1List,
                               confidenceInt = conf_Z_lm1List))        
          }else{
              ret <- c(res,list(
                               xvar = xvar,
                               yvar = yvar,
                               reg_int = unlist(cutterList),
                               linear_mod = lm1List[[1]],
                               Z = unlist(Z_lm1List),
                               se = unlist(SE_Z_lm1List),
                               confidenceInt = unlist(conf_Z_lm1List)))
          }

          if("M" %in% names(ret) && length(ret$M)==1){
              ret$FM <- lapply(ret$Z, function(x) x - ret$M)
          }
          names(ret)[names(ret) == "xvar"] <- xname
          names(ret)[names(ret) == "yvar"] <- yname
          class(ret) <- "catchCurve"

        # Calculate selection ogive from catch curve and add to ret
        if(calc_ogive & cumulative) stop(noquote("It is not possible to estimate the selection ogive for the cumulative catch curve."))
          if(calc_ogive){

              ## Assumption that Z of smallest selected individuals is most appropriate
              mini <- min(unlist(cutterList))
              temp <- lapply(cutterList, function(x) grep(mini,x))
              ind <- sapply(temp, function(x) length(x) > 0)
              cutter <- unlist(cutterList[ind])


              ## only use part of catch and t which is not fully exploited by the gear
              t_ogive <- xvar[1:(cutter[1]-1)]
              dt_ogive <- dt[1:(cutter[1]-1)]
              if("age" %in% names(res) == TRUE &
                 class(catch) == 'matrix' | class(catch) == 'data.frame'){
                catch_ogive <- catch[1:(cutter[1]-1)] ## catch.cohort
              }else catch_ogive <- catch[1:(cutter[1]-1)]


              # calculate observed selection ogive
              Sobs <- catch_ogive/(dt_ogive * exp(unlist(intercept_lm1List[ind]) - unlist(Z_lm1List[ind]) * t_ogive))

              # dependent vairable in following regression analysis
              ln_1_S_1 <- log((1/Sobs) - 1)

              # get rid of Inf
              ln_1_S_1[which(ln_1_S_1 == Inf)] <- NA
              t_ogive[which(t_ogive == Inf)] <- NA

              #regression analysis to caluclate T1 and T2
              mod_ogive <- lm(ln_1_S_1 ~ t_ogive, na.action = na.omit)
              sum_lm_ogive <- summary(mod_ogive)
              T1 <- sum_lm_ogive$coefficients[1]
              T2 <- abs(sum_lm_ogive$coefficients[2])

              # calculate estimated selection ogive
              Sest <- 1/(1+exp(T1 - T2*xvar))

              # selection parameters
              t50 <- T1/T2
              t75 <- (T1 + log(3))/T2
              t95 <-  (T1 - log((1 / 0.95) - 1)) / T2
              if((!is.null(res$Linf) & !is.null(res$K)) | ("par" %in% names(res) &&
                 (!is.null(res$par$Linf) & !is.null(res$par$K)))){
                if(is.null(res$t0)) t0 = 0
                L50 <- Linf*(1-exp(-K*(t50-t0)))
                L75 <- Linf*(1-exp(-K*(t75-t0)))
                L95 <- Linf*(1-exp(-K*(t95-t0)))
              }

              ret2 <- c(ret,list(
                intercept = intercept_lm1,
                linear_mod_sel = mod_ogive,
                Sobs = Sobs,
                ln_1_S_1 = ln_1_S_1,
                Sest = Sest,
                t50 = t50,
                t75 = t75,
                t95 = t95))
              if(exists("L50")) ret2$L50 = L50
              if(exists("L75")) ret2$L75 = L75
              if(exists("L95")) ret2$L95 = L95
              if(exists("L50")) names(ret2)[which(ret2 %in% L50)] <- "L50"
              if(exists("L75")) names(ret2)[which(ret2 %in% L75)] <- "L75"
              if(exists("L95")) names(ret2)[which(ret2 %in% L95)] <- "L95"

              class(ret2) <- "catchCurve"
              if(plot) plot(ret2, plot_selec=TRUE)
              return(ret2)
            }else {
              if(plot) plot(ret)
              return(ret)
            }
    }
}
