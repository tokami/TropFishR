#' @title Empirical formulas for the estimation of natural mortality
#
#' @description Functions to calculate the instantaneous natural mortality rate (M)
#'      according to 12 different empirical formulas.
#'
#' @param Linf infinite total length (TL) from a von Bertalanffy
#'    growth curve in cm.
#' @param Winf infinite weight form a von Bertalanffy growth curve
#'    in wet weight-grams.
#' @param K_l is the growth coefficient (per year) from a von Bertalanffy growth
#'    curve for length.
#' @param K_w is the growth coefficient (per year) from a von Bertalanffy growth
#'    curve for weight.
#' @param temp average annual temperature at the surface in degrees centigrade.
#' @param tmax the oldest age observed for the species.
#' @param tm50 age when 50\% of the population is mature [year]
#'      ("age of massive maturation").
#' @param GSI gonadosomatic index (wet ovary weight over wet body weight).
#' @param Wdry total dry weight in grams.
#' @param Wwet total wet weight at mean length in grams.
#' @param Bl vector with body lengths in cm for size dependent mortality estimates (method = "Gislason")
#' @param schooling logical; if TRUE it is accounted for the schooling behaviour of
#'      the species, only for Pauly's methods. Default is FALSE.
#' @param method vector of method names. Any combination of following methods can
#'    be employed: "AlversonCarney", "Gislason" (size dependent mortality estimates), "GundersonDygert", "Hoenig",
#'    "Lorenzen", "Pauly_Linf", "Pauly_Winf", "PetersonWroblewski",
#'    "RikhterEfanov", "Roff", "Then_growth", or "Then_tmax".
#'    Please refer to Details to see which input parameters
#'    are required by each method.
#' @param boot an object of class 'lfqBoot'
#' @param CI percentage for confidence intervals (Default: 95)
#'
#' @keywords function mortality M
#'
#' @examples
#' M_empirical(Linf = 80, K_l = 0.5, temp = 25, tmax = 30,
#'      method = c("Pauly_Linf","Hoenig"))
#'
#' ## bootstrapping application of M
#' ## coming soon
#'
#' @source https://cran.r-project.org/web/packages/fishmethods/index.html
#'
#' @details Function adapted from the mortality function of the fishmethods package
#'     by Gary A. Nelson
#'     (https://cran.r-project.org/web/packages/fishmethods/index.html). The function allows to
#'     be applied to a the results of the bootstrapping ELEFAN functions. This can be done by providing
#'     the 'lfqBoot' object by use of the argument \code{boot}. So far only the methods 'Pauly_Linf' and
#'     'Then_growth' can be used.
#'
#'
#' Depending on the method different input parameters are required:
#' \itemize{
#'    \item \code{"AlversonCarney"} requires \code{K_l} and \code{tmax},
#'    \item \code{"Gislason"} requires \code{Linf}, \code{K_l} and \code{Bl},
#'    \item \code{"GundersonDygert"} requires \code{GSI},
#'    \item \code{"Hoenig"} requires \code{tmax},
#'    \item \code{"Lorenzen"} requires \code{Wwet},
#'    \item \code{"Pauly_Linf"} requires \code{Linf}, \code{K_l} and \code{temp},
#'    \item \code{"Pauly_Winf"} requires \code{Winf}, \code{K_w} and \code{temp},
#'    \item \code{"PetersonWroblewski"} requires \code{Wdry},
#'    \item \code{"RikhterEfanov"} requires \code{tm50},
#'    \item \code{"Roff"} requires \code{K_l} and \code{tm50},
#'    \item \code{"Then_tmax"} requires \code{tmax},
#'    \item \code{"Then_growth"} requires \code{Linf} and \code{K_l}.
#' }
#' If accounting for schooling behaviour M is multiplied by 0.8 according to
#'    Pauly (1983).
#'
#' @importFrom ks kde
#'
#' @return A matrix of M estimates.
#'
#' @references
#' Alverson, D. L. and M. J. Carney. 1975. A graphic review of the growth and decay
#' of population cohorts. J. Cons. Int. Explor. Mer 36: 133-143.
#'
#' Gislason, H., N. Daan, J. C. Rice, and J. G. Pope. 2010. Size, growth,
#' temperature and the natural mortality of marine fish. Fish and Fisheries 11: 149-158.
#'
#' Gunderson, D. R. and P. H. Dygert. 1988. Reproductive effort as a predictor
#' of natural mortality rate. J. Cons. Int. Explor. Mer 44: 200-209.
#'
#' Hoenig, J. M. 1983. Empirical use of longevity data to estimate mortality rates.
#' Fish. Bull. 82: 898-903.
#'
#' Lorenzen, K. 1996. The relationship between body weight and natural mortality in
#' juvenile and adult fish: a comparison of natural ecosystems and aquaculture.
#' J. Fish. Biol. 49: 627-647.
#'
#' Pauly, D. 1980. On the interrelationships between natural mortality,
#' growth parameters, and mean environmental temperature in 175 fish stocks.
#' J. Cons. Int. Explor. Mer: 175-192.
#'
#' Pauly, D., 1983. Some simple methods for the assessment of tropical fish stocks.
#' \emph{FAO Fish.Tech.Pap.}, (234): 52p. Issued also in French and Spanish
#'
#' Peterson, I. and J. S. Wroblewski. 1984. Mortality rate of fishes in the
#' pelagic ecosystem. Can. J. Fish. Aquat. Sci. 41: 1117-1120.
#'
#'pp Rikhter, V.A., and V.N. Efanov, 1976. On one of the approaches to estimation of natural
#' mortality of fish populations. \emph{ICNAF Res.Doc.}, 76/VI/8: 12p.
#'
#' Roff, D. A. 1984. The evolution of life history parameters in teleosts.
#' Can. J. Fish. Aquat. Sci. 41: 989-1000.
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' Then, A. Y., J. M. Hoenig, N. G. Hall, D. A. Hewitt. 2015. Evaluating the predictive
#' performance of empirical estimators of natural mortality rate using information on over
#' 200 fish species. ICES J. Mar. Sci. 72: 82-92.
#'
#' @export

M_empirical <- function(Linf = NULL, Winf = NULL, K_l = NULL, K_w = NULL,
                        temp = NULL, tmax = NULL, tm50 = NULL, GSI = NULL,
                        Wdry = NULL, Wwet = NULL, Bl = NULL,
                        schooling = FALSE, method, boot = NULL, CI = 95){

    if(!is.null(boot) & class(boot) == "lfqBoot"){
        ## calculations with bootstrapped ELEFAN results
        bootOut <- boot
        bootRaw <- boot$bootRaw

        if (any(method == "Pauly_Linf") & any(is.null(bootRaw$Linf), is.null(bootRaw$K), is.null(temp)))
          stop("Pauly_Linf requires temp and a boot object with columns Linf and K")
        if (any(method == "Then_growth") & any(is.null(bootRaw$Linf), is.null(bootRaw$K)))
          stop("Then_growth requires a boot object with columns Linf and K")


        if(any(method == "Pauly_Linf")){
            Ms <- round(10^(-0.0066 - 0.279 * log10(bootRaw$Linf) +
                            0.6543 * log10(bootRaw$K) + 0.4634 * log10(temp)), 3)
            if(schooling == TRUE){
                Ms <- 0.8 * Ms
            }
            bootRaw[,ncol(bootRaw)+1]  <- Ms
            colnames(bootRaw) <- c(colnames(bootRaw)[-ncol(bootRaw)], "M_Pauly")
        }
        if (any(method == "Then_growth")) {
            bootRaw[,ncol(bootRaw)+1] <- round(4.118 * (bootRaw$K^0.73) * (bootRaw$Linf^-0.33), 3)
            colnames(bootRaw) <- c(colnames(bootRaw)[-ncol(bootRaw)], "M_Then")
        }

        tmp <- as.data.frame(bootRaw[,(ncol(boot$bootRaw)+1:length(method))])

        ## max density and CIS
        resMaxDen <- vector("numeric", ncol(tmp))
        resMed <- vector("numeric", ncol(tmp))                
        ciList <- vector("list", ncol(tmp))
        for(i in seq(length(method))){
            ## median
            resMed[i] <- median(tmp[,i], na.rm = TRUE)
            
            ## confidence intervals
            citmp <- (100-CI)/2/100
            ciList[[i]] <- quantile(tmp[,i],  probs = c(citmp, 1-citmp), na.rm = TRUE)
            
            ## max densities
            x <- try(ks::kde(as.numeric(na.omit(tmp[,i]))), TRUE)
            if(class(x) != "try-error"){
                ## max den
                ind <- which(x$estimate > x$cont["99%"])
                resMaxDen[i] <- mean(x$eval.points[ind])  ## why mean?
            }else{
                if((length(unique(as.character(tmp[,i]))) == 1 && all(!is.na(tmp[,i]))) |
                   (length(unique(as.character(na.omit(tmp[,i])))) == 1)){
                    resMaxDen[i] <- unique(as.numeric(na.omit(tmp[,i])))
                }else{
                    resMaxDen[i] <- NA
                }
            }
        }
        resCIs <- cbind(boot$CI,t(do.call(rbind,ciList)))
        colnames(resCIs) <- colnames(bootRaw)
        rownames(resCIs) <- c("lo","up")
        resMaxDen <- c(boot$maxDen, resMaxDen)
        names(resMaxDen) <- colnames(bootRaw)
        resMed <- c(boot$median, resMed)
        names(resMed) <- colnames(bootRaw)


        ret <- list()
        ret$bootRaw <- bootRaw
        ret$seed <- boot$seed
        ret$maxDen <- resMaxDen
        ret$median <- resMed        
        ret$CI <- resCIs
        if("multiCI" %in% names(boot)) ret$multiCI <- boot$multiCI
        class(ret) <- "lfqBoot"
        return(ret)

    }else if(!is.null(boot) & class(boot) != "lfqBoot"){
        stop("You provided an object for boot, but it does not have class 'lfqBoot'. Please check.")
    }else{

        if (any(method == "AlversonCarney") & any(is.null(tmax), is.null(K_l)))
          stop("AlversonCarney requires K_l and tmax")
        if (any(method == "Gislason") & any(is.null(Linf), is.null(K_l), is.null(Bl)))
          stop("Gislason requires Linf, K_l, and Bl")
        if (any(method == "GundersonDygert") & is.null(GSI))
          stop("GundersonDygert requires GSI")
        if (any(method == "Hoenig") & is.null(tmax))
          stop("Hoenig requires tmax")
        if (any(method == "Lorenzen") & is.null(Wwet))
          stop("Lorenzen requires Wwet")
        if (any(method == "Pauly_Linf") & any(is.null(Linf), is.null(K_l), is.null(temp)))
          stop("Pauly_Linf requires Linf, K_l, and temp")
        if (any(method == "Pauly_Winf") & any(is.null(Winf), is.null(K_w), is.null(temp)))
          stop("Pauly_Winf requires Winf, K_w, and temp")
        if (any(method == "PetersonWroblewski") & is.null(Wdry))
          stop("PetersonWroblewski requires Wdry")
        if (any(method == "RikhterEfanov") & any(is.null(tm50)))
          stop("RikhterEfanov requires K_l and tm50")
        if (any(method == "Roff") & any(is.null(tm50), is.null(K_l)))
          stop("Roff requires K_l and tm50")
        if (any(method == "Then_tmax") & any(is.null(tmax)))
          stop("Then_max requires tmax")
        if (any(method == "Then_growth") & any(is.null(Linf), is.null(K_l)))
          stop("Then_growth requires Linf and K_l")

        n <- length(method)
        if (any(method == "Hoenig"))
          n <- n + 1
        M_mat <- matrix(NA, n, 1L)
        dimnames(M_mat) <- list(rep(NA, n), c("M"))
        ind <- 0

        if(any(method == "AlversonCarney")){
          ind <- ind + 1
          # Alverson and Carney (1975)
          M_mat[ind, 1]  <- round((3 * K_l)/(exp(K_l * (0.38 * tmax)) - 1), 3)
          dimnames(M_mat)[[1]][ind] <- list("Alverson and Carney (1975)")
        }
        #if(any(method == "Gislason")){
        #  ind <- ind + 1
        #  # Gislason et al. (2010)
        #  M_mat[ind, 1]  <- round(exp(0.55 - 1.61 * log(Bl) + 1.44 * log(Linf) + log(K_l)), 3)
        #  dimnames(M_mat)[[1]][ind] <- list("Gislason et al. (2010)")
        #}
        if(any(method == "GundersonDygert")){
          ind <- ind + 1
          # Gunderson and Dygert (1988)
          M <- round(0.03 + 1.68 * GSI, 3)
          dimnames(M_mat)[[1]][ind] <- list("Gunderson and Dygert (1988)")
        }
        if(any(method == "Hoenig")){
          ind <- ind + 1
          # Hoenig (1983) - Joint Equation
          M_mat[ind, 1]  <- round(4.22/(tmax^0.982), 3)
          dimnames(M_mat)[[1]][ind] <- list("Hoenig (1983) - Joint Equation")

          ind <- ind + 1
          # Hoenig (1983) - Fish Equation
          M_mat[ind, 1]  <- round(exp(1.46 - 1.01 * log(tmax)), 3)
          dimnames(M_mat)[[1]][ind] <- list("Hoenig (1983) - Fish Equation")
        }
        if(any(method == "Lorenzen")){
          ind <- ind + 1
          # Lorenzen (1996)
          M_mat[ind, 1]  <- round(3 * (Wwet^-0.288), 3)
          dimnames(M_mat)[[1]][ind] <- list("Lorenzen (1996)")
        }
        if(any(method == "Pauly_Linf")){
          ind <- ind + 1
          M_mat[ind, 1]  <- round(10^(-0.0066 - 0.279 * log10(Linf) + 0.6543 * log10(K_l) + 0.4634 * log10(temp)), 3)  #exp( -0.0152 - 0.279 * log(Linf) + 0.6543 * log(K) + 0.463 * log(temp))
          dimnames(M_mat)[[1]][ind] <- list("Pauly (1980) - Length Equation")
          if(schooling == TRUE){
            M_mat[ind, 1] <- 0.8 * M_mat[ind, 1]
          }
        }
        if(any(method == "Pauly_Winf")){
          ind <- ind + 1
          M_mat[ind, 1]  <- round(10^(-0.2107 - 0.0824 * log10(Winf) + 0.6757 * log10(K_w) + 0.4627 * log10(temp)), 3)  #exp( -0.2107 - 0.0824 * log(Winf) + 0.6757 * log(K) + 0.4627 * log(temp))
          dimnames(M_mat)[[1]][ind] <- list("Pauly (1980) - Weight Equation")
          if(schooling == TRUE){
            M_mat[ind, 1] <- 0.8 * M_mat[ind, 1]
          }
        }
        if(any(method == "PetersonWroblewski")){
          ind <- ind + 1
          # Peterson and Wroblewski (1984)
          M_mat[ind, 1]  <- round(1.92 * (Wdry^-0.25), 3)
          dimnames(M_mat)[[1]][ind] <- list("Peterson and Wroblewski (1984)")
        }
        if(any(method == "RikhterEfanov")){
          ind <- ind + 1
          M_mat[ind, 1]  <- round(1.521 / ( tm50 ^ 0.720) - 0.155, 3)
          dimnames(M_mat)[[1]][ind] <- list("Rikhter and Efanov (1976)")
        }
        if(any(method == "Roff")){
          ind <- ind + 1
          # Roff (1984)
          M_mat[ind, 1]  <- round((3 * K_l)/(exp(K_l * tm50) - 1), 3)
          dimnames(M_mat)[[1]][ind] <- list("Roff (1984)")
        }
        if (any(method == "Then_tmax")) {
          ind <- ind + 1
          M_mat[ind, 1]  <- round(4.899 * tmax^-0.916, 3)
          dimnames(M_mat)[[1]][ind] <- list("Then (2015) - tmax")
        }
        if (any(method == "Then_growth")) {
          ind <- ind + 1
          M_mat[ind, 1]  <- round(4.118 * (K_l^0.73) * (Linf^-0.33), 3)
          dimnames(M_mat)[[1]][ind] <- list("Then (2015) - growth")
        }
        if (any(method == "Gislason")) {
          Ml <- round(exp(0.55 - 1.61 * log(Bl) + 1.44 *
                                       log(Linf) + log(K_l)), 3)
          M_mat <- as.data.frame(matrix(c(Bl,Ml),byrow = FALSE,ncol=2))
          colnames(M_mat) <- c("Bl","Ml")
        }
        return(M_mat)

    }

}






#' @title Empirical M formula by Then in bootstrap framework
#
#' @description Resamples from Then's original data and refits
#'     original model within bootstrap framework to estimate M with confidence intervals.
#' 
#' @param boot an object of class 'lfqBoot'
#' @param CI percentage for confidence intervals (Default: 95)
#'
#' @importFrom ks kde
#'
#' @references
#' Then, A. Y., J. M. Hoenig, N. G. Hall, D. A. Hewitt. 2015. Evaluating the predictive
#' performance of empirical estimators of natural mortality rate using information on over
#' 200 fish species. ICES J. Mar. Sci. 72: 82-92.
#'
#' @export

MstochThen <- function(boot, CI=95){

    Mboot <- function(preddat, datThenNA, modThen){
        sampi <- sample(x =  seq(nrow(datThenNA)), size = nrow(datThenNA), replace = TRUE)
        dfi <- datThenNA[sampi,]
        modTheni <- update(modThen, data = dfi)
        predi <- predict(modTheni, newdata = preddat)
        return(predi)
    }

    bootRaw <- boot$bootRaw

    ## load Then's data
    data("datThen")

    ## original model by Then
    modThen <- nls(M ~ a * K^b * Linf^c, data= datThen, start=list(a=4,b=0.8,c=-0.3))

    ## remove NA for resampling
    datThenNA <- na.omit(data.frame(M = datThen$M,K = datThen$K,Linf = datThen$Linf))

    Ms <- apply(bootRaw, 1, Mboot, datThenNA = datThenNA, modThen = modThen)

    bootRaw[,ncol(bootRaw)+1] <- Ms
    colnames(bootRaw) <- c(colnames(bootRaw)[-ncol(bootRaw)], "M_Then")
    

    ## max density and CIS
    ## median
    resMed <- median(Ms, na.rm = TRUE)

    ## confidence intervals
    citmp <- (100-CI)/2/100
    ciList <- quantile(Ms,  probs = c(citmp, 1-citmp), na.rm = TRUE)

    ## max densities
    x <- try(ks::kde(as.numeric(na.omit(Ms))), TRUE)
    if(class(x) != "try-error"){
        ## max den
        ind <- which(x$estimate > x$cont["99%"])
        resMaxDen <- mean(x$eval.points[ind])  ## why mean?
    }else{
        if((length(unique(as.character(Ms))) == 1 && all(!is.na(Ms))) |
           (length(unique(as.character(na.omit(Ms)))) == 1)){
            resMaxDen <- unique(as.numeric(na.omit(Ms)))
        }else{
            resMaxDen <- NA
        }
    }
    resCIs <- cbind(boot$CI,as.matrix(ciList))
    colnames(resCIs) <- colnames(bootRaw)
    rownames(resCIs) <- c("lo","up")
    resMaxDen <- c(boot$maxDen, resMaxDen)
    names(resMaxDen) <- colnames(bootRaw)
    resMed <- c(boot$median, resMed)
    names(resMed) <- colnames(bootRaw)
    

    ret <- list()
    ret$bootRaw <- bootRaw
    ret$seed <- boot$seed
    ret$maxDen <- resMaxDen
    ret$median <- resMed        
    ret$CI <- resCIs
    if("multiCI" %in% names(boot)) ret$multiCI <- boot$multiCI
    class(ret) <- "lfqBoot"
    return(ret)
}




