#' @title Empirical formulas for the estimation of natural mortality
#'
#' @description Functions to calculate the instantaneous natural mortality rate (M)
#'      according to 12 different empirical formulas.
#' @param lfq a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths}: midpoints of the length classes (length-frequency
#'   data; in cm) for size dependent mortality estimates (method = "Gislason")
#' \item \strong{par}: and following parameters (not all are required; dependent on method):
#'  \itemize{
#' \item \strong{Linf} infinite total length (TL) from a von Bertalanffy
#'    growth curve in cm.
#' \item \strong{Winf} infinite weight form a von Bertalanffy growth curve
#'    in wet weight-grams.
#' \item \strong{K} is the growth coefficient (per year) from a von Bertalanffy growth
#'    curve for length.
#' \item \strong{K_w} is the growth coefficient (per year) from a von Bertalanffy growth
#'    curve for weight.
#' \item \strong{temp} average annual temperature at the surface in degrees centigrade.
#' \item \strong{tmax} the oldest age observed for the species.
#' \item \strong{tm50} age when 50\% of the population is mature [year]
#'      ("age of massive maturation").
#' \item \strong{GSI} gonadosomatic index (wet ovary weight over wet body weight).
#' \item \strong{Wdry} total dry weight in grams.
#' \item \strong{Wwet} total wet weight at mean length in grams.
#'  }
#' }
#' @param schooling logical; if TRUE it is accounted for the schooling behaviour of
#'      the species, only for Pauly's methods. Default is FALSE.
#' @param method vector of method names. Any combination of following methods can
#'    be employed: "AlversonCarney", "Gislason" (size dependent mortality estimates), "GundersonDygert", "Hoenig",
#'    "Lorenzen", "Pauly_Linf", "Pauly_Winf", "PetersonWroblewski",
#'    "RikhterEfanov", "Roff", "Then_growth", or "Then_tmax".
#'    Please refer to Details to see which input parameters
#'    are required by each method.
#'
#' @keywords function mortality M
#'
#' @examples
#' ## with included data set
#' res <- M_empirical(synLFQ3, method="Then_growth")
#'
#' ## with specific parameters
#' pars <- list(par=list(Linf = 80, K = 0.5, temp = 25, tmax = 30))
#' M_empirical(pars, method = c("Pauly_Linf","Hoenig"))
#'
#' @source https://cran.r-project.org/web/packages/fishmethods/index.html
#'
#' @details Function adapted from the mortality function of the fishmethods package
#'     by Gary A. Nelson
#'     (https://cran.r-project.org/web/packages/fishmethods/index.html).
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
#' Rikhter, V.A., and V.N. Efanov, 1976. On one of the approaches to estimation of natural
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

M_empirical <- function(lfq, method, schooling = FALSE){

    res <- lfq
    if(!"par" %in% names(lfq)) stop(noquote("Please provide the required parameters in res$par!"))
    par <- lfq$par

    Linf <- ifelse("Linf" %in% names(par), par$Linf, NA)
    K_l <- ifelse("K" %in% names(par), par$K, NA)            ## old K_l
    K_w <- ifelse("K_w" %in% names(par), par$K_w, NA)
    Winf <- ifelse("Winf" %in% names(par), par$Winf, NA)
    temp <- ifelse("temp" %in% names(par), par$temp, NA)
    tmax <- ifelse("tmax" %in% names(par), par$tmax, NA)
    tm50 <- ifelse("tm50" %in% names(par), par$tm50, NA)
    GSI <- ifelse("GSI" %in% names(par), par$GSI, NA)
    Wdry <- ifelse("Wdry" %in% names(par), par$Wdry, NA)
    Wwet <- ifelse("Wwet" %in% names(par), par$Wwet, NA)
    Bl <- ifelse("midLengths" %in% names(res), res$midLengths, NA)


    if (any(method == "AlversonCarney") & any(is.na(tmax), is.na(K_l)))
        stop("AlversonCarney requires K and tmax in lfq$par")
    if (any(method == "Gislason") & any(is.na(Linf), is.na(K_l), is.na(Bl)))
        stop("Gislason requires Linf, K in lfq$par and midLengths in lfq")
    if (any(method == "GundersonDygert") & is.na(GSI))
        stop("GundersonDygert requires GSI in lfq$par")
    if (any(method == "Hoenig") & is.na(tmax))
        stop("Hoenig requires tmax in lfq$par")
    if (any(method == "Lorenzen") & is.na(Wwet))
        stop("Lorenzen requires Wwet in lfq$par")
    if (any(method == "Pauly_Linf") & any(is.na(Linf), is.na(K_l), is.na(temp)))
        stop("Pauly_Linf requires Linf, K, and temp in lfq$par")
    if (any(method == "Pauly_Winf") & any(is.na(Winf), is.na(K_w), is.na(temp)))
        stop("Pauly_Winf requires Winf, K_w, and temp in lfq$par")
    if (any(method == "PetersonWroblewski") & is.na(Wdry))
        stop("PetersonWroblewski requires Wdry in lfq$par")
    if (any(method == "RikhterEfanov") & any(is.na(tm50)))
        stop("RikhterEfanov requires K and tm50 in lfq$par")
    if (any(method == "Roff") & any(is.na(tm50), is.na(K_l)))
        stop("Roff requires K and tm50 in lfq$par")
    if (any(method == "Then_tmax") & any(is.na(tmax)))
        stop("Then_max requires tmax in lfq$par")
    if (any(method == "Then_growth") & any(is.na(Linf), is.na(K_l)))
        stop("Then_growth requires Linf and K in lfq$par")

    n <- length(method)
    if (any(method == "Hoenig"))
        n <- n + 1
    M_mat <- matrix(NA, n, 1L)
    dimnames(M_mat) <- list(rep(NA, n), c("M"))
    ind <- 0

    if(any(method == "AlversonCarney")){
        ind <- ind + 1
        ## Alverson and Carney (1975)
        M_mat[ind, 1]  <- round((3 * K_l)/(exp(K_l * (0.38 * tmax)) - 1), 3)
        dimnames(M_mat)[[1]][ind] <- list("Alverson and Carney (1975)")
    }
    ## if(any(method == "Gislason")){
    ##  ind <- ind + 1
    ##  ## Gislason et al. (2010)
    ##  M_mat[ind, 1]  <- round(exp(0.55 - 1.61 * log(Bl) + 1.44 * log(Linf) + log(K_l)), 3)
    ##  dimnames(M_mat)[[1]][ind] <- list("Gislason et al. (2010)")
    ## }
    if(any(method == "GundersonDygert")){
        ind <- ind + 1
        ## Gunderson and Dygert (1988)
        M <- round(0.03 + 1.68 * GSI, 3)
        dimnames(M_mat)[[1]][ind] <- list("Gunderson and Dygert (1988)")
    }
    if(any(method == "Hoenig")){
        ind <- ind + 1
        ## Hoenig (1983) - Joint Equation
        M_mat[ind, 1]  <- round(4.22/(tmax^0.982), 3)
        dimnames(M_mat)[[1]][ind] <- list("Hoenig (1983) - Joint Equation")

        ind <- ind + 1
        ## Hoenig (1983) - Fish Equation
        M_mat[ind, 1]  <- round(exp(1.46 - 1.01 * log(tmax)), 3)
        dimnames(M_mat)[[1]][ind] <- list("Hoenig (1983) - Fish Equation")
    }
    if(any(method == "Lorenzen")){
        ind <- ind + 1
        ## Lorenzen (1996)
        M_mat[ind, 1]  <- round(3 * (Wwet^-0.288), 3)
        dimnames(M_mat)[[1]][ind] <- list("Lorenzen (1996)")
    }
    if(any(method == "Pauly_Linf")){
        ind <- ind + 1
        M_mat[ind, 1]  <- round(10^(-0.0066 - 0.279 * log10(Linf) + 0.6543 * log10(K_l) + 0.4634 * log10(temp)), 3)  ## exp( -0.0152 - 0.279 * log(Linf) + 0.6543 * log(K) + 0.463 * log(temp))
        dimnames(M_mat)[[1]][ind] <- list("Pauly (1980) - Length Equation")
        if(schooling == TRUE){
            M_mat[ind, 1] <- 0.8 * M_mat[ind, 1]
        }
    }
    if(any(method == "Pauly_Winf")){
        ind <- ind + 1
        M_mat[ind, 1]  <- round(10^(-0.2107 - 0.0824 * log10(Winf) + 0.6757 * log10(K_w) + 0.4627 * log10(temp)), 3)  ## exp( -0.2107 - 0.0824 * log(Winf) + 0.6757 * log(K) + 0.4627 * log(temp))
        dimnames(M_mat)[[1]][ind] <- list("Pauly (1980) - Weight Equation")
        if(schooling == TRUE){
            M_mat[ind, 1] <- 0.8 * M_mat[ind, 1]
        }
    }
    if(any(method == "PetersonWroblewski")){
        ind <- ind + 1
        ## Peterson and Wroblewski (1984)
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
        ## Roff (1984)
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

    print(M_mat)

    ret <- res
    par$M <- as.numeric(M_mat)
    ret$par <- par

    return(ret)
}

#' @title Catch curve
#'
#' @description  This function applies the (length-converted) linearised catch
#'    curve to age composition and length-frequency data,
#'    respectively. It allows to estimate the instantaneous total mortality rate (Z).
#'  p  Optionally, the gear selectivity can be estimated and the cumulative catch
#'    curve cna be applied.
#'
#' @param lfq a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} or \strong{age}: midpoints of the length classes (length-frequency
#'   data) or ages (age composition data),
#'   \item \strong{catch}: catches, vector or matrix with catches of subsequent years if
#'   the catch curve with constat time intervals should be applied;
#'  \item \strong{par}: a list with growth paramters:
#'  \itemize{
#'     \item \strong{Linf}: infinite length for investigated species in cm [cm],
#'     \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'     \item \strong{ta}: time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'     \item \strong{t0}: theoretical time zero, at which individuals of this species hatch,
#'     \item \strong{C}: amplitude of growth oscillation of soVBGF (range: 0 to 1, default: 0),
#'     \item \strong{ts}: summer point of soVBGF (ts = WP - 0.5) (range: 0 to 1, default: 0);
#'
#'  }
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
#' @param plot logical; should a plot be displayed? Default = TRUE
#'
#'
#' @keywords function mortality Z catchCurve
#'
#' @examples
#' \donttest{
#' #_______________________________________________
#' ## Variable paramter system (with catch vector)
#' ## based on length frequency data
#' data(goatfish)
#' output <- catchCurve(goatfish)
#' summary(output$linear_mod)
#'
#' ## based on age composition data
#' data(whiting)
#' catchCurve(whiting, catch_columns = 1)
#'
#' #_______________________________________________
#' ## Constant parameter system based on age composition data (with catch matrix)
#' catchCurve(whiting)
#'
#' #_______________________________________________
#' ## Cumulative Catch Curve
#' ## based on length frequency data
#' data(goatfish)
#' catchCurve(goatfish, cumulative = TRUE)
#'
#' ## based on age composition data
#' data(synCAA2)
#' catchCurve(synCAA2, cumulative = TRUE)
#'
#' #_______________________________________________
#' ## Catch Curve with estimation of selection ogive
#' data(synLFQ3)
#' output <- catchCurve(synLFQ3, calc_ogive = TRUE)
#' summary(output$linear_mod_sel)
#'  }
#'
#' ### the same with predefined selection for regression line:
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

catchCurve <- function(lfq,
  catch_columns = NA,
  cumulative = FALSE,
  calc_ogive = FALSE,
  reg_int = NULL,
  reg_num = 1,
  auto = FALSE,
  plot = TRUE
  ){

    res <- lfq
    if("par" %in% names(res)){
        par <- res$par
    }else par <- list()

    if("midLengths" %in% names(res)) classes <- as.character(res$midLengths)
    if("age" %in% names(res)) classes <- as.character(res$age)
    ## create column without plus group (sign) if present
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

    ## non cumulative catch curve
    if(cumulative == FALSE){
        if(is.na(catch_columns[1]) & constant_dt == FALSE) catch <- res$catch
        if(!is.na(catch_columns[1])){
            catchmat <- res$catch[,(catch_columns)]
            if(length(catch_columns) > 1){
                catch <- rowSums(catchmat, na.rm = TRUE)
            }else catch <- catchmat
        }
    }
    ##  cumulative catch curve
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

    ## Error message if catch and age do not have same length
    ##   Linearised catch curve with constant time intervals
    if(constant_dt){
        if("midLengths" %in% names(res) == TRUE) stop(noquote(
                                                     "The catch curve with constant time interval is not applicable to length-frequency data. Please provide a catch vector."))

        ## if(length(classes) != length(catch[,1])) stop(noquote(
        ##  "Age/length classes and catch matrix do not have the same length!"))

        if(length(classes) != length(diag(as.matrix(catch)))) writeLines("Age/length classes and the real cohort in the catch matrix \ndo not have the same length. The missing age/length \nclasses will be omitted.")

        ## Aged based Catch curve
        if("age" %in% names(res) == TRUE){
            ## find cohort to analyse
            real.cohort <- diag(as.matrix(catch))
            catch <- c(real.cohort, rep(NA,length(classes.num) - length(real.cohort)))
        }

    }else if(class(catch) == 'numeric'){
        if(length(classes) != length(catch)) stop(noquote(
                                                 "Age/length classes and catch vector do not have the same length!"))
    }

    ## Length converted catch curve
    if("midLengths" %in% names(res) == TRUE){

        if(!"par" %in% names(lfq)) stop(noquote("Please provide the required parameters in res$par!"))
        par <- res$par

        if("par" %in% names(res)){
            Linf <- res$par$Linf
            K <- res$par$K
            ta <- ifelse("ta" %in% names(res$par), res$par$ta, 0)
            t0 <- ifelse("t0" %in% names(res$par), res$par$t0, 0)
            C <- ifelse("C" %in% names(res$par), res$par$C, 0)
            ts <- ifelse("ts" %in% names(res$par), res$par$ts, 0)
        }else{
            Linf <- res$Linf
            K <- res$K
            ta <- ifelse("ta" %in% names(res), res$ta, 0)
            t0 <- ifelse("t0" %in% names(res), res$t0, 0)
            C <- ifelse("C" %in% names(res), res$C, 0)
            ts <- ifelse("ts" %in% names(res), res$ts, 0)
        }

        if((is.null(Linf) | is.null(K))) stop(noquote(
                                             "You need to assign values to Linf and K in lfq$par for the catch curve based on length-frequency data!"))

                                        #calculate size class interval
        midLengths <- classes.num
        interval <- midLengths[2] - midLengths[1]

        ## L and t of lower length classes
        lowerLengths <- midLengths - (interval / 2)
        if("C" %in% names(par) & "ts" %in% names(par)){
            t_L1 <- VBGF(pars = list(Linf = Linf, K = K, t0 = t0, C=C, ts=ts), L = lowerLengths)
        }else{
            t_L1 <- VBGF(pars = list(Linf = Linf, K = K, t0 = t0), L = lowerLengths)
        }
### t0 - (1/K) * log(1 - (lowerLengths / Linf))

        ## delta t
        dt <- rep(NA,length(midLengths))
        for(x1 in 1:(length(dt)-1)){
            dt[x1] <- t_L1[x1+1] - t_L1[x1]
        }

        ## x varaible
                                        #ln (Linf - L)
        ln_Linf_L <- log(Linf - lowerLengths)
### t of midlengths
        if("C" %in% names(par) & "ts" %in% names(par)){
            t_midL <- VBGF(pars = list(Linf = Linf, K = K, t0 = t0, C=C, ts=ts), L = midLengths)
        }else{
            t_midL <- VBGF(pars = list(Linf = Linf, K = K, t0 = t0), L = midLengths)
        }
### t0 - (1/K) * log(1 - (midLengths / Linf))

        ## y variable
                                        #ln C(L1,Linf)

        lnC <- log(catch)
        ## ln( Catch / delta t)
        lnC_dt <- log(catch / dt)
        lnC_dt[which(lnC_dt == -Inf)] <- NA   #### OR zero???

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

    ## Aged based Catch curve
    if("age" %in% names(res) == TRUE){
        ## delta t
        if(constant_dt == FALSE){
            dt <- rep(NA,length(classes.num))
            for(x1 in 1:(length(dt)-1)){
                dt[x1] <- classes.num[x1+1] - classes.num[x1]
            }
        }
        if(constant_dt) dt <- rep(classes.num[2] - classes.num[1], length(classes.num))

        ## x variable
        ## (t + dt) / 2   ==   x
        if(cumulative == FALSE) tplusdt_2 <- classes.num + (dt / 2)
        if(cumulative) tplusdt_2 <- classes.num

        ## y variable
                                        #ln C(L1,Linf)
        lnC <- log(catch)
        ## ln( Catch / delta t)     ==    y
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

### plot previous regression lines when using multiple regression lines
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

### save results to list
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

### define result lists
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

### calculations + model
        df.CC <- as.data.frame(cbind(xvar,yvar))
        df.CC.cut <- df.CC[cutter[1]:cutter[2],]
        lm1 <- lm(yvar ~ xvar, data = df.CC.cut)
        sum_lm1 <- summary(lm1)
        r_lm1 <- sum_lm1$r.squared
        intercept_lm1 <- sum_lm1$coefficients[1]
        slope_lm1 <- sum_lm1$coefficients[2]
        se_slope_lm1 <- sum_lm1$coefficients[4]

### fit of regression line
        lm1.fit <- sum_lm1$r.squared
        Z_lm1 <- abs(slope_lm1)
        SE_Z_lm1 <- abs(se_slope_lm1)
        confi <-  abs(se_slope_lm1) * qt(0.975,sum_lm1$df[2])
        conf_Z_lm1 <- Z_lm1 + c(-confi,confi)

### special case when cumulative and length-frequency data
        if(cumulative & "midLengths" %in% names(res) == TRUE){
            Z_lm1 <- Z_lm1 * K
            SE_Z_lm1 <- SE_Z_lm1 * K
        }

### save results to lists
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
                         Z_se = SE_Z_lm1List,
                         confidenceInt = conf_Z_lm1List))
        par$Z <- Z_lm1List
    }else{
        ret <- c(res,list(
                         xvar = xvar,
                         yvar = yvar,
                         reg_int = unlist(cutterList),
                         linear_mod = lm1List[[1]],
                         Z_se = unlist(SE_Z_lm1List),
                         confidenceInt = unlist(conf_Z_lm1List)))
        par$Z <- unlist(Z_lm1List)
    }

    if("M" %in% names(par) && length(par$M)==1){
        par$FM <- lapply(par$Z, function(x) x - par$M)
        par$E <- lapply(par$FM, function(x) x/par$Z)
            if(length(par$FM) == 1){
                par$FM <- unlist(par$FM)
                par$E <- unlist(par$E)
            }
    }
    names(ret)[names(ret) == "xvar"] <- xname
    names(ret)[names(ret) == "yvar"] <- yname
    class(ret) <- "catchCurve"

    ## Calculate selection ogive from catch curve and add to ret
    if(calc_ogive & cumulative) stop(noquote("It is not possible to estimate the selection ogive for the cumulative catch curve."))
    if(calc_ogive){

### Assumption that Z of smallest selected individuals is most appropriate
        mini <- min(unlist(cutterList))
        temp <- lapply(cutterList, function(x) grep(mini,x))
        ind <- sapply(temp, function(x) length(x) > 0)
        cutter <- unlist(cutterList[ind])


### only use part of catch and t which is not fully exploited by the gear
        t_ogive <- xvar[1:(cutter[1]-1)]
        dt_ogive <- dt[1:(cutter[1]-1)]
        if("age" %in% names(res) == TRUE &
           class(catch) == 'matrix' | class(catch) == 'data.frame'){
            catch_ogive <- catch[1:(cutter[1]-1)] ### catch.cohort
        }else catch_ogive <- catch[1:(cutter[1]-1)]


        ## calculate observed selection ogive
        Sobs <- catch_ogive/(dt_ogive * exp(unlist(intercept_lm1List[ind]) - unlist(Z_lm1List[ind]) * t_ogive))

        ## dependent vairable in following regression analysis
        ln_1_S_1 <- log((1/Sobs) - 1)

        ## get rid of Inf
        ln_1_S_1[which(ln_1_S_1 == Inf)] <- NA
        t_ogive[which(t_ogive == Inf)] <- NA

                                        #regression analysis to caluclate T1 and T2
        mod_ogive <- lm(ln_1_S_1 ~ t_ogive, na.action = na.omit)
        sum_lm_ogive <- summary(mod_ogive)
        T1 <- sum_lm_ogive$coefficients[1]
        T2 <- abs(sum_lm_ogive$coefficients[2])

        ## calculate estimated selection ogive
        Sest <- 1/(1+exp(T1 - T2*xvar))

        ## selection parameters
        t50 <- T1/T2
        t75 <- (T1 + log(3))/T2
        t95 <-  (T1 - log((1 / 0.95) - 1)) / T2
        if(!is.null(res$par$Linf) & !is.null(res$par$K)){
            if(is.null(res$par$t0)) t0 = 0
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
        ret2$par <- par

        class(ret2) <- "catchCurve"
        if(plot) plot(ret2, plot_selec=TRUE)
        return(ret2)
    }else {
        ret$par <- par
        if(plot) plot(ret)
        return(ret)
    }
}

#' @title Beverton & Holt's Z-Equations
                                        #
#' @description A method to estimate the instantaneous total mortality rate (Z) based
#'    on a method derived by Beverton and Holt (1956).
#'
#' @param lfq a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} or \strong{age}: midpoints of the length classes (length-frequency
#'   data) or ages (age composition data),
#'   \item \strong{catch}: catches, vector or matrix with catches of subsequent years if
#'   the catch curve with constat time intervals should be applied;
#'  \item \strong{par}: a list with growth paramters:
#'  \itemize{
#'     \item \strong{Linf}: infinite length for investigated species in cm [cm],
#'     \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'     \item \strong{ta}: time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'     \item \strong{t0}: theoretical time zero, at which individuals of this species hatch,
#'     \item \strong{C}: amplitude of growth oscillation of soVBGF (range: 0 to 1, default: 0),
#'     \item \strong{ts}: summer point of soVBGF (ts = WP - 0.5) (range: 0 to 1, default: 0);
#'
#'  }
#' }
#' @param catch_columns optional; in case catch is a matrix or data.frame, a number or vector
#'    indicating which column(s) of the matrix should be analysed (Default: \code{NA}).
#' @param Lprime_tprime length or age prime, above which all fish are under full exploitation as
#'    mid length or age class.
#'
#' @keywords function mortality Z
#'
#' @examples
#' ## based on length-frequency data
#' data(synLFQ2)
#' Z_BevertonHolt(synLFQ2, catch_columns = 2, Lprime_tprime = 47.5)
#'
#' ## based on age composition data
#' data(synCAA1)
#' Z_BevertonHolt(synCAA1, catch_columns = 3, Lprime_tprime = 2.5)
#'
#' @details  The first length group or age class within the list object \code{midLengths} or
#'    \code{age} will be used as the Lprim or tprime (length of recruitment to fishery).
#'
#' @return A list with the input parameters and following objects:
#' \itemize{
#'   \item \strong{tmean} or \strong{Lmean}: mean age or length of fish,
#'   \item \strong{tprime} or \strong{Lprime}: some age or length for which all fish of
#'      that length and longer are under full exploitation,
#'   \item \strong{Z}: total mortality.
#' }
#'
#' @references
#' Beverton R.J.H and S.J. Holt, 1956. A review of methods of estimating mortality rates
#' in exploited fish populations, with special reference to sources of bias in catch
#' sampling. \emph{Rapp.P.-v.Reun.CIEM}, 140:67-83
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

Z_BevertonHolt <- function(lfq,
                           catch_columns = NA,
                           Lprime_tprime){

    res <- lfq
    if("par" %in% names(res)){
        par <- res$par
    }else par <- list()
    catch <- res$catch

    if(class(catch) == "data.frame" | class(catch) == "matrix"){
        if(is.na(catch_columns[1])) stop("Please provide numbers indicating which column of the catch matrix should be analysed!")
        catchmat <- res$catch[,(catch_columns)]
        if(length(catch_columns) > 1){
            catch <- rowSums(catchmat, na.rm = TRUE)
        }else catch <- catchmat
    }

    ##   Length based equation
    if("midLengths" %in% names(res)){

        if(!"par" %in% names(res)) stop(noquote("Please provide the required parameters in res$par!"))

        classes <- as.character(res$midLengths)
        ## create column without plus group (sign) if present
        classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
        classes.num <- as.numeric(classes.num[,1])
        Lprime_tprime_ind <- which.min(abs(classes.num - Lprime_tprime))

        Linf <- res$par$Linf
        K <- res$par$K
        ta <- ifelse("ta" %in% names(res$par), res$par$ta, 0)
        t0 <- ifelse("t0" %in% names(res$par), res$par$t0, 0)
        C <- ifelse("C" %in% names(res$par), res$par$C, 0)
        ts <- ifelse("ts" %in% names(res$par), res$par$ts, 0)

        if((is.null(Linf) | is.null(K))) stop(noquote("You need to assign values to Linf and K to lfq$par for Z_BevertonHolt!"))


        ## Error message if catch and age do not have same length
        if(class(catch) == 'numeric'){
            if(length(classes) != length(catch)) stop("Ages and catch do not have the same length!")
        }

        ## calculate L prime
        ##Lprime <- classes.num[1] -
        ##  ((classes.num[2] - classes.num[1]) / 2)
        interval <- (classes.num[2] - classes.num[1]) / 2
        Lprime_tprime <- Lprime_tprime - interval


        ## calculate  C * (L1 + L2) / 2
        c_midlength <- catch * classes.num

        ## calculate L mean
        c_midlength_for_Lmean <- c_midlength[Lprime_tprime_ind:length(c_midlength)]
        catch_for_Lmean <- catch[Lprime_tprime_ind:length(catch)]
        Lmean <- sum(c_midlength_for_Lmean, na.rm = TRUE) / sum(catch_for_Lmean, na.rm = TRUE)

        Z <- K * (Linf - Lmean) / (Lmean - Lprime_tprime)

        ##save all in list
        ret <- c(res,list(Lmean = Lmean,Lprime = Lprime_tprime))
        par$Z <- Z
        if("M" %in% names(par) && length(par$M)==1){
            par$FM <- lapply(par$Z, function(x) x - par$M)
            par$E <- lapply(par$FM, function(x) x/par$Z)
            if(length(par$FM) == 1){
                par$FM <- unlist(par$FM)
                par$E <- unlist(par$E)
            }
        }
        ret$par <- par
        return(ret)
    }

    ##     Aged based equation
    if("midAge" %in% names(res) | "age" %in% names(res)){

        if("midAge" %in% names(res)) classes <- as.character(res$midAge)
        if("age" %in% names(res)) classes <- as.character(res$age)
        ## create column without plus group (sign) if present
        classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
        classes.num <- as.numeric(classes.num[,1])
        Lprime_tprime_ind <- which.min(abs(classes.num - Lprime_tprime))

        ## Error message if catch and age do not have same length
        if(class(catch) == 'numeric'){
            if(length(classes) != length(catch)) stop("Ages and catch do not have the same length!")
        }

        interval <- (classes.num[2] - classes.num[1]) / 2
        Lprime_tprime <- Lprime_tprime - interval
                                        #tprime <- classes.num[1] - interval
        catch_for_tprime <- catch[Lprime_tprime_ind:length(catch)]
        classes.num_for_tprime <- classes.num[Lprime_tprime_ind:length(classes.num)]
        sample.size <- sum(catch_for_tprime,na.rm=TRUE)
        sum.age.number <- sum((catch_for_tprime * classes.num_for_tprime), na.rm=TRUE)
        tmean <- sum.age.number/sample.size

        Z.BH <- 1 / (tmean - Lprime_tprime)

                                        #save all in list
        ret <- c(res,list(
                         tmean = tmean,
                         tprime = Lprime_tprime
                     ))
        par$Z <- Z.BH
        if("M" %in% names(par) && length(par$M)==1){
            par$FM <- lapply(par$Z, function(x) x - par$M)
            par$E <- lapply(par$FM, function(x) x/par$Z)
            if(length(par$FM) == 1){
                par$FM <- unlist(par$FM)
                par$E <- unlist(par$E)
            }
        }
        ret$par <- par
        return(ret)
    }

}

#' @title Estimate Z from CPUE data
#'
#' @description Method to estimate the instantaneous total mortality rate (Z) from
#'    catch per unit of effort (CPUE) data according to standard, Heincke's, or
#'    Robson & Chapman's method.
#'
#' @param cpue a list consisting of following parameters:
#' \itemize{
#'   \item \code{cohort}: a vector with with a cohort label,
#'   \item \code{age}: a vector with ages,
#'   \item \code{CPUE}: a vector with CPUE values;
#' }
#' @param method a character string indicating which assessment method should be used:
#'    \code{"standard"}, \code{"Heincke"}, or \code{"RobsonChapman"}.
#' @param omit_age1 logical; if \code{TRUE} the first age group is
#'    omitted (Default \code{FALSE}).
#'
#' @keywords function mortality Z CPUE
#'
#' @examples
#' ## load data
#' data(synCPUE)
#'
#' ## run model with standard method
#' Z_CPUE(synCPUE, method = "standard")
#'
#' ## run model with Heincke's method
#' Z_CPUE(synCPUE, method = "Heincke")
#'
#' ## run model with Robson and Chapman's method
#' Z_CPUE(synCPUE, method = "RobsonChapman", omit_age1 = TRUE)
#'
#' @details In Heincke's and RobsonChapman's method age groups older than 4 are lumped,
#'   because age groups older than 3 or 4 are said to be hard to seperate (Ricker, 1975).
#'   Sparre and Venema (1998) recommend to omit the first age group in case it is
#'   not fully exploited by the fishery.
#'
#' @return A list with input parameters and a Z value or matrix depending on the method.
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' Sparre, P., Venema, S.C., 1999. Introduction to tropical fish stock assessment.
#' Part 2. Excercises. \emph{FAO Fisheries Technical Paper}, (306.2, Rev. 2). 94 p.
#'
#' Ricker, W.E., 1975. Computation and interpretation of biological statistics of fish
#' populations. \emph{Bull.Fish.Res.Board Can.}, (191):382 p.
#'
#' @export

Z_CPUE <- function(cpue, method = "standard", omit_age1 = FALSE){

    res <- cpue
    if("par" %in% names(res)){
        par <- res$par
    }else par <- list()
    cohort <- res$cohort
    classes <- as.character(res$age)
    CPUE <- res$CPUE

    ## create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    switch(method,

           "standard" ={
               df.HZ <- data.frame(cohort = cohort[1:(length(cohort)-1)])
               result_Z <- list()
               for(i in 2:length(classes.num)){
                   Zi <- rep(NA,(length(classes.num)-1))
                   for(k in 1:(i-1)){
                       Zi[k] <- round((1 / (classes.num[i] - classes.num[k])) *
                                      (log(CPUE[k] /CPUE[i])),digits = 2)
                   }
                   result_Z[[i-1]] <- Zi
               }
               for(i in 1:length(result_Z)){
                   df.HZ[[paste(cohort[i+1])]] <- result_Z[[i]]
               }
               ret <- res
               par$Z <- df.HZ
               ret$par <- par
               return(ret)
           },

           "Heincke" ={
               CPUE[4] <- sum(CPUE[4:length(CPUE)])
               if(omit_age1) CPUE <- CPUE[-1]
               CPUE.H.n <- CPUE[2:length(CPUE)]
               CPUE.H.d <- CPUE[1:length(CPUE)]
               Z.H = - log( (sum(CPUE.H.n)) /
                            (sum(CPUE.H.d)))
               ret <- res
               par$Z <- Z.H
               if("M" %in% names(par) && length(par$M)==1){
                   par$FM <- lapply(par$Z, function(x) x - par$M)
                   par$E <- lapply(par$FM, function(x) x/par$Z)
                   if(length(par$FM) == 1){
                       par$FM <- unlist(par$FM)
                       par$E <- unlist(par$E)
                   }
               }
               ret$par <- par
               return(ret)
           },
           "RobsonChapman" ={
               CPUE[4] <- sum(CPUE[4:length(CPUE)])
               if(omit_age1) CPUE <- CPUE[-1]
               sum_CPUE.H.n <- sum(CPUE[2:length(CPUE)])
               sum_CPUE.H.d <- sum(CPUE[1:length(CPUE)])

               Z.H = - log((sum_CPUE.H.n) /
                           (sum_CPUE.H.d + sum_CPUE.H.n - 1))

               ret <- res
               par$Z <- Z.H
               if("M" %in% names(par) && length(par$M)==1){
                   par$FM <- lapply(par$Z, function(x) x - par$M)
                   par$E <- lapply(par$FM, function(x) x/par$Z)
                   if(length(par$FM) == 1){
                       par$FM <- unlist(par$FM)
                       par$E <- unlist(par$E)
                   }
               }
               ret$par <- par
               return(ret)
           })
}



#' Prepare length-frequency data for catch curve analysis
#'
#' @description The function can apply various methods of converting
#'   length-frequency data into aggregate numbers by relative age and
#'   (optionally) cohort. The method LCCC (length converted catch curve)
#'   does not account for seasonal oscillations in growth when determining age.
#'   Each length bin corresponds to given age according to the VBGF,
#'   and numbers are aggregated across length bins.
#'   With the GOTCHA method (Pauly, 1990), lfq data is sliced along growth
#'   trajectories to identify 'pseudocohorts', whose numbers are aggregated.
#'   As demonstrated by Pauly (1990), the method can be applied to cases where
#'   seasonally-oscillating growth creates inconsistencies of age at length
#'   depending on the time of year. The GOTCHA method aggregates by
#'   pseudocohort, as opposed to by length class in LCCC, and therefore does
#'   not require the correction of dividing numbers by time spent in a given
#'   length bin (n/dt).
#'   SLICC (slice-based catch curve) is an adaptation of GOTCHA whereby yearly
#'   cohorts are assigned but no aggregations are done. With this approach, one
#'   has the ability to regress numbers (i.e. individual bin counts) against
#'   both relative age (rel.age) and a factor relating to cohort. In principle,
#'   this should allow for the calculation of a single total mortality value
#'   (Z), while allowing for variable intercepts for each cohort, as might be
#'   assumed under conditions of variable recruitment.
#'
#' @param lfq an object of class lfq
#' @param method method of conversion (LCCC, GOTCHA, SLICC)
#' @param agemax integer. Maximum age plus group
#' @param n.cohort integer. Number of pseudocohort to derive from lfq.
#'   Only used GOTCHA method (Defaults to number of length bins for
#'   consistency with LCCC).
#' @param n.cohort.per.yr integer. Number of pseudocohorts per year.
#'   Overrides `n.cohort`. Can be applied to GOTCHA and SLICC. Default
#' @param use.ndt logical. Should numbers be divided by the time required to
#'   grow through bin. This argument was used for method development and is
#'   not expected to be defined by the user. For LCCC, this is set to TRUE,
#'   and set to FALSE for other methods.
#'
#' @return list.
#' @export
#'
#' @references
#' Pauly, D. (1990). Length-converted catch curves and the seasonal growth
#'   of fishes. Fishbyte, 8(3), 33-38.
#'
#' @examples
#'
#' data("synLFQ4")
#' lfq <- synLFQ4
#' lfq$par <- list(Linf = 80, K = 0.5, C = 0.75, ts = 0.5, ta = 0.25)
#' lfq <- lfqModify(lfq, year = 2006, bin_size = 2)
#' lfq <- lfqRestructure(lfq)
#' plot(lfq)
#'
#' # LCCC
#' res <- catchCurvePrep(lfq = lfq, method = "LCCC")
#' plot(log(n) ~ rel.age, res$tab)
#' incl <- which(res$tab$rel.age > 1)
#' points(log(n) ~ rel.age, res$tab[incl,], pch = 20 )
#' fit <- lm(log(n) ~ rel.age, res$tab[incl,])
#' abline(fit)
#' -coef(fit)[2] # true: Z = 1.0
#'
#' # GOTCHA (default settings)
#' res <- catchCurvePrep(lfq = lfq, method = "GOTCHA")
#' plot(log(n) ~ rel.age, res$tab)
#' incl <- which(res$tab$rel.age > 2 & res$tab$rel.age < 6)
#' points(log(n) ~ rel.age, res$tab[incl,], pch = 20 )
#' fit <- lm(log(n) ~ rel.age, res$tab[incl,])
#' abline(fit)
#' -coef(fit)[2] # true: Z = 1.0
#'
#' # GOTCHA (1 cohort per year)
#' res <- catchCurvePrep(lfq = lfq, method = "GOTCHA", n.cohort.per.yr = 1)
#' plot(log(n) ~ rel.age, res$tab)
#' incl <- which(res$tab$rel.age > 1.5 & res$tab$rel.age < 6)
#' points(log(n) ~ rel.age, res$tab[incl,], pch = 20 )
#' fit <- lm(log(n) ~ rel.age, res$tab[incl,])
#' abline(fit)
#' -coef(fit)[2] # true: Z = 1.0
#'
#' # SLICC
#' res <- catchCurvePrep(lfq = lfq, method = "SLICC")
#' plot(log(n) ~ rel.age, res$tab, col = res$tab$cohort)
#' incl <- which(res$tab$rel.age > 1 & res$tab$rel.age < 6)
#' points(log(n) ~ rel.age, res$tab[incl,], pch = 20,
#'   col = res$tab$cohort[incl])
#' fit0 <- lm(log(n) ~ rel.age + cohort, res$tab[incl,])
#' fit1 <- lm(log(n) ~ rel.age, res$tab[incl,])
#' fit <- get(c("fit0", "fit1")[which.min(AIC(fit0, fit1)$AIC)])
#' abline(fit)
#' -coef(fit)[2] # true: Z = 1.0
#'
#'
#'
#'
catchCurvePrep <- function(
  lfq,
  method = "LCCC",
  agemax = NULL,
  n.cohort = NULL, # only applies to GOTCHA1
  n.cohort.per.yr = NULL, # this should override n.cohort in GOTCHA
  use.ndt = NULL
  ){

  # replace zeros in catch mat
  lfq$catch[which(is.na(lfq$catch))] <- 0

  if(is.null(agemax)) agemax <- ceiling(VBGF(pars = lfq$par, L = lfq$par$Linf*0.95))

  if(is.null(use.ndt)){
    if(method == "LCCC"){
      use.ndt <- TRUE
    }else{
      use.ndt <- FALSE
    }
  }

  if(method == "LCCC"){
    lfqx <- lfq
    lfqx$par$C <- 0 # remove seasonality
    if(is.null(n.cohort.per.yr)){n.cohort.per.yr <- 36}
    lfqx <- lfqCohort(lfqx, n.cohort.per.yr = n.cohort.per.yr, calc_dt = TRUE, agemax = agemax)
    sumtab = data.frame(rel.age = apply(lfqx$rel.age, 1, mean, na.rm = TRUE))
    if(use.ndt){
      sumtab$n <- apply(lfqx$catch/lfqx$dt, 1, sum, na.rm = TRUE)
    }else{
      sumtab$n <- apply(lfqx$catch, 1, sum, na.rm = TRUE)
    }
    sumtab$length <- apply(array(lfqx$midLengths, dim = dim(lfqx$catch)), 1, mean, na.rm = TRUE)
    sumtab <- sumtab[which(sumtab$n > 0),] # remove rows where n == 0 | NA
  }

  if(method == "GOTCHA"){
    # if n.cohort.per.yr not NULL, then use this criteria for slicing
    if(!is.null(n.cohort.per.yr)) n.cohort <- NULL
    # if NULL, then use thin slices (n=36 per year) and aggregate later to n.cohort bins
    if(is.null(n.cohort.per.yr)){
      n.cohort.per.yr <- 36
      if(is.null(n.cohort)) n.cohort <- dim(lfq$catch)[1]
    }
    lfqx <- lfq
    lfqx <- lfqCohort(lfqx, calc_dt = TRUE, n.cohort.per.yr = n.cohort.per.yr, agemax = agemax)
    df <- data.frame(
      length = rep(lfqx$midLengths, times = length(lfqx$dates)),
      rel.age = c(lfqx$rel.age),
      bday = c(lfqx$bday),
      cohort = factor(c(lfqx$cohort)),
      dt = c(lfqx$dt))
    if(use.ndt){
      df$n <- c(lfqx$catch)/c(lfqx$dt)
    }else{
      df$n <- c(lfqx$catch)
    }
    df <- df[which(df$n > 0),] # remove rows where n == 0 | NA
    # create new slice aggregates if n.cohort is defined
    if(!is.null(n.cohort)){
      breaks <- seq(min(df$bday), max(df$bday), length.out = n.cohort)
      mids <- breaks[-1] - (diff(breaks)[1]/2)
      df$bday.cat <- cut(df$bday, breaks = breaks)
      df$bday2 <- mids[df$bday.cat]
      df$rel.age2 <- max(date2yeardec(lfqx$dates)) - df$bday2
      agg <- aggregate(n ~ rel.age2, df, sum, na.rm = TRUE)
      agg2 <- aggregate(length ~ rel.age2, df, max, na.rm = TRUE)
      sumtab <- data.frame(rel.age = agg$rel.age2, n = agg$n, length = agg2$length)
    }
    if(is.null(n.cohort)){
      agg <- aggregate(n ~ bday, df, sum, na.rm = TRUE)
      agg$rel.age <- max(date2yeardec(lfqx$dates)) - agg$bday
      agg2 <- aggregate(length ~ bday, df, max, na.rm = TRUE)
      sumtab <- data.frame(rel.age = agg$rel.age, n = agg$n, length = agg2$length)
    }
  }

  if(method == "SLICC"){
    lfqx <- lfq
    if(!is.null(n.cohort)) n.cohort <- NULL # always set n.cohort to NULL
    if(is.null(n.cohort.per.yr)) n.cohort.per.yr <- 1 # set to 1 cohort per year as default
    lfqx <- lfqCohort(lfqx, calc_dt = TRUE, n.cohort.per.yr = n.cohort.per.yr, agemax = agemax)
    df <- data.frame(
      length = rep(lfqx$midLengths, times = length(lfqx$dates)),
      rel.age = c(lfqx$rel.age),
      bday = c(lfqx$bday),
      cohort = factor(c(lfqx$cohort)),
      dt = c(lfqx$dt))
    if(use.ndt){
      df$n <- c(lfqx$catch)/c(lfqx$dt)
    }else{
      df$n <- c(lfqx$catch)
    }
    df <- df[which(df$n > 0),] # remove rows where n == 0 | NA
    agg <- aggregate(n ~ cohort + rel.age, df, sum, na.rm = TRUE)
    agg2 <- aggregate(length ~ cohort + rel.age, df, mean, na.rm = TRUE)
    sumtab <- data.frame(rel.age = agg$rel.age, cohort = agg$cohort, n = agg$n, length = agg2$length)
  }

  # order sumtab by rel.age
  sumtab <- sumtab[order(sumtab$rel.age),]

  res <- list(tab = sumtab, method = method, agemax = agemax,
    n.cohort = NULL, n.cohort.per.yr = n.cohort.per.yr,
    use.ndt = use.ndt)

  return(res)

}






#' Automated determination of values to include in catch curve
#'
#' @description Procedure described by Pauly (1990) to automate the selection
#'   of points used in a catch curve (\code{log(n) ~ rel.age}). The sequential
#'   data indices that result in the highest F-ststistic in a linear regression
#'   are returned. Results should still be visualized to prevent suboptimal
#'   results from 'pathological' datasets.
#'
#' @param fmla formula. Default is `formula(log(n) ~ rel.age)`. It is assumed
#'   that the main explanatory variable is named `rel.age`. Other formulas
#'   are allowed (e.g. `formula(log(n) ~ rel.age + cohort)`), but only the
#'   extent of rel.age is tested.
#' @param df data.frame holding the data.
#' @param minN numeric. Minimum number of values to include in the regression
#'   (Default: `minN = round(max(3,0.35*nrow(df)))`).
#' @param minP numeric. A value used as a trhreshhold for eliminating
#'   models with insignificant slopes (i.e. the statistical significance
#'   of the term `rel.age`)
#'
#' @return list containing the indices and F-statistic of the highest scoring
#' subset of (sequantial) data points.
#'
#' @export
#'
#' @references
#' Pauly, D. (1990). Length-converted catch curves and the seasonal growth
#'   of fishes. Fishbyte, 8(3), 33-38.
#'
#' @examples
#' # generate example data
#' set.seed(1)
#' n <- 10
#' rel.age <- seq(n)
#' logn <- 10 + -2*rel.age + rnorm(n, sd = 0.5)
#' logn[1:2] <- logn[1:2] * c(0.5, 0.75)
#' n <- exp(logn)
#' df <- data.frame(n = n, rel.age = rel.age)
#' plot(log(n) ~ rel.age, df)
#'
#' # When lm can reduce to 3 points, this is chosen as the 'best'
#' (tmp <- bestCC(df = df, minN = 3))
#' plot(f~n, data = tmp$bestbyn, t = "b", log="y")
#' plot(log(n)~rel.age, df)
#' points(log(n)~rel.age, df[tmp$best$samples,], pch = 16)
#'
#' # Any min number of points above this threshhold increases the best n a lot
#' (tmp <- bestCC(df = df))
#' plot(log(n)~rel.age, df)
#' points(log(n)~rel.age, df[tmp$best$samples,], pch = 16)
#' fit <- lm(tmp$fmla, data = df[tmp$best$samples,])
#' abline(fit)
#' -coef(fit)[2] # true = 2
#'
bestCC <- function(fmla = formula(log(n)~rel.age),
  df, minN = round(max(3,0.35*nrow(df))), minP = 0.01){

  if(minN < 3){stop("'minN' must be greater than or equal to 3")}

  # test if order of rel.age is sequential
  if( !all(diff(order(df$rel.age)) == 1) ){
    stop("df$rel.age values must be ordered sequentially")}

  # return unique combinations of sequential data points
  res <- vector("list", length(minN:nrow(df)))
  for(i in seq(minN:nrow(df))){
    len <- (minN:nrow(df))[i]
    res.i <- vector("list", length(1:(nrow(df)-len+1)))
    for(j in seq(1:(nrow(df)-len+1))){
      start <- (1:(nrow(df)-len+1))[j]
      res.i[[j]] <- seq(from = start, by = 1, length.out = len)
    }
    res[[i]] <- res.i
  }
  res <- do.call("c", res)
  res <- unique(res)

  # fit lm to each comb and return the indices and f-statistic
  res2 <- vector("list", length(res))
  pb <- txtProgressBar(min = 0, max = length(res), style = 3)
  for(i in seq(res2)){
    res2[[i]]$n <- length(res[[i]])
    res2[[i]]$samples <- res[[i]]
    df.i <- df[res[[i]],]
    fit <- lm(fmla, data = df.i)
    S <- summary(fit)
    res2[[i]]$f <- 0
    # if slope is negative and significant at minP, then record f-statistic
    if(coef(fit)[2] < 0 & S$coefficients[2,4] < minP){
      res2[[i]]$f <- S$fstatistic[1]
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)

  bestbyn <- do.call("rbind",
    lapply(res2, FUN = function(x){data.frame(n = x$n,f = x$f)}))
  bestbyn <- aggregate(f ~ n, data = bestbyn, FUN = max)
  # plot(f ~ n, bestbyn)

  # 'best' model is that with highest f-statistic
  best <- which.max(do.call("c", lapply(res2, FUN = function(x){x$f})))

  return(list(
    fmla = fmla,
    df = df,
    minN = minN,
    minP = minP,
    best = res2[[best]],
    bestbyn = bestbyn
  ))

}
