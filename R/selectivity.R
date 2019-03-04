#' @title Selectivity model
#
#' @description  This function estimates the selecitvity of gillnets and
#'    trawl nets from experimental catches.
#'
#' @param param a list with following parameters:
#' \itemize{
#'   \item \code{type} a string indicating which type of gear was used (options:
#'   "gillnet" or "trawl_net"),
#'   \item \code{midLengths} midlengths of size classes
#'   \item \code{meshSizes} a vector with mesh sizes in increasing order,
#'   \item \code{CatchPerNet_mat} a matrix with the catches per net in corresponding
#'   order of mesh sizes
#' }
#' @param plot logical; should the results be displayed graphically?
#'
#' @examples
#' # Gillnet selectivity
#' data(tilapia)
#' out <- select(param = tilapia)
#' plot(out)
#'
#' # Trawl selectivity
#' data(bream)
#' select(bream)
#'
#' @return A list with the input parameters and following list objects for type = "gillnet":
#' #' \itemize{
#'   \item \strong{classes.num}: numeric vector with length classes without a plus group,
#'   \item \strong{SNet1}: selection ogive net 1,
#'   \item \strong{SNet2}: selection ogive net 2,
#'   \item \strong{lnNet2_Net1}: logarithm ratio between nets,
#'   \item \strong{linear_mod}: linear model,
#'   \item \strong{LmNet1}: optimum length net 1,
#'   \item \strong{LmNet2}: optimum length net 2,
#'   \item \strong{SF}: selection factor,
#'   \item \strong{stand.dev}: standard deviation of selection factor.
#' }
#' and following objects for type = "trawl_net":
#' \itemize{
#'   \item \strong{classes.num}: numeric vector with length classes,
#'   \item \strong{SLobs}: observed selection ogive,
#'   \item \strong{SLest}: estimated selection ogive,
#'   \item \strong{lnSL}: logarithm of observed selection,
#'   \item \strong{linear_mod}: linear model,
#'   \item \strong{S1}: constant of selection curve,
#'   \item \strong{S2}: another constant of selection curve,
#'   \item \strong{L25}: length at which 25% of the fish are retained in the codend,
#'   \item \strong{L50}: length at which 50% of the fish are retained in the codend,
#'   \item \strong{L75}: length at which 75% of the fish are retained in the codend,
#'   \item \strong{SF}: selection factor.
#' }
#'
#' @details This function estimates the fractions retained by each net, the
#'   optimum lengths for each net, the selection factor (SF), and the standard deviation
#'   of the factor (stand.dev). Calculations are based on a normal distribution with common spread.
#'   Assumptions of this method are, that (i) the optimum length Lm is proportional to the mesh
#'   size (Lm = SF * m), (ii) the selection curves are normally distributed with a common
#'   standard deviation, (iii) the nets have the same fishing power (same dimensions and material).
#'   Requirements for the experimental set-up are: selection curves corresponding to the two
#'   mesh sizes have to overlap, and the nets have to be set in the same area, during the
#'   same time.
#'   To calculate selection factor (SF), L25, L50 and L75 for trawl nets /fisheries.
#'
#' @importFrom graphics plot
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

select <- function(param, plot = FALSE){
  res <- param
  if(is.null(res$midLengths)) stop("The method requires midpoints of length classes (midLengths).")
  classes <- as.character(res$midLengths)

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
  classes.num <- as.numeric(classes.num[,1])

  #  Gillnet
  if(is.null(res$type)) stop(noquote("Please indicate if a gillnet or trawl selectivity experiment was peroformed (see 'param$type')."))
  if(res$type == "gillnet"){
    numNet1 <- res$CatchPerNet_mat[,1]
    numNet2 <- res$CatchPerNet_mat[,2]
    msNet1 <- res$meshSizes[1]
    msNet2 <- res$meshSizes[2]

    # log ratios y = ln(numNet2 / numNet1)
    lnNet2_Net1 <- log(numNet2 / numNet1)
    lnNet2_Net1[which(numNet1 == 0 | numNet2 == 0 |
                        lnNet2_Net1 == Inf | lnNet2_Net1 == -Inf)] <- NA

    #regression analysis
    mod <- lm(lnNet2_Net1 ~ classes.num)
    sum.mod <- summary(mod)
    a <- sum.mod$coefficients[1]
    b <- sum.mod$coefficients[2]

    # SF
    SF <- (- 2 * a) / (b * (msNet1 + msNet2))
    LmNet1 <- SF * msNet1
    LmNet2 <- SF * msNet2

    # standard deviation
    s2 <- SF * ((msNet2 - msNet1) / b)
    s <- sqrt(s2)

    # points on selection curves
    SNet1 <- exp(-((classes.num - LmNet1)^2 / (2 * s2)))
    SNet2 <- exp(-((classes.num - LmNet2)^2 / (2 * s2)))

    # numbers in population
    NNet1 <- numNet1 / SNet1
    NNet2 <- numNet2 / SNet2

    res2 <- list(classes.num=classes.num,
                 SNet1=SNet1,
                 SNet2=SNet2,
                 lnNet2_Net1 = lnNet2_Net1,
                 linear_mod = mod,
                 LmNet1 = LmNet1,  #Lm optimum length
                 LmNet2 = LmNet2,
                 SF = SF,
                 stand.dev = s)
  }

  # Trawl
  if(res$type == "trawl_net"){
    numCodend <- res$CatchPerNet_mat[,2]
    numCover <- res$CatchPerNet_mat[,1]
    meshsizeCodend <- res$meshSizes[2]

    # calculate fraction retained (SL obs)
    SLobs <- numCodend/(numCodend + numCover)

    # ln( 1 / SL - 1)
    lnSL <- log(1/SLobs - 1)

    #excluding point where no (0) or full (1) retention was obtained and all values beyond those points, even if they are between 0 and 1
    lnSL[which(classes.num >= classes.num[which(SLobs == 1)[1]])] <- NA   ##many ones in data, use first one?
    lnSL[which(lnSL == Inf | lnSL == -Inf)] <- NA

    #model
    mod <- lm(lnSL ~ classes.num)
    sum.mod <- summary(mod)
    S1 <- sum.mod$coefficients[1]
    S2 <- abs(sum.mod$coefficients[2])

    #L25,L50,L75
    L25 <- (S1 - log(3)) / S2
    L50 <- S1 / S2
    L75 <- (S1 + log(3)) / S2

    #estimated SL (SL est)
    SLest <- 1 / (1 + exp(S1 - S2 * classes.num))

    #Selection factor
    SF <- L50/meshsizeCodend

    res2 <- list(classes.num = classes.num,
                 SLobs = SLobs,
                 SLest = SLest,
                 lnSL = lnSL,
                 linear_mod = mod,
                 S1 = S1,
                 S2 = S2,
                 L25 = L25,
                 L50 = L50,
                 L75 = L75,
                 SF = SF)
  }
  ret <- c(res,res2)
  class(ret) = 'select'
  if(plot) plot(ret)
  return(ret)
}

#' @title Selectivity patterns
#'
#' @description  Based on a few parameters, this function estimates the fraction per
#'    length group retained in the net. Different selection curves can be used for the
#'    estimation.
#'
#' @param s_list a list with selectivity parameters dependent on the type of
#'    selection curve:
#' \itemize{
#'   \item \code{selecType} type of
#'    selection curve used for estimation (options:
#'    "knife_edge",
#'    "trawl_ogive",
#'    "lognormal",
#'    "normal_fixed"),
#'   \item \code{Lc} length-at-first-capture (also called L50),
#'   \item \code{meshSizes} a vector with mesh sizes in increasing order,
#'   \item \code{select_p1} selectivity parameter 1 (see Millar and Holst (1997)),
#'   \item \code{select_p2} selectivity parameter 2 (see Millar and Holst (1997)),
#'   \item \code{L75} length at which individuals are caught with a
#'    probability of 75%.
#' }
#' @param Lt a vector with lengths corresponding to age classes
#' @param Lc length-at-first-capture (Default: NA)
#'
#' @examples
#' # create list with selectivity information
#' select.list <- list(selecType = 'knife_edge',
#'    Lc = 34, L75 = 37, tc = 5, meshSizes = c(60,80),
#'    select_p1 = 2.7977, select_p2 = 0.1175)
#'
#'
#' # create vector with mid lengths
#' Lt <- seq(5, 50, 0.01)
#'
#' # knife edge selectivity
#' sel_ke <- select_ogive(select.list, Lt)
#'
#' # trawl ogive selectivity
#' select.list$selecType = "trawl_ogive"
#' sel_to <- select_ogive(select.list, Lt)
#'
#' plot(Lt, sel_ke, type = 'l')
#' lines(Lt, sel_to, col = 'blue')
#'
#'
#' # Gillnet selectivity ("lognormal" and "normal_fixed")
#' select.list$selecType <- "lognormal"
#' sel_log <- select_ogive(select.list, Lt)
#'
#' select.list$selecType <- "normal_fixed"
#' select.list$select_p1 <- 0.2
#' select.list$select_p2 <- 1.5
#' sel_nf <- select_ogive(select.list, Lt)
#'
#' plot(Lt, sel_log, type = 'l')
#' lines(Lt, sel_nf, col = 'blue')
#'
#' @details This function is embedded within \code{\link{predict_mod}}. \code{selecType}
#'    "knife_edge" only requires a Lc value. "trawl_ogive" requires a Lc (L50) and
#'    a L75 value. "lognormal" requires two mesh sizes, an estimate of mu and of sigma.
#'    "normal_fixed" requires two mesh sizes with an estimate of the selection factor (SF) and an
#'    estimate of sigma.
#'
#' @references
#' Millar, R. B., Holst, R. (1997). Estimation of gillnet and hook selectivity using log-linear
#' models. ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

select_ogive <- function(s_list, Lt, Lc = NA){
  type <- s_list$selecType

  switch(type,

         'knife_edge' = {
           if(is.na(Lc)){
             if("L50" %in% names(s_list)){L50 <- s_list$L50
             }else L50 <- s_list$Lc
           }else L50 <- Lc
           sel <- rep(0, length(Lt))
           sel[as.numeric(as.character(Lt)) >= L50] <- 1
         },

         'trawl_ogive' = {
           # if no new Lc (L50) value is provided:
           if(is.na(Lc)){
             if("L50" %in% names(s_list)){L50 <- s_list$L50
             }else L50 <- s_list$Lc

             L75 <- s_list$L75
           }
           # if new Lc (L50) value is provided:
           if(!is.na(Lc)){
             if("L50" %in% names(s_list)){L50old <- s_list$L50
             }else L50old <- s_list$Lc
             L75old <- s_list$L75
             L50 <- Lc
             # new L75 based on relation of L50 values
             L75 <- L75old * (L50/L50old)
           }

           sel <- 1 / (1 + exp(- (Lt - L50)/
                                 ((2*(L75 - L50))/(log(0.75/(1-0.75))-
                                                     log(0.25/(1-0.25)))))) #or: sel <- 1 / (1 + exp(S1.TS - S2.TS * Lt))

         },

         "lognormal" = {
           mu <- s_list$select_p1
           sigma <- s_list$select_p2
           meshSizes <- s_list$meshSizes
           sel <- (1/Lt) * exp(mu + log(meshSizes[2]/meshSizes[1]) -
                                 (sigma^2/2) - (((log(Lt) - mu -
                                      log(meshSizes[2]/meshSizes[1]))^2)/(2 *sigma^2)))
         },

         "normal_fixed" = {
           SF <- s_list$select_p1
           sigma <- s_list$select_p2
           meshSizes <- s_list$meshSizes
           sel <- exp(-((Lt - meshSizes[1] * SF)^2/(2 * sigma^2)))
         },

         stop(paste("\n",type, "not recognised, possible curve types are \n",
                    "\"knife_edge\", \"trawl_ogive\", \"lognormal\" \n",
                    " and \"normal_fixed\""))
  )
  return(sel)
}


#' @title Millar's selectivity model
#
#' @description  This model estimates the selecitvity of different gears from experimental
#'    catches.
#'
#' @param param A list with following parameters: vector with midlengths of size classes
#'    (\code{midLengths}), vector with meshSizes in increasing order (\code{meshSizes}),
#'    and a matrix with the number of individuals caught with each sized mesh
#'    (\code{CatchPerNet_mat}).
#' @param x0 A string of initial values for the parameters to be optimized over when applying the
#'    function \code{\link[stats]{optim}}. When rtype = "norm.loc", "norm.sca", or "lognorm", initial values
#'    may be ommitted (\code{x0 = NULL}), and starting vallues will be estimated using the earlier approach of
#'    Millar and Holst (1997).
#' @param rtype A character string indicating which method for estimating selection curves
#'    should be used:
#'    \code{"norm.loc"} for a normal curve with common spread,
#'    \code{"norm.sca"} for a normal curve with variable spread,
#'    \code{"lognorm"} for a lognormal curve,
#'    \code{"binorm.sca"} for a bi-normal curve,
#'    \code{"bilognorm"} for a bi-lognormal curve,
#'    \code{"tt.logistic"} for a control and logistic curve
#' @param rel.power A string indicating the relative power of different meshSizes,
#'    must have same length as \code{meshSizes} (Default: \code{rel.power = NULL}).
#' @param plot logical; should a plot be printed?
#'
#'
#' @examples
#' data(haddock)
#'
#' output <- select_Millar(haddock, x0 = c(-10,0.3,0),
#'    rtype = "tt.logistic")
#'
#' plot(output, plotlens=seq(25,35,0.1), deviance_plot = FALSE)
#' legend("topleft",c("Control","Experimental"), lty=1:2, col=1:2)
#'
#'
#' # Gillnet
#' data(gillnet)
#'
#' # Using inital estimates from old method
#' select_Millar(gillnet, x0 = NULL, rtype = "norm.loc")$value
#' select_Millar(gillnet, x0 = NULL, rtype = "norm.sca")$value
#' select_Millar(gillnet, x0 = NULL, rtype = "lognorm")$value
#'
#' # Two rtypes which require starting values
#' select_Millar(gillnet, x0 = c(55,4,65,4,3), rtype="binorm.sca")
#' select_Millar(gillnet, x0 = c(4,0.2,4.2,0.1,2), rtype="bilognorm")
#'
#' # Calculation with finer length resolution
#' output <- select_Millar(gillnet, x0 = NULL, rtype = "lognorm")
#' plot(output, plotlens=seq(40,90,0.1))
#'
#' # Use alternate plot settings
#' output <- select_Millar(gillnet, x0 = NULL, rtype = "lognorm")
#' ncolor <- length(output$meshSizes)
#' plot(output, plotlens=seq(40,90,0.1), deviance_plot = FALSE,
#'  lty=1, col=rainbow(ncolor))
#' legend("topleft", col=rainbow(ncolor), legend=output$meshSizes,
#'  lty=1, title="Mesh size [cm]")
#'
#'  # deviance plot only
#'  plot(output, plotlens=seq(40,90,0.1), selectivity_plot = FALSE)
#'
#'
#' # Stacked trammel net
#' # The data come from two experiments using different mesh sizes
#' # This analysis assumes common retention curve in both experiments.
#' # Note that summary function does not produce residual plot
#' # since lengths are not unique
#' data(trammelnet)
#'
#' output <- select_Millar(trammelnet, x0 = c(25,4),
#'    rtype="norm.loc", rel.power = rep(1,6))
#'
#' plot(output,plotlens=seq(10,40,0.1))
#'
#' @source https://www.stat.auckland.ac.nz/~millar/selectware/
#'
#' @details Model adapted from the selectivity functions provided by Prof. Dr. Russell Millar
#'   (https://www.stat.auckland.ac.nz/~millar/). In the deviance plot open circles correspond to negative,
#'   closed to positive residuals. The size of the circles is proportional to the square of the residuals.
#'   To assess the model fit by the deviance plot it requires some experience, in general the pattern should
#'   be random and the sizes not too big. Please refer to Millar's publications and other publications for
#'   comparison. The model can produce errors if the starting values (\code{x0}) for the \link{optim}
#'   function are not realistic. Please be aware that if the method is changed the outcoming parameters
#'   can greatly vary. Simliarly the starting values have to be adapted when changing the method (\code{rtype}).
#'
#' @importFrom graphics plot
#' @importFrom stats deviance optim
#' @import msm
#'
#' @references
#'  Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#'  54(3):471-477
#'
#'  Holt, S. J. 1963. A method for determining gear selectivity and its application.
#'  \emph{ICNAF Special Publication}, 5: 106-115.
#'
#'
#' @export

select_Millar <- function(param,
                          x0 = NULL,
                          rtype = "norm.loc",
                          rel.power = NULL, plot = TRUE
                          ){
# Adapted R code from Russell Millar (https://www.stat.auckland.ac.nz/~millar/selectware/)

  res <- param
  meshSizes <- res$meshSizes
  mSize1 <- meshSizes[1]
  classes <- res$midLengths
  CatchPerNet_mat <- res$CatchPerNet_mat

  if(is.null(CatchPerNet_mat)) CatchPerNet_mat <- cbind(res$numCover,res$numCodend)
  res$CatchPerNet_mat <- CatchPerNet_mat

  # NetFit
  if(sum(sort(meshSizes)==meshSizes) != length(meshSizes))
    stop("Mesh size must be in ascending order!")

  if(is.null(rel.power)) rel.power <- rep(1,length(meshSizes))

  if(!is.null(meshSizes) & ncol(CatchPerNet_mat) != length(meshSizes))
    stop("Number of mesh sizes should be ", ncol(CatchPerNet_mat))

  CountPropns <- CatchPerNet_mat/apply(CatchPerNet_mat,1,sum,na.rm=TRUE)

  fullfit.l <- sum(CatchPerNet_mat*log(CountPropns),na.rm=TRUE)

  r <- rtypes_Millar(rtype) #Get selection curve function


  # nllhood
  nllhood <- function(theta, classes, CatchPerNet_mat, meshSizes, r, rel.power){
    rmatrix <- outer(classes, meshSizes, r, theta)
    rmatrix[is.na(CatchPerNet_mat)] <- NA #No fitted retention for missing meshSizes
    rmatrix <- t(t(rmatrix) * rel.power)
    phi <- rmatrix/apply(rmatrix,1,sum,na.rm=TRUE)
    nll <- - sum(CatchPerNet_mat * log(phi),na.rm=TRUE)
    return(nll)
  }

  # fit old gillnet selectivity approach of Millar for initial values and pass to x0 if needed
  if(is.null(x0) & rtype %in% c("norm.loc", "norm.sca", "lognorm")){
    dat.tmp <- cbind(classes, CatchPerNet_mat)
    gillnetfit.res <- gillnetfit(data=dat.tmp, meshsizes=meshSizes, rtype=rtype, rel.power=rel.power,
      plotlens=NULL, details=TRUE)
    if(rtype %in% c("norm.loc", "norm.sca")){
      x0 <- unname(gillnetfit.res$gear.pars[3:4,1])
    }
    if(rtype %in% c("lognorm")){
      x0 <- unname(gillnetfit.res$gear.pars[1:2,1])
    }
  }else if(is.null(x0) & !(rtype %in% c("norm.loc", "norm.sca", "lognorm"))){
    stop(paste0("You have to provide starting values for ",rtype))
  }

  res2 <- optim(x0, nllhood, classes = classes, CatchPerNet_mat = CatchPerNet_mat,
               meshSizes = meshSizes, r = r, rel.power = rel.power,
               hessian = TRUE, control = list(trace = FALSE))

    #res2

    #cat("Parameters =", res2$par, ",  Deviance =", 2*(fullfit.l + res2$value), "\n")

    res3 <- c(res, res2, deviance=deviance, rtype=rtype, rel.power=list(rel.power))  # invisible


    # Estimates
    x <- res3$par
    if(det(res3$hess) == 0) stop("It is impossible to invert the matrix. Please change starting parameters (x0) or the method (rtype) and try again.")
    varx <- solve(res3$hess)
    names <- c("Mode(mesh1)","Std dev.(mesh1)")
    switch(res3$rtype,
           "norm.loc"={
             pars <- x
             varpars <- varx
             },
           "norm.sca"={
             pars <- x
             varpars <- varx
             },
           "lognorm"={
             pars <- c(exp(x[1]-x[2]^2),
                       sqrt(exp(2*x[1]+x[2]^2)*(exp(x[2]^2)-1)))
             varpars <- msm::deltamethod(list(~exp(x1-x2^2),
                                      ~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1))), x, varx, ses=FALSE)
             },
           "binorm.sca"={
             pars <- c(x[1:4],exp(x[5])/(1+exp(x[5])))
             names <- c("Mode1(mesh1)","Std dev.1(mesh1)",
                     "Mode2(mesh1)","Std dev.2(mesh1)","P(mode1)")
             varpars <- msm::deltamethod(list(~x1,~x2,~x3,~x4,~exp(x5)/(1+exp(x5))),
                                 x,varx,ses=F)
             },
           "bilognorm"={
             pars <- c(exp(x[1]-x[2]^2),sqrt(exp(2*x[1]+x[2]^2)*(exp(x[2]^2)-1)),
                    exp(x[3]-x[4]^2),sqrt(exp(2*x[3]+x[4]^2)*(exp(x[4]^2)-1)),
                    exp(x[5])/(1+exp(x[5])))
             names <- c("Mode1(mesh1)","Std dev.1(mesh1)",
                     "Mode2(mesh1)","Std dev.2(mesh1)","P(mode1)")
             varpars <- msm::deltamethod(
               list(~exp(x1-x2^2),~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1)),
                    ~exp(x3-x4^2),~sqrt(exp(2*x3+x4^2)*(exp(x4^2)-1)),
                    ~exp(x5)/(1+exp(x5))),x,varx,ses=FALSE)
             },
           "tt.logistic"={
             pars <- c(-x[1]/x[2],2*(log(3))/x[2],exp(x[3])/(1+exp(x[3])))
             names <- c("L50","SR","p")
             varpars <- msm::deltamethod(list(~-x1/x2,~2*log(3)/x2,~exp(x3)/(1+exp(x3))),
                                 x,varx,ses=FALSE)
             },
           stop(paste("\n",res3$rtype, "not recognised, possible curve types are \n",
                      "\"norm.loc\", \"norm.sca\", \"lognorm\" \n",
                      "\"binorm.sca\", \"bilognorm\", and \"tt.logistic\""))
    )

    estimates <- cbind(pars,sqrt(diag(varpars)))
    colnames(estimates) <- c("par", "s.e.")
    rownames(estimates) <- names

    # Summary
    nlens <- length(classes)
    meshSizes <- res3$meshSizes
    nmeshes <- length(meshSizes)

    O <- res3$CatchPerNet_mat #Matrix of observed classes

    rmatrix <- outer(classes,meshSizes,r,res3$par)
    rmatrix[is.na(O)] <- NA #No fitted retention for missing meshSizes
    rmatrix <- t(t(rmatrix)*res3$rel.power)
    phi <- rmatrix/apply(rmatrix,1,sum,na.rm=TRUE)
    E <- apply(O,1,sum,na.rm=TRUE)*phi #Matrix of expected classes
    Pearson.resids <- (O-E)/sqrt(E)
    Pearson.chisq <- sum(Pearson.resids^2,na.rm=TRUE)
    wk <- O*log(O/E); wk[is.na(wk)]=0
    Dev.resids <- sign(O-E)*sqrt(2*(E-O+wk))
    Deviance <- sum(Dev.resids^2,na.rm=TRUE)
    full.l <- sum(-O+O*log(O),na.rm=TRUE)
    null.E <- matrix(apply(O,1,mean,na.rm=TRUE),nrow=nlens,ncol=nmeshes)
    null.l <- sum(-null.E+O*log(null.E),na.rm=TRUE)
    model.l <- sum(-E+O*log(E),na.rm=TRUE)
    NonZeroDat <- O[apply(O,1,sum,na.rm=TRUE)>0,]
    d.o.f. <- nrow(NonZeroDat)*(nmeshes-1)-length(res3$par)-sum(is.na(NonZeroDat))
    out <- rbind(null.l, model.l, full.l, Deviance, Pearson.chisq, d.o.f.)
    AreLensUnique <- (length(classes)==length(unique(classes)))

    lenrmatrix <- cbind(classes,rmatrix) ##cbind(plotlens,rmatrix)
    colnames(lenrmatrix) <- c("Length",meshSizes)  #invisible(lenrmatrix)

    # return results
    ret <- c(res3, list(Dev.resids = Dev.resids,
      selection_ogive_mat = lenrmatrix,
      out = out,
      estimates = estimates)
    )


  class(ret) <- "select_Millar"
  if(plot) plot(ret)
  return(ret)
}



#' @title Millar's selectivity types
#
#' @description  This function returns a function corresponding to the type of curve
#'    which was selected to represent the selectivity of nets or hooks.
#'
#' @param rtype a character string indicating which method for the estimation of selection
#'    curves should be used:
#'    \code{"norm.loc"} for normal with common spread method,
#'    \code{"norm.sca"} for normal with variable spread method,
#'    \code{"lognorm"} for lognormal method,
#'    \code{"binorm.sca"} for bi-normal method,
#'    \code{"bilognorm"} for bi-lognormal method,
#'    \code{"tt.logistic"} for control and logistic method,
#'    \code{"gamma"} for gamma method.
#'
#' @source https://www.stat.auckland.ac.nz/~millar/selectware/
#'
#' @details Function adapted from the selectivity functions provided by
#'   Prof. Dr. Russell Millar (https://www.stat.auckland.ac.nz/~millar/).
#'   Until now following curves are incorporated:
#'   \code{"norm.loc"} for a normal curve with common spread,
#'    \code{"norm.sca"} for a normal curve with variable spread,
#'    \code{"lognorm"} for a lognormal curve,
#'    \code{"binorm.sca"} for a bi-normal curve,
#'    \code{"bilognorm"} for a bi-lognormal curve,
#'    \code{"tt.logistic"} for a control and logistic curve.
#'
#' @references
#'  Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#'  54(3), 471-477.

rtypes_Millar <- function(rtype) {
# Adapted R code from Russell Millar (https://www.stat.auckland.ac.nz/~millar/selectware/)
  switch(rtype,

         "norm.loc" = {
           r <- function(classes, meshSizes, th) {
             relsize <- meshSizes / meshSizes[1]
             seln <- exp(-(classes - th[1] * relsize)^2 / (2 * th[2]^2))
             return(seln)
           }
         },

         "norm.sca" = {
           r <- function(classes, meshSizes, th) {
             relsize <- meshSizes / meshSizes[1]
             seln <- exp(-(classes - th[1]*relsize)^2 / (2 * th[2]^2 * relsize^2))
             return(seln)
             }
           },

         "lognorm" = {
           r <- function(classes, meshSizes, th) {
             relsize <- meshSizes / meshSizes[1]
             seln <- (relsize / classes) * exp(th[1] - th[2]^2 / 2)
             seln <- seln * exp( -(log(classes) - th[1] - log(relsize))^2 / (2*th[2]^2))
             return(seln)
             }
           },

         "binorm.sca" = {
           r <- function(classes, meshSizes, th) {
             relsize <- meshSizes / meshSizes[1]
             seln1 <- exp(-(classes - th[1] * relsize)^2 / (2 * th[2]^2 * relsize^2))
             seln2 <- exp(-(classes - th[3] * relsize)^2 / (2 * th[4]^2 * relsize^2))
             p <- exp(th[5]) / (1 + exp(th[5])) #i.e., th[5]=logit(p)
             seln <- p * seln1 + (1 - p) * seln2
             return(seln)
             }
           },

         "bilognorm" = {
           r <- function(classes, meshSizes, th) {
             relsize <- meshSizes / meshSizes[1]
             seln1 <- (relsize / classes) * exp(th[1] - th[2]^2 / 2)
             seln1 <- seln1*exp( -(log(classes) - th[1] - log(relsize))^2 / (2 * th[2]^2))
             seln2 <- (relsize / classes) * exp(th[3] - th[4]^2 / 2)
             seln2 <- seln2 * exp( -(log(classes) - th[3] - log(relsize))^2 / (2 * th[4]^2))
             p <- exp(th[5]) / (1 + exp(th[5])) #i.e., th[5]=logit(p)
             seln <- p * seln1 + (1 - p) * seln2
             return(seln)
             }
           },

         "tt.logistic" = {
           r <- function(classes, meshSizes, th) {
             control <- (meshSizes == meshSizes[1])
             p <- exp(th[3]) / (1 + exp(th[3])) #i.e., th[3]=logit(p)
             wk <- exp(th[1] + th[2] * classes)
             lselect <- wk / (1 + wk)
             seln <- (1 - p) * control + p * lselect * (1 - control)
             return(seln)
             }
           },
         "gamma" = {
           r <- function(classes, meshSizes, th){
             seln <- (classes / ((th[1] - 1) * th[2] * meshSizes)) ^ (th[1] - 1) *
               exp(th[1] - 1 - classes / (th[2] * meshSizes))
             return(seln)
           }
         },

         stop(paste("\n",rtype, "not recognised, possible curve types are \n",
                    "\"norm.loc\", \"norm.sca\", \"lognorm\" \n",
                    "\"binorm.sca\", \"bilognorm\", \"tt.logistic\" and \"gamma\""))
  )
  return(r)
}



#' Predict gillnet selectivity (old Millar method)
#'
#' @param type ("norm.loc", "norm.sca", "gamma", "lognorm")
#' @param meshsizes mesh sizes
#' @param rel relative powers (one for each mesh size)
#' @param pars selection curve parameters
#' @param plotlens vector of new mesh sizes for selectivity prediction
#'
#' @source https://www.stat.auckland.ac.nz/~millar/selectware/
#'
#' @references
#'  Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#'  54(3):471-477
#'
#' @return selectivities

rcurves_Millar <- function(type, meshsizes, rel, pars, plotlens){
# Adapted R code from Russell Millar (https://www.stat.auckland.ac.nz/~millar/selectware/)

  lens <- rep(plotlens,length(meshsizes))
  relsizes <- meshsizes/meshsizes[1]

  if(is.null(rel)){
    releff <- 1
  }else releff <- rep(rel,rep(length(plotlens),length(meshsizes)))

  msizes <- rep(meshsizes,rep(length(plotlens),length(meshsizes)))
  relsizes <- rep(relsizes,rep(length(plotlens),length(relsizes)))

  switch(type,
   "norm.loc"={ k=pars[1]; mean1=pars[3]; sigma=pars[2]
     rselect <- exp(-((lens-k*msizes)^2)/(2*sigma^2)) },
   "norm.sca"={ k1=pars[1]; k2=pars[2]
     rselect <- exp(-((lens-k1*msizes)^2)/(2*k2*msizes^2)) },
   "gamma"={ alpha=pars[1]; k=pars[2]
     rselect <- (lens/((alpha-1)*k*msizes))^(alpha-1)*exp(alpha-1-lens/(k*msizes)) },
   "lognorm"={ mu1=pars[1]; sigma=pars[2]
     rselect <- (1/lens)*exp(mu1+log(relsizes)-sigma^2/2-
                    (log(lens)-mu1-log(relsizes))^2/(2*sigma^2)) }
   )

  rselect <- releff*rselect/max(releff)
  rselect <- matrix(rselect,ncol=length(meshsizes))
  return(rselect)
}


#' @title Millar's original gillnet selectivity fitting function
#'
#' @description Function to estimate selectivity parameters from experimental data.
#'    This function is applied within \code{\link{select_Millar}} to derive starting
#'    parameters. \code{\link{select_Millar}} is the recommended function for
#'    selectivity estimation.
#'
#' @param data matrix with the number of individuals caught with each sized mesh
#'    (\code{CatchPerNet_mat}).
#' @param meshsizes vector with meshSizes in increasing order (\code{meshSizes}),
#' @param rtype A character string indicating which method for estimating selection curves
#'    should be used:
#'    \code{"norm.loc"} for a normal curve with common spread,
#'    \code{"norm.sca"} for a normal curve with variable spread,
#'    \code{"lognorm"} for a lognormal curve,
#'    \code{"gamma"} for a gamma curve.
#' @param rel.power A string indicating the relative power of different meshSizes,
#'    must have same length as \code{meshSizes} (Default: \code{rel.power = NULL}).
#' @param plotlens lengths which should be used for graphical output,
#'    for more detailed curves. Default : NULL
#' @param details logical; should details be included in the output?
#'
#' @source https://www.stat.auckland.ac.nz/~millar/selectware/
#'
#' @return list of fitted parameters
#'
#' @importFrom stats as.formula coef confint glm poisson resid
#'
#' @examples
#' data(gillnet)
#'
#' dat <- matrix(c(gillnet$midLengths, gillnet$CatchPerNet_mat),
#'          byrow = FALSE, ncol=(dim(gillnet$CatchPerNet_mat)[2]+1))
#'
#' gillnetfit(data = dat, meshsizes = gillnet$meshSizes)
#'
#'
#' @references
#' Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#'  54(3):471-477
#'
#' @export

gillnetfit <- function(data, meshsizes,
                       rtype="norm.loc",
                       rel.power=NULL,
                       plotlens=NULL,
                       details=FALSE){
  # Adapted R code from Russell Millar (https://www.stat.auckland.ac.nz/~millar/selectware/)


 if(sum(sort(meshsizes)==meshsizes)!=length(meshsizes))
   stop("Mesh sizes must be ascending order")
 lens=rep(data[,1],ncol(data[,-1]))
 msizes=rep(meshsizes,rep(nrow(data),ncol(data[,-1])))
 msize1=msizes[1]
 dat=as.vector(data[,-1])
 var1=lens*msizes; var2=msizes^2; var3=(lens/msizes)^2
 var4=lens/msizes; var5=-log(msizes); var6=log(msizes/msizes[1])
 var7=var6*log(lens) - 0.5*var6*var6; var8=lens*lens
 var9=msizes/lens

 if(is.null(plotlens)) plotlens=data[,1]
 if(is.null(rel.power)) os=0
  else os=rep(log(rel.power),rep(nrow(data),ncol(data[,-1])))
 switch(rtype,
  "norm.loc"={
  if(missing(rel.power))
   fit=glm(dat ~ -1 + var1 + var2 + as.factor(lens),family=poisson)
  else
   fit=glm(dat ~ -1 + var1 + var2 + as.factor(lens) + offset(os),family=poisson)
  x=coef(fit)[c("var1","var2")]
  varx=summary(fit)$cov.unscaled[1:2,1:2]
  k=-2*x[2]/x[1]; sigma=sqrt(-2*x[2]/(x[1]^2))
  vartemp=msm::deltamethod(list(~-2*x2/x1,~sqrt(-2*x2/(x1^2))),x,varx,ses=F)
  pars=c(k,sigma,k*msizes[1],sigma)
  form1=as.formula(sprintf("~x1*%f",msize1)) #Deltamethod quirk
  varpars=msm::deltamethod(list(~x1,~x2,form1,~x2),c(k,sigma),vartemp,ses=F)
  gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars)))
  rownames(gear.pars)=c("k","sigma","mode(mesh1)","std_dev(all meshes)") },
  "norm.sca"={
  if(missing(rel.power))
   fit=glm(dat ~ -1 + var3 + var4 + as.factor(lens),family=poisson)
  else
   fit=glm(dat ~ -1 + var3 + var4 + as.factor(lens) + offset(os),family=poisson)
  x=coef(fit)[c("var3","var4")]
  varx=summary(fit)$cov.unscaled[1:2,1:2]
  k1=-x[2]/(2*x[1]); k2=-1/(2*x[1])
  vartemp=msm::deltamethod(list(~-x2/(2*x1),~-1/(2*x1)),x,varx,ses=F)
  pars=c(k1,k2,k1*msizes[1],sqrt(k2*msizes[1]^2))
  form1=as.formula(sprintf("~x1*%f",msize1)) #Deltamethod quirk
  form2=as.formula(sprintf("~sqrt(x2*%f^2)",msize1)) #Deltamethod quirk
  varpars=msm::deltamethod(list(~x1,~x2,form1,form2),c(k1,k2),vartemp,ses=F)
  gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars)))
  rownames(gear.pars)=c("k1","k2","mode(mesh1)","std_dev(mesh1)") },
  "gamma"={
  if(missing(rel.power))
   fit=glm(dat ~ -1 + var4 + var5 + as.factor(lens),family=poisson)
  else
   fit=glm(dat ~ -1 + var4 + var5 + as.factor(lens) + offset(os),family=poisson)
  x=coef(fit)[c("var4","var5")]
  varx=summary(fit)$cov.unscaled[1:2,1:2]
  alpha=x[2]+1; k=-1/x[1]
  vartemp=msm::deltamethod(list(~x2+1,~-1/x1),x,varx,ses=F)
  pars=c(alpha,k,(alpha-1)*k*msizes[1],sqrt(alpha*(k*msizes[1])^2))
  form1=as.formula(sprintf("~(x1-1)*x2*%f",msize1)) #Deltamethod quirk
  form2=as.formula(sprintf("~sqrt(x1*(x2*%f)^2)",msize1)) #Deltamethod quirk
  varpars=msm::deltamethod(list(~x1,~x2,form1,form2),c(alpha,k),vartemp,ses=F)
  gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars)))
  rownames(gear.pars)=c("alpha","k","mode(mesh1)","std_dev(mesh1)")  },
  "lognorm"={
  if(missing(rel.power))
   fit=glm(dat ~ -1 + var6 + var7 + as.factor(lens),family=poisson)
  else
   fit=glm(dat ~ -1 + var6 + var7 + as.factor(lens) + offset(os),family=poisson)
  x=coef(fit)[c("var6","var7")]
  varx=summary(fit)$cov.unscaled[1:2,1:2]
  mu1=-(x[1]-1)/x[2]; sigma=sqrt(1/x[2])
  vartemp=msm::deltamethod(list(~-(x1-1)/x2,~sqrt(1/x2)),x,varx,ses=F)
  pars=c(mu1,sigma,exp(mu1-sigma^2),sqrt(exp(2*mu1+sigma^2)*(exp(sigma^2)-1)))
  varpars=msm::deltamethod(list(~x1,~x2,~exp(x1-x2^2),
                      ~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1))),c(mu1,sigma),vartemp,ses=F)
  gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars)))
  rownames(gear.pars)=c("mu1(mode log-scale, mesh1)","sigma(std_dev log scale)",
                "mode(mesh1)","std_dev(mesh1)")  },
 stop(paste("\n",rtype, "not recognised, possible curve types are ",
        "\"norm.loc\", \"norm.sca\", \"gamma\", and \"lognorm\"")))
 rselect=rcurves_Millar(rtype,meshsizes,rel.power,pars,plotlens)
 devres=matrix(resid(fit,type="deviance"),nrow(data),ncol(data[,-1]))
 # if(plots[1]) plot.curves(type,plotlens,rselect)
 # if(plots[2]) plot.resids(devres,meshsizes,data[,1])
 g.o.f=c(deviance(fit),sum(resid(fit,type="pearson")^2),fit$df.res,fit$null)
 names(g.o.f)=c("model_dev","Pearson chi-sq","dof","null_dev")
 fit.type=paste(paste(rtype,ifelse(is.null(rel.power),"",": with unequal mesh efficiencies")))
 if(details==FALSE)
  return(list(fit.type=fit.type,gear.pars=gear.pars,fit.stats=g.o.f))
 else
  return(list(
   fit.type=fit.type,gear.pars=gear.pars,fit.stats=g.o.f,devres=devres,rselect=rselect))
}

