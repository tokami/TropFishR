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
#' # Calculation with finer length resolution
#' output <- select_Millar(gillnet, x0 = NULL, rtype = "lognorm")
#' plot(output, plotlens=seq(40,90,0.1))
#'
#' # Use alternate plot settings
#' output <- select_Millar(gillnet, x0 = NULL, rtype = "lognorm")
#' ncolor <- length(output$meshSizes)
#' plot(output, plotlens=seq(40,90,0.1), deviance_plot = FALSE, lty=1, col=rainbow(ncolor))
#' legend("topleft", col=rainbow(8), legend=output$meshSizes, lty=1, title="Mesh size [cm]")
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
#' }
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
    gillnetfit.res <- gillnetfit(data=dat.tmp, meshsizes=meshSizes, type=rtype, rel=rel.power,
      plots=c(FALSE, FALSE), plotlens=NULL, details=TRUE)
    if(rtype %in% c("norm.loc", "norm.sca")){
      x0 <- unname(gillnetfit.res$gear.pars[3:4,1])
    }
    if(rtype %in% c("lognorm")){
      x0 <- unname(gillnetfit.res$gear.pars[1:2,1])
    }
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


