#' @title Millar's selectivity model
#
#' @description  This model estimates the selecitvity of different gears from experimental
#'    catches.
#'
#' @param param A list with following parameters: vector with midlengths of size classes
#'    (\code{$midLengths}), vector with meshSizes in increasing order (\code{$meshSizes}),
#'    and a matrix with the number of individuals caught with each sized mesh
#'    (\code{$CatchPerNet_mat}).
#' @param x0 A string of initial values for the parameters to be optimized over when applying the
#'    function \code{\link[stats]{optim}}.
#' @param rtype A character string indicating which method for estimating selection curves
#'    should be used: \code{"norm.loc"} for a normal curve with common spread,
#'    \code{"norm.sca"} for a normal curve with variable spread,
#'    \code{"lognorm"} for a lognormal curve,
#'    \code{"binorm.sca"} for a bi-normal curve,
#'    \code{"bilognorm"} for a bi-lognormal curve,
#'    \code{"tt.logistic"} for a control and logistic curve
#' @param rel.power A string indicating the relative power of different meshSizes,
#'    must have same length as \code{meshSizes} (Default: \code{rel.power = NULL}).
#'
#'
#' @examples
#' # Trouser trawl net
#' # load data
#' data(haddock)
#'
#' # run model
#' output <- select_Millar(haddock, x0 = c(-10,0.3,0),
#'    rtype = "tt.logistic", rel.power = NULL)
#'
#' # plot results
#' plot(output, plotlens=seq(25,35,0.1))
#' legend(32,0.5,c("Control","Experimental"),lty=1:2,col=1:2,cex=0.8,bg="white")
#'
#'
#' # Gillnet
#' # load data
#' data(gillnet)
#'
#' # run model
#' output <- select_Millar(gillnet, x0 = c(60,4), rel.power = rep(1,8),
#'    rtype = "norm.loc")
#'
#' # investigate results
#' output
#'
#' # plot results
#' plot(output,plotlens=seq(40,90,0.1))
#'
#'
#' # Stacked trammel net
#' #Here, the data are from two expts using different mesh sizes
#' #This analysis assumes common retention curve in both expts.
#' #Note that summary function does not produce residual plot
#' #since lengths are not unique
#' # load data
#' data(trammelnet)
#'
#' # run model
#' output <- select_Millar(trammelnet, x0 = c(25,4),
#'    rtype="norm.loc", rel.power = rep(1,6))
#'
#' # plot results
#' plot(output,plotlens=seq(10,40,0.1))
#'
#'
#' @source https://www.stat.auckland.ac.nz/~millar/selectware/
#'
#' @details Model adapted from the selectivity functions provided by Prof. Dr. Russell Millar
#'   (https://www.stat.auckland.ac.nz/~millar/).
#'
#'
#' @references
#'  Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
#'  Holt, S. J. 1963. A method for determining gear selectivity and its application.
#'  ICNAF Special Publication, 5: 106-115.
#'
#'
#' @export

select_Millar <- function(param,
                          x0,
                          rtype = "norm.loc",
                          rel.power = NULL
                          ){
  res <- param
  meshSizes <- res$meshSizes
  classes <- res$midLengths
  CatchPerNet_mat <- res$CatchPerNet_mat

  if(is.null(CatchPerNet_mat)) CatchPerNet_mat <- cbind(res$numCover,res$numCodend)
  res$CatchPerNet_mat <- CatchPerNet_mat

  #  NetFit  #
  if(sum(sort(meshSizes)==meshSizes)!=length(meshSizes))
    stop("Mesh size must be in ascending order!")

  if(is.null(rel.power)) rel.power=rep(1,length(meshSizes))

  if(!is.null(meshSizes) & ncol(CatchPerNet_mat) != length(meshSizes))
    stop("Number of mesh sizes should be ", ncol(CatchPerNet_mat))

  CountPropns <- CatchPerNet_mat/apply(CatchPerNet_mat,1,sum,na.rm=TRUE)

  fullfit.l <- sum(CatchPerNet_mat*log(CountPropns),na.rm=TRUE)

  r <- rtypes_Millar(rtype) #Get selection curve function

  #  nllhood  #
  nllhood <- function(theta,classes,CatchPerNet_mat,meshSizes,r,rel.power){
    rmatrix <- outer(classes,meshSizes,r,theta)
    rmatrix[is.na(CatchPerNet_mat)] <- NA #No fitted retention for missing meshSizes
    rmatrix <- t(t(rmatrix) * rel.power)
    phi <- rmatrix/apply(rmatrix,1,sum,na.rm=TRUE)
    nll <- - sum(CatchPerNet_mat * log(phi),na.rm=TRUE)
    return(nll)
  }

  res2 <- optim(x0, nllhood, classes = classes, CatchPerNet_mat = CatchPerNet_mat,
               meshSizes = meshSizes, r = r, rel.power = rel.power,
               hessian = T, control = list(trace = F))

  cat("Parameters=",fit$par,",  Deviance=",2*(fullfit.l + fit$value),"\n")

  res3 <- c(res,res2,deviance=deviance,rtype=rtype,rel.power=list(rel.power))  # invisible


  #  Estimates  #

  x <- res3$par; varx <- solve(res3$hess)
  names <- c("Mode(mesh1)","Std dev.(mesh1)")
  switch(res3$rtype,
         "norm.loc"={ pars=x; varpars=varx },
         "norm.sca"={ pars=x; varpars=varx },
         "lognorm"={
           pars=c(exp(x[1]-x[2]^2),sqrt(exp(2*x[1]+x[2]^2)*(exp(x[2]^2)-1)))
           varpars=msm::deltamethod(list(~exp(x1-x2^2),
                                    ~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1))),x,varx,ses=F)},
         "binorm.sca"={
           pars=c(x[1:4],exp(x[5])/(1+exp(x[5])))
           names=c("Mode1(mesh1)","Std dev.1(mesh1)",
                   "Mode2(mesh1)","Std dev.2(mesh1)","P(mode1)")
           varpars=msm::deltamethod(list(~x1,~x2,~x3,~x4,~exp(x5)/(1+exp(x5))),
                               x,varx,ses=F)},
         "bilognorm"={
           pars=c(exp(x[1]-x[2]^2),sqrt(exp(2*x[1]+x[2]^2)*(exp(x[2]^2)-1)),
                  exp(x[3]-x[4]^2),sqrt(exp(2*x[3]+x[4]^2)*(exp(x[4]^2)-1)),
                  exp(x[5])/(1+exp(x[5])))
           names=c("Mode1(mesh1)","Std dev.1(mesh1)",
                   "Mode2(mesh1)","Std dev.2(mesh1)","P(mode1)")
           varpars=msm::deltamethod(
             list(~exp(x1-x2^2),~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1)),
                  ~exp(x3-x4^2),~sqrt(exp(2*x3+x4^2)*(exp(x4^2)-1)),
                  ~exp(x5)/(1+exp(x5))),x,varx,ses=F)},
         "tt.logistic"={
           pars=c(-x[1]/x[2],2*(log(3))/x[2],exp(x[3])/(1+exp(x[3])))
           names=c("L50","SR","p")
           varpars=msm::deltamethod(list(~-x1/x2,~2*log(3)/x2,~exp(x3)/(1+exp(x3))),
                               x,varx,ses=F)},
         stop(paste("\n",res3$rtype, "not recognised, possible curve types are \n",
                    "\"norm.loc\", \"norm.sca\", \"lognorm\" \n",
                    "\"binorm.sca\", \"bilognorm\", and \"tt.logistic\""))
  )#End of switch
  estimates=cbind(pars,sqrt(diag(varpars)))
  colnames(estimates)=c("par","s.e.")
  rownames(estimates)=names
   # return estimates


  #  Summary  #

  nlens=length(classes)
  meshSizes=res3$meshSizes; nmeshes=length(meshSizes)

  O=res3$CatchPerNet_mat; #Matrix of observed classes

  rmatrix=outer(classes,meshSizes,r,res3$par)
  rmatrix[is.na(O)]=NA #No fitted retention for missing meshSizes
  rmatrix=t(t(rmatrix)*res3$rel.power)
  phi=rmatrix/apply(rmatrix,1,sum,na.rm=TRUE)
  E=apply(O,1,sum,na.rm=TRUE)*phi #Matrix of expected classes
  Pearson.resids=(O-E)/sqrt(E)
  Pearson.chisq=sum(Pearson.resids^2,na.rm=TRUE)
  wk=O*log(O/E); wk[is.na(wk)]=0
  Dev.resids=sign(O-E)*sqrt(2*(E-O+wk))
  Deviance=sum(Dev.resids^2,na.rm=TRUE)
  full.l=sum(-O+O*log(O),na.rm=TRUE)
  null.E=matrix(apply(O,1,mean,na.rm=TRUE),nrow=nlens,ncol=nmeshes)
  null.l=sum(-null.E+O*log(null.E),na.rm=TRUE)
  model.l=sum(-E+O*log(E),na.rm=TRUE)
  NonZeroDat=O[apply(O,1,sum,na.rm=TRUE)>0,]
  d.o.f.=nrow(NonZeroDat)*(nmeshes-1)-length(res3$par)-sum(is.na(NonZeroDat))
  out=rbind(null.l,model.l,full.l,Deviance,Pearson.chisq,d.o.f.)
  AreLensUnique=(length(classes)==length(unique(classes)))
  op <- par(mfrow=c(2,1), mar=c(4,4,2,2),oma=c(1,1,1,1))
  if(nmeshes > 2 & AreLensUnique) {
    plot(1,1,xlim=range(classes),xlab="Length (cm)",ylab="Mesh size (cm)",
         ylim=range(meshSizes)+(1/50)*c(-1,1)*(max(meshSizes)-min(meshSizes)), # (cex/50)
         yaxt="n",type="n",main="Deviance residuals")
    axis(2,meshSizes,meshSizes,las=1)
    for(i in 1:nlens)
      for(j in 1:nmeshes)
        points(classes[i],meshSizes[j],pch=ifelse(Dev.resids[i,j]>0,16,1),
               cex=3*abs(Dev.resids[i,j])*1/(abs(max(Dev.resids))))   # cex / (abs...)
    }else
    if(nmeshes == 2) {
      Dev.resids.len=sign(Dev.resids[,2])*sqrt(apply(Dev.resids^2,1,sum))
      plot(classes,Dev.resids.len,type=ifelse(AreLensUnique,"h","p"),las=1,
           main="Deviance residuals",xlab="Length (cm)",ylab="Mesh size (cm)",cex=1)   # cex = cex
      abline(h=0) }
  #return(out)


  #   PlotCurves  #

  r=rtypes_Millar(res3$rtype) #Get selection curve function
  if(is.null(plotlens)) plotlens=classes
  if(is.null(meshSizes)) meshSizes=res3$meshSizes
  plot.title=switch(res3$rtype,
                    "norm.loc"="Normal (common spread)",
                    "norm.sca"="Normal",
                    "lognorm"="Lognormal",
                    "binorm.sca"="Bi-normal",
                    "bilognorm"="Bi-lognormal",
                    "tt.logistic"="Control and logistic","")
  rmatrix=outer(plotlens,meshSizes,r,res3$par)
  rmatrix=t(t(rmatrix)*res3$rel.power)
  rmatrix=rmatrix/max(rmatrix)
  matplot(plotlens,rmatrix,type="l",las=1,ylim=c(0,1),
          xlab="Length (cm)",ylab="Relative retention")
  #abline(h=seq(0,1,0.25),lty=3)
  lenrmatrix=cbind(plotlens,rmatrix)
  colnames(lenrmatrix)=c("Length",meshSizes)
  invisible(lenrmatrix)
  par(op)

  ret <- res3
  class(ret) <- "select_Millar"
  return(ret)
}
