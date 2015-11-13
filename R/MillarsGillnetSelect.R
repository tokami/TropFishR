#' @title Millar'S gillnet selctivity
#
#' @description  This function estimates the selecitvity of a gillnet from an experimental catch with two gillnets with different mesh sizes.
#'
#' @param param A list with following parameters: vector with midlengths of size classes
#'    (\code{$midLengths}), vector with meshSizes in increasing order (\code{$meshSizes}),
#'    and a matrix with the number of individuals caught with each sized mesh
#'    (\code{$CatchPerNet_mat}).
#' @param model A character string indicating which type fo model for the estimation
#'   of the selection curves should be used: \code{"fixed_normal"}, \code{"normal"},
#'   \code{"gamma"}, or \code{"lognormal"}.
#' @param plotlens Vector of lengths for which values of relative retention are required
#' @param rel Vector of realtive efficencies
#'
#' @examples
#' data(data_MillarsGillnet)
#' output <- MillarsGillnetSelect(data_MillarsGillnet, model = "normal_fixed",
#'    plotlens = NULL, rel = NULL)
#' output
#' plot(output)
#'
#' @details This function estimates the fractions retained by each net (SNet1 & SNet2), the
#'   optimum lengths for each net, the selection factor (SF), and the standard deviation
#'   of the factor (stand.dev).
#'   Assumptions of this method are, that (i) the optimum length Lm is proportional to the mesh
#'   size (Lm = SF * m), (ii) the selection curves are normally distributed with a common
#'   standard deviation, (iii) the nets have the same fishing power (same dimensions and material).
#'   Requirements for the experimental set-up are: selection curves corresponding to the two
#'   mesh sizes have to overlap, and the nets have to be set in the same area, during the
#'   same time.
#'
#' @references
#'  Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
#'  Holt, S. J. 1963. A method for determining gear selectivity and its application.
#'  ICNAF Special Publication, 5: 106–115.
#'
#'
#' @export

MillarsGillnetSelect <- function(param, model = "normal_fixed", rel = NULL,
                           plots=c(T,T),plotlens=NULL){

  res <- param
  classes <- as.character(res$midLengths)
  catch_mat <- res$CatchPerNet_mat
  meshSizes <- res$meshSizes

  if(sum(sort(meshSizes)==meshSizes)!=length(meshSizes))
    stop("Mesh sizes must be ascending order")

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
  classes.num <- as.numeric(classes.num[,1])

  lens <- rep(classes.num, ncol(catch_mat))
  msizes <- rep(meshSizes, rep(nrow(catch_mat), ncol(catch_mat)))
  msize1 <- msizes[1]
  dat <- as.vector(catch_mat)

  # variables
  var1 <- lens * msizes
  var2 <- msizes ^ 2
  var3 <- (lens / msizes) ^ 2
  var4 <- lens / msizes
  var5 <- -log(msizes)
  var6 <- log(msizes / msizes[1])
  var7 <- var6 * log(lens) - 0.5 * var6 * var6
  var8 <- lens * lens
  var9 <- msizes / lens


  #necessary??????
  if(is.null(plotlens)) plotlens <- classes.num
  if(is.null(rel)){
    os=0
  }else os=rep(log(rel),rep(nrow(catch_mat),ncol(catch_mat)))


  switch(model,

         "normal_fixed"={
           if(is.null(rel)){
             fit <- glm(dat ~ -1 + var1 + var2 + as.factor(lens), family = poisson)
           } else{
             fit <- glm(dat ~ -1 + var1 + var2 + as.factor(lens) +
                          offset(os), family = poisson)
             }

           x <- coef(fit)[c("var1","var2")]

           varx <- summary(fit)$cov.unscaled[1:2,1:2]

           k <- -2*x[2]/x[1]
           sigma <- sqrt(-2*x[2]/(x[1]^2))

           vartemp <- msm::deltamethod(list(~-2*x2/x1,~sqrt(-2*x2/(x1^2))),x,varx,ses=F)
           pars <- c(k,sigma,k*msizes[1],sigma)
           form1 <- as.formula(sprintf("~x1*%f",msize1)) #Deltamethod quirk
           varpars <- msm::deltamethod(list(~x1,~x2,form1,~x2),c(k,sigma),vartemp,ses=F)
           gear.pars <- cbind(estimate=pars,s.e.=sqrt(diag(varpars)))
           rownames(gear.pars)=c("k","sigma","mode(mesh1)","std_dev(all meshes)")
           },

         "normal"={

           if(missing(rel)){
             fit <- glm(dat ~ -1 + var3 + var4 + as.factor(lens), family = poisson)
           }else {
             fit <- glm(dat ~ -1 + var3 + var4 + as.factor(lens) +
                          offset(os), family = poisson)
           }

           x <- coef(fit)[c("var3","var4")]

           varx <- summary(fit)$cov.unscaled[1:2,1:2]
           k1 <- -x[2]/(2*x[1])
           k2 <- -1/(2*x[1])

           vartemp <- msm::deltamethod(list(~-x2/(2*x1),~-1/(2*x1)),x,varx,ses=F)
           pars <- c(k1,k2,k1*msizes[1],sqrt(k2*msizes[1]^2))
           form1 <- as.formula(sprintf("~x1*%f",msize1)) #Deltamethod quirk
           form2 <- as.formula(sprintf("~sqrt(x2*%f^2)",msize1)) #Deltamethod quirk
           varpars <- msm::deltamethod(list(~x1,~x2,form1,form2),c(k1,k2),vartemp,ses=F)
           gear.pars <- cbind(estimate = pars, s.e. = sqrt(diag(varpars)))
           rownames(gear.pars) <- c("k1","k2","mode(mesh1)","std_dev(mesh1)")
           },

         "gamma"={
           if(missing(rel)){
             fit <- glm(dat ~ -1 + var4 + var5 + as.factor(lens),family=poisson)
           }else{
             fit <- glm(dat ~ -1 + var4 + var5 + as.factor(lens) +
                          offset(os), family = poisson)
           }

           x <- coef(fit)[c("var4","var5")]

           varx <- summary(fit)$cov.unscaled[1:2,1:2]

           alpha <- x[2]+1
           k <- -1/x[1]

           vartemp <- msm::deltamethod(list(~x2+1,~-1/x1),x,varx,ses=F)

           pars <- c(alpha,k,(alpha-1)*k*msizes[1],sqrt(alpha*(k*msizes[1])^2))
           form1 <- as.formula(sprintf("~(x1-1)*x2*%f",msize1)) #Deltamethod quirk
           form2 <- as.formula(sprintf("~sqrt(x1*(x2*%f)^2)",msize1)) #Deltamethod quirk
           varpars <- msm::deltamethod(list(~x1,~x2,form1,form2),c(alpha,k),vartemp,ses=F)
           gear.pars <- cbind(estimate = pars,s.e. = sqrt(diag(varpars)))
           rownames(gear.pars) <- c("alpha","k","mode(mesh1)","std_dev(mesh1)")  },

         "lognormal"={
           if(missing(rel)){
             fit <- glm(dat ~ -1 + var6 + var7 + as.factor(lens),family=poisson)
           }else{
             fit <- glm(dat ~ -1 + var6 + var7 + as.factor(lens) +
                          offset(os), family = poisson)
           }

           x <- coef(fit)[c("var6","var7")]

           varx <- summary(fit)$cov.unscaled[1:2,1:2]
           mu1 <- -(x[1]-1)/x[2]
           sigma <- sqrt(1/x[2])

           vartemp <- msm::deltamethod(list(~-(x1-1)/x2,~sqrt(1/x2)),x,varx,ses=F)
           pars <- c(mu1,sigma,exp(mu1-sigma^2),sqrt(exp(2*mu1+sigma^2)*(exp(sigma^2)-1)))
           varpars <- msm::deltamethod(list(~x1,~x2,~exp(x1-x2^2),
                                    ~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1))),
                                  c(mu1,sigma),vartemp,ses=F)
           gear.pars <- cbind(estimate=pars,s.e.=sqrt(diag(varpars)))
           rownames(gear.pars)=c("mu1(mode log-scale, mesh1)","sigma(std_dev log scale)",
                                 "mode(mesh1)","std_dev(mesh1)")  },

         stop(paste("\n",model, "not recognised, possible curve types are ",
                    "\"normal_fixed\", \"normal\", \"gamma\", and \"lognormal\"")))

  ### other function rcurves
  rsel.lens <- rep(plotlens,length(meshSizes))

  relSizes <- meshSizes / meshSizes[1]

  if(is.null(rel)) { releff <- 1
  }else releff <- rep(rel,rep(length(plotlens),length(meshSizes)))

  msizes <- rep(meshSizes,rep(length(plotlens),length(meshSizes)))
  relSizes <- rep(relSizes,rep(length(plotlens),length(relSizes)))


  switch(model,
         "normal_fixed"={ k=pars[1]; mean1=pars[3]; sigma=pars[2]
         rselect=exp(-((rsel.lens-k*msizes)^2)/(2*sigma^2)) },
         "normal"={ k1=pars[1]; k2=pars[2]
         rselect=exp(-((rsel.lens-k1*msizes)^2)/(2*k2*msizes^2)) },
         "gamma"={ alpha=pars[1]; k=pars[2]
         rselect=(rsel.lens/((alpha-1)*k*msizes))^(alpha-1)*exp(alpha-1-rsel.lens/(k*msizes)) },
         "lognormal"={ mu1=pars[1]; sigma=pars[2]
         rselect=
           (1/rsel.lens)*exp(mu1+log(relSizes)-sigma^2/2-
                          (log(rsel.lens)-mu1-log(relSizes))^2/(2*sigma^2)) })

  rselect <- releff * rselect / max(releff)
  rselect <- matrix(rselect, ncol=length(meshSizes))
  ###

  devres <- matrix(resid(fit,type="deviance"),nrow(catch_mat),ncol(catch_mat))

  #Plot the relative selection curves
  plot.title <- switch(model,
                       "normal_fixed"="Normal (common spread)",
                       "normal"="Normal",
                       "gamma"="Gamma",
                       "lognormal"="Log-normal")
  plot.title <- paste(plot.title,"retention curve")

  matplot(plotlens,rselect,type="l",lty=1,las=1,ylim=c(0,1),
          xlab="Length (cm)",ylab="Relative retention",main=plot.title)

  plot(c(min(classes.num),max(classes.num)),range(meshSizes), xlab="Length (cm)",
       ylab="Mesh size",
       ylim=range(meshSizes)+(1/25)*c(-1,1)*(max(meshSizes)-min(meshSizes)),
       type="n", main = "Deviance residuals")

  for(i in 1:nrow(devres)){
    for(j in 1:ncol(devres)){
      points(classes.num[i],meshSizes[j],pch=ifelse(devres[i,j] > 0,16,1),
             cex=3*abs(devres[i,j])*1/(abs(max(devres))))
    }
  }

  g.o.f <- c(deviance(fit),sum(resid(fit,type="pearson")^2),fit$df.res,fit$null)
  names(g.o.f) <- c("model_dev","Pearson chi-sq","dof","null_dev")
  fit.type <- paste(paste(model,ifelse(is.null(rel),"",": with unequal mesh efficiencies")))


  res2 <- list(model = model, rselect = rselect, plotlens = plotlens, msizes = msizes,
               devres = devres, gear.pars = gear.pars, fit.stats = g.o.f)
  ret <- c(res,res2)
  class(ret) <- "MillarsGillnet"
  return(ret)
}
