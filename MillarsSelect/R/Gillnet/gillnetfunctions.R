#Welcome to:
# 
######  ####### #       #######  #####  #######       ######
#     # #       #       #       #     #    #          #     #
#       #       #       #       #          #          #     #
 #####  #####   #       #####   #          #     ###  ######
      # #       #       #       #          #          #     #
#     # #       #       #       #     #    #          #     #
 #####  ####### ####### #######  #####     #          #     #
#
#R CODE for fitting SELECT models to gillnet data
#*******************************************************************************
#Copyright Russell Millar, *****************************************************
#          Department of Statistics, *******************************************
#          University of Auckland, *********************************************
#          P.O. Box 92019, Auckland.  08 Jan 1998*******************************
#
#Adapted for R from original Splus code:       Russell Millar,  8 Sept 2002
#Updated to produce plots and better output:   Russell Millar,  1 Nov 2003
#Documented (see gillnetfunctions.pdf):        Russell Millar,  1 Nov 2003
#Variances added and return object modified:   Russell Millar, 16 Nov 2009    
#===============================================================================
#THIS SOFTWARE IS FREE TO USE.
#
#THIS SOFTWARE MAY BE REDISTRIBUTED, BUT MUST BE REDISTRIBUTED FREE OF CHARGE.
#
#PLEASE ACKNOWLEDGE THE USE OF THIS SOFTWARE WHERE APPROPRIATE.
#
#NO WARRANTIES ARE GIVEN, AND NO RESPONSIBILITY IS ACCEPTED FOR 
#ANY FAULTS OR CONSEQUENCES ARISING FROM THE USE OF THIS SOFTWARE.
#===============================================================================
#
####################################
#Function for gillnet data analyses#
####################################
gillnetfit=function(data,meshsizes,type="norm.loc",rel=NULL,
                    plots=c(T,T),plotlens=NULL,details=F) {
 require(msm)
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
 if(is.null(rel)) os=0
  else os=rep(log(rel),rep(nrow(data),ncol(data[,-1])))
 switch(type,
  "norm.loc"={
  if(missing(rel))
   fit=glm(dat ~ -1 + var1 + var2 + as.factor(lens),family=poisson)
  else
   fit=glm(dat ~ -1 + var1 + var2 + as.factor(lens) + offset(os),family=poisson)
  x=coef(fit)[c("var1","var2")]
  varx=summary(fit)$cov.unscaled[1:2,1:2]
  k=-2*x[2]/x[1]; sigma=sqrt(-2*x[2]/(x[1]^2))
  vartemp=deltamethod(list(~-2*x2/x1,~sqrt(-2*x2/(x1^2))),x,varx,ses=F)
  pars=c(k,sigma,k*msizes[1],sigma)
  form1=as.formula(sprintf("~x1*%f",msize1)) #Deltamethod quirk
  varpars=deltamethod(list(~x1,~x2,form1,~x2),c(k,sigma),vartemp,ses=F)
  gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars))) 
  rownames(gear.pars)=c("k","sigma","mode(mesh1)","std_dev(all meshes)") },
  "norm.sca"={
  if(missing(rel))
   fit=glm(dat ~ -1 + var3 + var4 + as.factor(lens),family=poisson)
  else
   fit=glm(dat ~ -1 + var3 + var4 + as.factor(lens) + offset(os),family=poisson)
  x=coef(fit)[c("var3","var4")]
  varx=summary(fit)$cov.unscaled[1:2,1:2]
  k1=-x[2]/(2*x[1]); k2=-1/(2*x[1])
  vartemp=deltamethod(list(~-x2/(2*x1),~-1/(2*x1)),x,varx,ses=F)
  pars=c(k1,k2,k1*msizes[1],sqrt(k2*msizes[1]^2))
  form1=as.formula(sprintf("~x1*%f",msize1)) #Deltamethod quirk
  form2=as.formula(sprintf("~sqrt(x2*%f^2)",msize1)) #Deltamethod quirk
  varpars=deltamethod(list(~x1,~x2,form1,form2),c(k1,k2),vartemp,ses=F)
  gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars))) 
  rownames(gear.pars)=c("k1","k2","mode(mesh1)","std_dev(mesh1)") },
  "gamma"={
  if(missing(rel))
   fit=glm(dat ~ -1 + var4 + var5 + as.factor(lens),family=poisson)
  else
   fit=glm(dat ~ -1 + var4 + var5 + as.factor(lens) + offset(os),family=poisson)
  x=coef(fit)[c("var4","var5")] 
  varx=summary(fit)$cov.unscaled[1:2,1:2]
  alpha=x[2]+1; k=-1/x[1]
  vartemp=deltamethod(list(~x2+1,~-1/x1),x,varx,ses=F)
  pars=c(alpha,k,(alpha-1)*k*msizes[1],sqrt(alpha*(k*msizes[1])^2))
  form1=as.formula(sprintf("~(x1-1)*x2*%f",msize1)) #Deltamethod quirk
  form2=as.formula(sprintf("~sqrt(x1*(x2*%f)^2)",msize1)) #Deltamethod quirk
  varpars=deltamethod(list(~x1,~x2,form1,form2),c(alpha,k),vartemp,ses=F)
  gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars)))  
  rownames(gear.pars)=c("alpha","k","mode(mesh1)","std_dev(mesh1)")  },
  "lognorm"={
  if(missing(rel))
   fit=glm(dat ~ -1 + var6 + var7 + as.factor(lens),family=poisson)
  else
   fit=glm(dat ~ -1 + var6 + var7 + as.factor(lens) + offset(os),family=poisson)
  x=coef(fit)[c("var6","var7")] 
  varx=summary(fit)$cov.unscaled[1:2,1:2]
  mu1=-(x[1]-1)/x[2]; sigma=sqrt(1/x[2])
  vartemp=deltamethod(list(~-(x1-1)/x2,~sqrt(1/x2)),x,varx,ses=F)
  pars=c(mu1,sigma,exp(mu1-sigma^2),sqrt(exp(2*mu1+sigma^2)*(exp(sigma^2)-1)))
  varpars=deltamethod(list(~x1,~x2,~exp(x1-x2^2),
                      ~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1))),c(mu1,sigma),vartemp,ses=F)
  gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars)))  
  rownames(gear.pars)=c("mu1(mode log-scale, mesh1)","sigma(std_dev log scale)",
                "mode(mesh1)","std_dev(mesh1)")  },
 stop(paste("\n",type, "not recognised, possible curve types are ", 
        "\"norm.loc\", \"norm.sca\", \"gamma\", and \"lognorm\"")))
 rselect=rcurves(type,meshsizes,rel,pars,plotlens) 
 devres=matrix(resid(fit,type="deviance"),nrow(data),ncol(data[,-1]))
 if(plots[1]) plot.curves(type,plotlens,rselect)
 if(plots[2]) plot.resids(devres,meshsizes,data[,1])
 g.o.f=c(deviance(fit),sum(resid(fit,type="pearson")^2),fit$df.res,fit$null)
 names(g.o.f)=c("model_dev","Pearson chi-sq","dof","null_dev")
 fit.type=paste(paste(type,ifelse(is.null(rel),"",": with unequal mesh efficiencies")))
 if(details==F)
  return(list(fit.type=fit.type,gear.pars=gear.pars,fit.stats=g.o.f))     
 else 
  return(list(
   fit.type=fit.type,gear.pars=gear.pars,fit.stats=g.o.f,devres=devres,rselect=rselect)) }


########################Other functions used by gillnetfit#########################
#Calculate the relative selection curves
rcurves=function(type,meshsizes,rel,pars,plotlens) {
  lens=rep(plotlens,length(meshsizes))
  relsizes=meshsizes/meshsizes[1]
  if(is.null(rel)) releff=1
   else releff=rep(rel,rep(length(plotlens),length(meshsizes)))
  msizes=rep(meshsizes,rep(length(plotlens),length(meshsizes)))
  relsizes=rep(relsizes,rep(length(plotlens),length(relsizes)))
  switch(type,
   "norm.loc"={ k=pars[1]; mean1=pars[3]; sigma=pars[2]
     rselect=exp(-((lens-k*msizes)^2)/(2*sigma^2)) },
   "norm.sca"={ k1=pars[1]; k2=pars[2]
     rselect=exp(-((lens-k1*msizes)^2)/(2*k2*msizes^2)) },
   "gamma"={ alpha=pars[1]; k=pars[2]
     rselect=(lens/((alpha-1)*k*msizes))^(alpha-1)*exp(alpha-1-lens/(k*msizes)) },
   "lognorm"={ mu1=pars[1]; sigma=pars[2]
     rselect=
      (1/lens)*exp(mu1+log(relsizes)-sigma^2/2-
                    (log(lens)-mu1-log(relsizes))^2/(2*sigma^2)) })   
  rselect=releff*rselect/max(releff)
  rselect=matrix(rselect,ncol=length(meshsizes))
  return(rselect)
  }

#Plot the relative selection curves
plot.curves=function(type,plotlens,rselect) {
  plot.title=switch(type,
   "norm.loc"="Normal (common spread)",
   "norm.sca"="Normal",
   "gamma"="Gamma",
   "lognorm"="Log-normal")
  plot.title=paste(plot.title,"retention curve")
  matplot(plotlens,rselect,type="l",lty=1,las=1,ylim=c(0,1),
          xlab="Length (cm)",ylab="Relative retention",main=plot.title) }

#Plot the deviance residuals matrix
plot.resids=
 function(residuals,msizes,lens,cex=1,title="Deviance residuals",...) {
 if(missing(lens)) lens=1:nrow(residuals)
 if(missing(msizes)) msizes=1:ncol(residuals)
 plot(c(min(lens),max(lens)),range(msizes),xlab="Length (cm)",ylab="Mesh size",
     ylim=range(msizes)+(cex/25)*c(-1,1)*(max(msizes)-min(msizes)),
     type="n",main=title,...)
 for(i in 1:nrow(residuals))
  for(j in 1:ncol(residuals))
   points(lens[i],msizes[j],pch=ifelse(residuals[i,j]>0,16,1),
                  cex=3*abs(residuals[i,j])*cex/(abs(max(residuals)))) }
  
#################################EXAMPLES OF USE######################################
#Example of use: From Millar and Holst, ICES J. Mar. Sci. (1997) 54: 471-477
#holt.dat=matrix(scan("holt.dat"),byrow=T,ncol=9)
#meshsizes=c(13.5,14.0,14.8,15.4,15.9,16.6,17.8,19.0)
#The fits implemented by Millar and Holst
#par(mfrow=c(2,2),las=1)
#fit1=gillnetfit(holt.dat,meshsizes,type="norm.loc",plotlens=40:100)
#fit2=gillnetfit(holt.dat,meshsizes,type="norm.sca",plotlens=40:100)
#fit3=gillnetfit(holt.dat,meshsizes,type="gamma",plotlens=40:100)
#fit4=gillnetfit(holt.dat,meshsizes,type="lognorm",plotlens=40:100)

#Specify relative efficiency proportional to meshsize, and produce detailed output
#fit5=gillnetfit(holt.dat,meshsizes,type="norm.loc",rel=meshsizes,details=T)
#Plot retention curves but not residual plot
#gillnetfit(holt.dat,meshsizes,type="norm.loc",plotlens=40:100,plots=c(T,F))

#Example of manual use of the plot.resids function
#par(
#plot.resids(fit5$devres,meshsizes,holt.dat[,1],cex=2,title="")
#plot.resids(fit5$devres,meshsizes,holt.dat[,1],cex=1.5,title="")
#plot.resids(fit5$devres,meshsizes,holt.dat[,1],cex=1,title="")



