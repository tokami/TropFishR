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
#R CODE for fitting SELECT models to covered codend and alternate hauls data
#*******************************************************************************
#Copyright Russell Millar, *****************************************************
#          Department of Statistics, *******************************************
#          University of Auckland, *********************************************
#          P.O. Box 92019, Auckland.  08 Jan 1998*******************************
#
#Adapted for R from original Splus code:    Russell Millar, 17 April 2002
#Updated to pass data matrix           :    Russell Millar,  9 July 2003
#Rep.ttfit added                       :    Russell Millar, 27 May 2004

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
#****The core functions are ttfit() for alternate haul or trouser trawl data****
#****and ccfit() for covered codend data****************************************
#****and Rep.ttfit() for calculation of replicate estimate of overdispersion****
#*****from individual haul alternate or trouser trawl data**********************
#
#***********The sample program code in haddock.R, developed 9 July 2003*********
#*******************************************************************************
#*******************************************************************************
#The data are assumed to be in a matrix with 3 columns:
#Column 1:  Midpoint of lengthclass
#Column 3:  Numbers in control codend (trouser trawl) or cover 
#Column 2:  Numbers in experimental codend
#
#EXAMPLES:  If the data matrix is called  catch  then:
#For a logistic fit to covered codend    ccfit(catch)
#For a Richards fit to covered codend    ccfit(catch,type="rich")
#For a logistic fit to trouser trawl     ttfit(catch)
#For a Richards fit to trouser trawl     ttfit(catch,type="rich")
#To fix split parameter to 0.5, say      ttfit(catch,psplit=0.5)
#
#Plots of fit and deviance resids are produced unless  plots=F  is specified.
#******************************************************************************

######################################
#Function for alternate haul analyses#
######################################
ttfit=function(catch=catchdat,type="logit",probs=c(0.25,0.5,0.75),
             psplit=NULL,x0=c(-10,0.3,0.5),delta=1.0,suff.big=3,nullfit=F,
             plots=T,cex=0.8,mkh=0.07,error.bars=FALSE,plotlens=NULL,details=F){
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]; ntotal=nfine+nwide
  nobs=length(ntotal[ntotal>0])
  if(!is.null(psplit)) cat("\n"," Fixed split, p= ",psplit)
  fullfithood=sum ( nwide*log(ifelse(nwide>0.001,nwide/ntotal,1)) +
          nfine*log(ifelse(nfine>0.001,nfine/ntotal,1)) )
  nullfithood=sum(nwide)*log(sum(nwide)/sum(ntotal)) +
               sum(nfine)*log(sum(nfine)/sum(ntotal))
  #cat("\n"," Likelihood of null model is ",format(nullfithood))
  #cat("\n"," Likelihood of full model is ",format(fullfithood))

  if(plots) {
   propn=ifelse(ntotal>0.001,nwide/ntotal,2)
   xyticks=c(length(lenclass)-1,10,7)
   uniquelens=sort(unique(lenclass))
   AreLensUnique=(length(lenclass)==length(uniquelens))
   if(!error.bars) 
    plot(lenclass[propn!=2], propn[propn!=2], pch = 5, mkh=mkh,
           type=ifelse(AreLensUnique,"b","p"),
           lab = xyticks, xlab = "", ylab = "", xlim = c(lenclass[1],
           lenclass[length(lenclass)]), ylim = c(0,1),cex=cex)
   else {
    lower.bnds=pmax(propn[propn!=2] + qnorm(0.1/2)*0.5/sqrt(ntotal)[propn!=2],0)
    upper.bnds=pmin(propn[propn!=2] - qnorm(0.1/2)*0.5/sqrt(ntotal)[propn!=2],1)
    error.bar(lenclass[propn!=2],propn[propn!=2],lower.bnds,upper.bnds,incr=F,
          pch = 5,lab = xyticks,
          xlab = "Length (cm)",ylab = "Propn retained in codend",mkh=0.07,
          xlim = c(lenclass[1], lenclass[length(lenclass)]), ylim = c(0,1)) }
   title(xlab = "Length (cm)", ylab = "Propn in large mesh codend",
                 main="Proportion of catch in large mesh codend",cex=cex)}

  if(nullfit) 
    {
    if(!is.null(psplit)) cselect=rep(psplit,length(lenclass))
      else cselect=rep(sum(nwide)/sum(ntotal),length(lenclass))
    }
  else
  {
  cat("\n","***NOTE: warning messages may occur as normal part of optimization***","\n") 
  if(type=="logit") 
   {
   if(!is.null(psplit)) {
    Tfit=nlm(hood2par,x0[1:2],iterlim=200,catch=catch,psplit=psplit); 
    Pars=Tfit$est;
    Tcov=cov2par(Pars,catch,type="logit",p=psplit)
    select=lselect(Pars,lenclass)
    if(is.null(plotlens)) r=select else r=lselect(Pars,plotlens)  
    cselect=psplit*select/(psplit*select + 1-psplit)
    Tlens=retentionlens(Pars,cov=Tcov$covar,probs=probs) 
    p=cbind(psplit,NA) }
   else {
    Tfit=nlm(hood3par,x0,iterlim=300,catch=catch); Pars=Tfit$est;
    Tcov=cov3par(Pars,catch,type="logit") 
    select=lselect(Pars,lenclass) 
    if(is.null(plotlens)) r=select else r=lselect(Pars,plotlens)  
    cselect=Pars[3]*select/(Pars[3]*select + 1-Pars[3])
    Tlens=retentionlens(Pars,cov=Tcov$covar[1:2,1:2],probs=probs) 
    p=cbind(Pars[3],sqrt(Tcov$covar[3,3])) }
   }
  if(type=="rich") 
   {
   if(!is.null(psplit)) {
    Tfit=nlm(richhood2par,c(x0[1:2],delta),iterlim=300,catch=catch,psplit=psplit)
    Pars=Tfit$est
    cat("\n"," Likelihood of fitted model is ",
             format(-richhood2par(Pars,catch,psplit)),"\n")
    select=lselect(Pars[1:2],lenclass)^(1/Pars[3])
    if(is.null(plotlens)) r=select else r=lselect(Pars[1:2],plotlens)^(1/Pars[3])
    cselect=psplit*select/(psplit*select + 1-psplit)
    Tcov=covrich(Pars,catch,npars=2,p=psplit)
    Tlens=retentionlens(Pars,cov=Tcov$covar[1:3,1:3],type="rich",probs=probs)
    p=cbind(psplit,NA)  }
   else {
    Tfit=nlm(richhood,c(x0[1:2],delta,x0[3]),iterlim=500,catch=catch); Pars=Tfit$est
    cat("\n"," Likelihood of fitted model is ",
             format(-richhood(Pars,catch)),"\n")
    select=lselect(Pars[1:2],lenclass)^(1/Pars[3])
    if(is.null(plotlens)) r=select else r=lselect(Pars[1:2],plotlens)^(1/Pars[3])
    cselect=Pars[4]*select/(1-Pars[4] + Pars[4]*select) 
    Tcov=covrich(Pars,catch)
    Tlens=retentionlens(Pars,cov=Tcov$covar[1:3,1:3],type="rich",probs=probs)
    p=cbind(Pars[4],sqrt(Tcov$covar[4,4])) }
     }
   modelhood=-Tfit$min
   lhoods=cbind(c(modelhood,nullfithood,fullfithood),nobs-c(length(Pars),1,nobs))
   }
   Tdevres=devres(nwide,ntotal*cselect,ntotal,suff.big)
   
   if(plots){
   lines(lenclass[order(lenclass)],cselect[order(lenclass)],type="l",lty=2)
   plot(lenclass,Tdevres$devres,type=ifelse(AreLensUnique,"h","p"),
          lab=xyticks,xlab="",ylab="",cex=cex)
   abline(h=0)
   title(xlab="Length (cm)",ylab="Deviance residual",
              main="Deviance residuals",cex=cex) }
  if(nullfit)
   {
    if(!is.null(psplit)) list(p=psplit)
      else list(p=sum(nwide)/sum(ntotal))
    }
  else
  if(details)
   list(converged=Tfit$code,x=Pars,l=lhoods,lens=Tlens$lens,sr=Tlens$sr,p=p,
       xcovar=Tcov$covar,lensr.covar=Tlens$covar,
       r=r,devres=Tdevres$devres,suff.dat=Tdevres$suff.dat) 
  else
   list(converged=Tfit$code,x=Pars,l=lhoods,lens=Tlens$lens,sr=Tlens$sr,p=p)}

######################################
#Function for covered codend analyses#
######################################
ccfit=function(catch=catchdat,type="logit",probs=c(0.25,0.5,0.75),x0=c(-10,0.3),
               delta=1.0,plots=T,suff.big=3,error.bars=F,plotlens=NULL,details=F) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]; ntotal=nfine+nwide
  nobs=length(ntotal[ntotal>0])
  fullfithood=sum ( nwide*log(ifelse(nwide>0.001,nwide/ntotal,1)) +
          nfine*log(ifelse(nfine>0.001,nfine/ntotal,1)) )
  nullfithood=sum(nwide)*log(sum(nwide)/sum(ntotal)) +
              sum(nfine)*log(sum(nfine)/sum(ntotal))
  #cat("\n"," Likelihood of full model is ",format(fullfithood))
  #cat("\n"," Likelihood of null model is ",format(nullfithood))
  if(plots) 
   {
   propn=ifelse(ntotal>0.001,nwide/ntotal,2) 
   xyticks=c(length(lenclass)-1,10,7)
   uniquelens=sort(unique(lenclass))
   AreLensUnique=(length(lenclass)==length(uniquelens))
   if(!error.bars) {
    plot(lenclass[propn!=2], propn[propn!=2],pch = 5,lab = xyticks,
          type=ifelse(AreLensUnique,"b","p"), 
          xlab = "Length (cm)", ylab = "Propn retained in codend",mkh=0.07,las=1, 
          xlim = c(lenclass[1], lenclass[length(lenclass)]), ylim = c(0,1)) }
   else {
    lower.bnds=pmax(propn[propn!=2] + qnorm(0.05)*0.5/sqrt(ntotal)[propn!=2],0)
    upper.bnds=pmin(propn[propn!=2] - qnorm(0.05)*0.5/sqrt(ntotal)[propn!=2],1)
    error.bar(lenclass[propn!=2],propn[propn!=2],lower.bnds,upper.bnds,incr=F,
          pch = 5,lab = xyticks,
          xlab = "Length (cm)",ylab = "Propn retained in codend",mkh=0.07,
          xlim = c(lenclass[1], lenclass[length(lenclass)]), ylim = c(0,1)) }
    title(main="Proportion in codend") 
   }
  if(type=="logit")
   {
    Tfit=nlm(cchood,x0,iterlim=100,catch=catch); Pars=Tfit$est;
    cat("\n"," Log-likelihood of logistic model is ",
          format(-cchood(Pars,catch)))
    Tcov=cccov(Pars,catch)
    Tlens=retentionlens(Pars,Tcov$covar)
    select=lselect(Pars,lenclass)
    if(is.null(plotlens)) r=select else r=lselect(Pars,plotlens)
   }
  if(type=="rich")
   {
    Tfit=nlm(ccrichhood,c(x0,delta),iterlim=200,catch=catch); Pars=Tfit$est;
    cat("\n"," Log-likelihood of Richard's model is ",
          format(-ccrichhood(Pars,catch)))
    Tcov=cccovrich(Pars,catch)
    Tlens=retentionlens(Pars,Tcov$covar,type="rich")
    select=lselect(Pars[1:2],lenclass)^(1/Pars[3])
    if(is.null(plotlens)) r=select else r=lselect(Pars[1:2],plotlens)^(1/Pars[3])
   }
  Tdevres=devres(nwide,ntotal*select,ntotal,suff.big) 
  modelhood=-Tfit$min
  lhoods=cbind(c(modelhood,nullfithood,fullfithood),nobs-c(length(Pars),1,nobs))
  if(plots) {
    lines(lenclass[order(lenclass)],select[order(lenclass)],type="l",lty=2)
    abline(h=c(0.25,0.5,0.75),lty=3)
    plot(lenclass,Tdevres$devres,type=ifelse(AreLensUnique,"h","p"),
          xlab="Length (cm)",lab=xyticks,las=1,ylab="Deviance residual")
    title("Residual plot") 
    abline(h=0) }
  if(details)
   list(converged=Tfit$code,x=Pars,l=lhoods,lens=Tlens$lens,sr=Tlens$sr,
       xcovar=Tcov$covar,lensr.covar=Tlens$covar,r=r,devres=Tdevres$devres,
       suff.dat=Tdevres$suff.dat) 
  else
   list(converged=Tfit$code,x=Pars,l=lhoods,lens=Tlens$lens,sr=Tlens$sr) }

######################################################################
#Function for REP calculation for individual hauls alternate haul data
######################################################################
#Replicate tows fits with individually estimated p's (if ind.psplits=T)
#Also fit null model (no selectivity)
Rep.ttfit=function(catch=catchdat,ntows,numlens,type="logit",x0=c(-10,0.3,0.5),
                   delta=1.0,ind.psplit=T,suff.big=3,details=F) {
  cat("\n","COMBINED TOWS FIT:")
  lenclass=catch[1:numlens,1]
  nobs=nrow(catch)
  if(ntows*numlens != nobs) {
    cat("\n","ERROR: #lengths or #tows incorrectly specified")
    cat("\n","All tows must use same range of lengthclasses")
    return() }
  if(length(unique(catch[,1])) != numlens)
    cat("\n","***WARNING: #unique lengthclasses differs from numlens***")
  nwide=catch[,2]
  nfine=catch[,3]
  propninwide=sum(nwide)/sum(nwide+nfine)
  lenorder=order(catch[,1]) #Use below to get plot right
  combfit=ttfit(catch[lenorder,],x0=x0,delta=delta,type=type,suff.big=suff.big)
  x=combfit$x
  cat("\n"," Combined haul parameters:") 
  if(type=="logit") {
    r=lselect(x[1:2],lenclass)
    cat("\n"," a, b and p:",x) }
  if(type=="rich") {
    r=lselect(x[1:2],lenclass)^(1/x[3])
    cat("\n"," a, b, delta, and p:",x) }
  
  if(type=="rich") cat("\n","***Method is approximate for non-logistic curves***")
  if(ind.psplit)
    cat("\n","FIT OF COMBINED-HAULS CURVE TO REPLICATE TOWS WITH INDIVIDUAL SPLITS")
  else cat("\n","FIT OF COMBINED-HAULS CURVE TO REPLICATE TOWS WITH COMMON SPLIT")
  
  suff.dat=rep(T,nobs)
  psplits=rep(0,ntows)

  for(nonsel in c(F,T)){
  if(nonsel) {
    if(ind.psplit) 
      cat("\n","NONSELECTIVE FIT WITH INDIVIDUAL TOW SPLITS")
    else cat("\n","NONSELECTIVE FIT USING A COMMON SPLIT")
    r=1
    }

  Pearson=0; DevRes=0; llhood=0; LensUsed=0; TowsUsed=0
  #Loop over tow
  for(tow in 1:ntows) { 
    indices=(tow-1)*numlens+(1:numlens)
    nwide=catch[indices,2]
    nfine=catch[indices,3]
    n=nwide+nfine; p=ifelse(n>0,nwide/n,0); y=nwide
  #Calculate p and phi
    if(sum(n)>0) {
      if(ind.psplit) psplit=sum(nwide)/sum(nwide+r*nfine)
        else if(nonsel) psplit=propninwide
          else psplit=x[3]
      phi=psplit*r/((1-psplit)+psplit*r)
      yhat=n*phi
      PearsonStat=(y-yhat)^2/(n*phi*(1-phi))
      l=y*log(phi)+(n-y)*log(1-phi)
      dev=ifelse(p>0,y*(log(p)-log(phi)),0)+
        ifelse(p<1,(n-y)*(log(1-p)-log(1-phi)),0)
      if(!nonsel) suff.dat[indices]=(yhat>suff.big & (n-yhat)>suff.big)
      suff=suff.dat[indices]
      if(sum(suff)>0) {
        TowsUsed=TowsUsed+1
        llhood=llhood+sum(l[suff])
        Pearson=Pearson+sum(PearsonStat[suff])
        DevRes=DevRes+sum(2*dev[suff])
        LensUsed=LensUsed+sum(suff) }
    }
    if(!nonsel) psplits[tow]=psplit
  }
  cat("\n\n Using suff=",suff.big,": #lens used=",LensUsed,"#tows used=",TowsUsed,"and")
  cat("\n  Pearson Chisq=",round(Pearson,4),", Model dev=",round(DevRes,4),
              ", l=",round(llhood,4),"\n\n")
  if(ind.psplit) dof=LensUsed-TowsUsed-2
        else if(nonsel) dof=LensUsed-1
          else dof=LensUsed-3
  if(!nonsel) {
    DevCF=DevRes/dof
    cat(" Using deviance, REP factor is ",
          round(DevRes,4),"/",dof,"=",round(DevCF,4),
        ", P-value=",1-pchisq(DevRes,dof),"\n")
    cat(" With correction applied: ","\n")    
    cat(" L25:",combfit$lens[1,]*c(1,sqrt(DevCF)),
      "\n"," L50:",combfit$lens[2,]*c(1,sqrt(DevCF)),
      "\n"," L75:",combfit$lens[3,]*c(1,sqrt(DevCF)),
      "\n","  SR:",combfit$sr*c(1,sqrt(DevCF)),
      "\n","   p:",combfit$p*c(1,sqrt(DevCF)),"\n")
  }
}
if(details)
  list(DevCF=DevCF,psplits=psplits)
else
  list(DevCF=DevCF)
}

###############################################################################
##########END OF USER DESIGNED FUNCTIONS#######################################
###############################################################################

#******************************************************************************
#These next functions are called by ttfit(), ccfit(), and Rep.ttfit() 
#and would not normally be called by the user.
#******************************************************************************

lselect=function(x,lenclass) {
   expo=exp(pmin(500,x[1]+x[2]*lenclass))
   expo/(1+expo) }

cchood=function(x,catch) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  select=lselect(x,lenclass)
  -sum(nwide*log(select) + nfine*log(1-select)) }

ccrichhood=function(x,catch) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  select=lselect(x[1:2],lenclass)^(1/x[3])
  -sum( nwide*ifelse(select>0,log(select),-1e+06) +
      nfine*ifelse(select<1,log(1-select),-1e+06) ) }

hood2par=function(x,catch,psplit) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  expo=exp(x[1]+x[2]*lenclass)
  cselect=psplit*expo/( 1-psplit + expo)
  -sum( nwide*log(cselect) + nfine*log(1-cselect) ) }

hood3par=function(x,catch) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  expo=exp(x[1]+x[2]*lenclass)
  cselect=x[3]*expo/( (1-x[3]) + expo)
  -sum( nwide*log(cselect) + nfine*log(1-cselect) ) }

richhood=function(x,catch) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  select=lselect(x[1:2],lenclass)^(1/x[3])
  cselect=x[4]*select/(1-x[4] + x[4]*select)
  -sum(nwide*log(cselect) + nfine*log(1-cselect)) }

richhood2par=function(x,catch,psplit) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  select=lselect(x[1:2],lenclass)^(1/x[3])
  cselect=psplit*select/(1-psplit + psplit*select)
  -sum(nwide*log(cselect) + nfine*log(1-cselect)) }

#Returns the Pearson and deviance residuals
devres=function(y,yhat,n,suff.big=3) {
  if( any((n==0)&(y>0 | yhat>0)) ) stop("Wrong data in function devres")
  if( any((yhat==0 & y >0) | (yhat==n & y<n)) )
    stop("Impossibility in function devres")
  p=ifelse(n>0,y/n,0.5); phat=ifelse(n>0,yhat/n,0.5)
  sign=ifelse(y>=yhat,1,-1)
  Pearson=(y-yhat)/ifelse(n*phat*(1-phat)>0,sqrt(n*phat*(1-phat)),1)
  l=ifelse(y>0,y*(log(p)-log(phat)),0) +
     ifelse(y<n,(n-y)*(log(1-p)-log(1-phat)),0)
  cat("\n"," Pearson Chisq=",round(sum(Pearson^2),4),
      ", Dev=",round(2*sum(l),4)," #lens used=",sum(n>0),sep="")
  suff.dat=(yhat>suff.big & (n-yhat)>suff.big)
  cat("\n"," Pearson Chisq=",round(sum(Pearson[suff.dat]^2),4),
      ", Dev=",round(2*sum(l[suff.dat]),4),
      " #lens used=",sum(suff.dat)," (Expected count >",suff.big,")","\n",sep="")
  list(Pearson=Pearson,devres=sign*sqrt(2*l),suff.dat=suff.dat) }

retentionlens=function(x,covar,type="logit",probs=c(0.25,0.5,0.75)) {
  np=length(probs)
  if(type=="logit") {
    rlens=( log(probs/(1-probs)) -x[1] ) / x[2]
    srange=2*log(3)/x[2]
    if(!missing(covar)) {
      derivs=matrix(0,nrow=2,ncol=np+1)
      for(i in 1:np) derivs[,i]=c(-1/x[2],-rlens[i]/x[2])
      derivs[,np+1]=c(0,-srange/x[2])
      rlencovar= t(derivs) %*% covar %*% derivs 
      lens=matrix(c(rlens,sqrt(diag(rlencovar)[1:np])),nrow=np,byrow=F)
      sr=c(srange,sqrt(rlencovar[np+1,np+1]))
      return(list(lens=lens,sr=sr,covar=rlencovar)) } }
  if(type=="rich") {
    work=(probs^x[3])/(1-probs^x[3]); rlens=(log(work)-x[1])/x[2]
    worksr=(c(0.25,0.75)^x[3])/(1-c(0.25,0.75)^x[3])
    srange=(log(worksr[2])-log(worksr[1]))/x[2]
    if(!missing(covar)) {
     derivs=matrix(0,nrow=3,ncol=np+1)
     for(i in 1:np)
      derivs[,i]=c(-1/x[2],-rlens[i]/x[2],log(probs[i])/(x[2]*(1-probs[i]^x[3])))
      derivs[,np+1]=c(0,-srange/x[2],
                     (log(0.75)/(1-0.75^x[3])-log(0.25)/(1-0.25^x[3]))/x[2])
     rlencovar= t(derivs) %*% covar %*% derivs
     lens=matrix(c(rlens,sqrt(diag(rlencovar)[1:np])),nrow=np,byrow=F)
     sr=c(srange,sqrt(rlencovar[np+1,np+1]))
     return(list(lens=lens,sr=sr,covar=rlencovar)) } }
   }

#**************************************************************************
#***************************COVARIANCE FUNCTIONS***************************
#**************************************************************************
#***If problems occur with trouser trawl cov functions then numerical******
#***accuracy enhancements used in cccovrich will need to be added**********

#Information matrix and diagnostic calculations for fixed split fit.
cov2par=function(x,catch,type="logit",p=0.5) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  ntotal=nfine+nwide
  cat("\n","cov2par: Fixed split p= ",p)
  eta=x[1] + x[2]*lenclass
  if(type=="logit") {
    mu=p*exp(eta)/(1-p+exp(eta))
    dmudeta=p*(1-p)*exp(eta)/( (1-p+exp(eta))^2 ) }
  if(type=="cl") {
    dblexpo=exp(-exp(eta))
    mu=1-(1-p)/(1-p*dblexpo)
    dmudeta=(1-p)*p*exp(eta)*dblexpo/ ( (1-p*dblexpo)^2 ) }
  if(type=="logit" | type=="cl") {
    info=matrix(0,nrow=2,ncol=2)
    info[1,1]=sum( ntotal*dmudeta^2 / (mu*(1-mu)) )
    info[1,2]=sum( ntotal*lenclass*dmudeta^2 / (mu*(1-mu)) )
    info[2,1]=info[1,2]
    info[2,2]=sum( ntotal*(lenclass*dmudeta)^2 / (mu*(1-mu)) )
    covar=solve(info) 
    list(covar=covar)  }
  else return(NA)  }

#Information matrix and diagnostic calculations for estimated p fit.
cov3par=function(x,catch,type="logit") {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  ntotal=nfine+nwide
  eta=x[1] + x[2]*lenclass
  if(type=="logit") {
    mu=(x[3]*exp(eta))/(1-x[3]+exp(eta))
    dmudeta= x[3]*(1-x[3])*exp(eta) / ( (1-x[3]+exp(eta))^2 )
    dmudp= exp(eta)*(1+exp(eta)) / ( (1-x[3]+exp(eta))^2 ) }
  if(type=="cl") {
    dblexpo=exp(-exp(eta))
    mu=1-(1-x[3])/(1-x[3]*dblexpo)
    dmudeta=(1-x[3])*x[3]*exp(eta)*dblexpo/ ( (1-x[3]*dblexpo)^2 )
    dmudp=(1-dblexpo)/( (1-x[3]*dblexpo)^2 ) }
  if(type=="logit" | type=="cl") {
    info=matrix(0,nrow=3,ncol=3)
    info[1,1]=sum( ntotal*dmudeta^2 / (mu*(1-mu)) )
    info[1,2]=sum( ntotal*lenclass*dmudeta^2 / (mu*(1-mu)) )
    info[1,3]=sum( ntotal*dmudp*dmudeta / (mu*(1-mu)) )
    info[2,1]=info[1,2]
    info[2,2]=sum( ntotal*(lenclass*dmudeta)^2 / (mu*(1-mu)) )
    info[2,3]=sum( ntotal*lenclass*dmudeta*dmudp / (mu*(1-mu)) )
    info[3,1]=info[1,3]
    info[3,2]=info[2,3]
    info[3,3]=sum( ntotal*dmudp^2 / (mu*(1-mu)) )
    covar=solve(info,tol=1e-12) 
    list(covar=covar)  }
  else return(NA)  }


#Covariance matrix for logistic fit to covered codend data
cccov=function(x,catch) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  ntotal=nfine+nwide
  eta=x[1] + x[2]*lenclass
  mu=exp(eta)/(1+exp(eta))
  dmudeta=mu/(1+exp(eta))
  info=matrix(0,nrow=2,ncol=2)
  info[1,1]=sum( ntotal*dmudeta^2 / (mu*(1-mu)) )
  info[1,2]=sum( ntotal*lenclass*dmudeta^2 / (mu*(1-mu)) )
  info[2,1]=info[1,2]
  info[2,2]=sum( ntotal*(lenclass*dmudeta)^2 / (mu*(1-mu)) )
  covar=solve(info)
  list(covar=covar) }

#Covariance matrix for Richard's fit to covered codend data
cccovrich=function(x,catch) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  ntotal=nfine+nwide
  eta=x[1] + x[2]*lenclass
  nu=exp(eta)/(1+exp(eta)); mu=nu^(1/x[3])
  eta <- ifelse(mu == 1, Inf, eta) # To avoid numerical problems below
  dmudeta=(1/x[3])*mu/(1+exp(eta))
  dmudgamma=-x[3]^(-2)*mu*log(nu)
  info=matrix(0,nrow=3,ncol=3)
  info[1,1]=sum( ntotal*dmudeta^2 / (mu*(1-mu)) ,na.rm=T)
  info[1,2]=sum( ntotal*lenclass*dmudeta^2 / (mu*(1-mu)) ,na.rm=T)
  info[1,3]=sum( ntotal*dmudgamma*dmudeta / (mu*(1-mu)) ,na.rm=T)
  info[2,1]=info[1,2]
  info[2,2]=sum( ntotal*(lenclass*dmudeta)^2 / (mu*(1-mu)) ,na.rm=T)
  info[2,3]=sum( ntotal*lenclass*dmudeta*dmudgamma / (mu*(1-mu)) ,na.rm=T)
  info[3,1]=info[1,3]
  info[3,2]=info[2,3]
  info[3,3]=sum( ntotal*dmudgamma^2 / (mu*(1-mu)) ,na.rm=T)
  covar=solve(info)
  list(covar=covar) }

#Covariance matrix for Richard's fit to trouser trawl data
covrich=function(x,catch,npars=3,p=0.5) {
  lenclass=catch[,1]; nwide=catch[,2]; nfine=catch[,3]
  ntotal=nfine+nwide
  eta=x[1] + x[2]*lenclass 
  if(npars==2) cat("\n","covrich: Fixed split p= ",p) else p=x[4]
  nu=exp(eta)/(1+exp(eta)); r=nu^(1/x[3]); mu=p*r/(1-p+p*r)
  dmudr=p*(1-p)/((1-p+p*r)^2)
  drdeta=(1/x[3])*r/(1+exp(eta)); dmudeta=dmudr*drdeta
  drdgamma=-x[3]^(-2)*r*log(nu);  dmudgamma=dmudr*drdgamma
  dmudp=r/((1-p+p*r)^2)
  info=matrix(0,nrow=4,ncol=4)
  info[1,1]=sum( ntotal*dmudeta^2 / (mu*(1-mu)) )
  info[1,2]=sum( ntotal*lenclass*dmudeta^2 / (mu*(1-mu)) )
  info[1,3]=sum( ntotal*dmudgamma*dmudeta / (mu*(1-mu)) )
  info[1,4]=sum( ntotal*dmudeta*dmudp / (mu*(1-mu)) )
  info[2,1]=info[1,2]
  info[2,2]=sum( ntotal*(lenclass*dmudeta)^2 / (mu*(1-mu)) )
  info[2,3]=sum( ntotal*lenclass*dmudeta*dmudgamma / (mu*(1-mu)) )
  info[2,4]=sum( ntotal*lenclass*dmudeta*dmudp / (mu*(1-mu)) )
  info[3,1]=info[1,3]
  info[3,2]=info[2,3]
  info[3,3]=sum( ntotal*dmudgamma^2 / (mu*(1-mu)) )
  info[3,4]=sum( ntotal*dmudgamma*dmudp / (mu*(1-mu)) )
  info[4,1]=info[1,4]
  info[4,2]=info[2,4]
  info[4,3]=info[3,4]
  info[4,4]=sum( ntotal*dmudp*dmudp / (mu*(1-mu)) )
  if(npars==3) covar=solve(info) else covar=solve(info[1:3,1:3])
  list(covar=covar) }


