source("NextGeneration.R")
source("SelnCurveDefinitions.R") #These can be extended by the user

###################################
####Analyse trouser trawl data#####
###################################
#NOTE: the horizontal line in the PlotCurves figure is the standardized 
#selectivity of the control gear
par(mfrow=c(2,1))
haddock.dat=read.table("haddock.dat",head=T) #control is 2nd col
fit=NetFit(haddock.dat,c(0,130),c(-10,0.3,0),rtype="tt.logistic")
Estimates(fit); Summary(fit)
PlotCurves(fit,plotlens=seq(25,35,0.1))
legend(32,0.5,c("Control","Experimental"),lty=1:2,col=1:2,cex=0.8,bg="white")

###################################
########Analyse gillnet data#######
###################################
holt.dat=read.table("holt.dat",head=F); 
Meshsize=c(13.5,14,14.8,15.4,15.9,16.6,17.8,19) 
#Equal fishing power
pwr=rep(1,8)
#Or use this for fishing power proportional to meshsize
#pwr=Meshsize

par(mfrow=c(5,2),mar=c(4.1,4.1,1,1))
fit=NetFit(holt.dat,Meshsize,c(60,4),rtype="norm.loc",rel.power=pwr)
Estimates(fit); Summary(fit); PlotCurves(fit,plotlens=seq(40,90,0.1))

fit=NetFit(holt.dat,Meshsize,c(60,4),rtype="norm.sca",rel.power=pwr)
Estimates(fit); Summary(fit); PlotCurves(fit,plotlens=seq(40,90,0.1))

#fit=NetFit(holt.dat,Meshsize,c(4,0.1),rtype="gamma") #Gamma not yet implemented

fit=NetFit(holt.dat,Meshsize,c(4,0.1),rtype="lognorm",rel.power=pwr)
Estimates(fit); Summary(fit); PlotCurves(fit,plotlens=seq(40,90,0.1))

fit=NetFit(holt.dat,Meshsize,c(55,4,65,4,2),rtype="binorm.sca",rel.power=pwr)
Estimates(fit); Summary(fit); PlotCurves(fit,plotlens=seq(40,90,0.1))

fit=NetFit(holt.dat,Meshsize,c(4,0.2,4.2,0.1,2),rtype="bilognorm",rel.power=pwr)
Estimates(fit); Summary(fit); PlotCurves(fit,plotlens=seq(40,90,0.1))

#NOTE: The bilognorm (mixture of lognormals) is preferred here.

########################################
####Analyse stacked trammel net data####
########################################
#Here, the data are from two expts using different mesh sizes
#This analysis assumes common retention curve in both expts.
#Note that summary function does not produce residual plot
#since lengths are not unique
T.dat=read.table("trammel.dat",head=F); 
#T.dat=as.matrix(T.dat) #Needed only if a complete column of NA's
Meshsize=c(6.1,7.6,7.9,9.1,10.6,13) 
#Equal fishing power (within each of the two expts)
pwr=rep(1,6)
#Or use this for fishing power proportional to meshsize
#pwr=Meshsize

par(mfrow=c(3,2),mar=c(4.1,4.1,1,1))
fit=NetFit(T.dat,Meshsize,c(25,4),rtype="norm.loc",rel.power=pwr)
Estimates(fit); Summary(fit); PlotCurves(fit,plotlens=seq(10,40,0.1))

fit=NetFit(T.dat,Meshsize,c(25,4),rtype="norm.sca",rel.power=pwr)
Estimates(fit); Summary(fit); PlotCurves(fit,plotlens=seq(10,40,0.1))

#fit=NetFit(T.dat,Meshsize,c(3,0.1),rtype="gamma") #Gamma not yet implemented

fit=NetFit(T.dat,Meshsize,c(3,0.1),rtype="lognorm",rel.power=pwr)
Estimates(fit); Summary(fit); PlotCurves(fit,plotlens=seq(10,40,0.1))

fit=NetFit(T.dat,Meshsize,c(18,4,30,4,2),rtype="binorm.sca",rel.power=pwr)
Estimates(fit); Summary(fit); PlotCurves(fit,plotlens=seq(10,40,0.1))

fit=NetFit(T.dat,Meshsize,c(3,0.2,3.2,0.1,2),rtype="bilognorm",rel.power=pwr)
Estimates(fit); Summary(fit); PlotCurves(fit,plotlens=seq(10,40,0.1))

