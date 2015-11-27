###################################################################################
################################ ttfit example ####################################
###################################################################################
#Code to plot the fitted selection curve is provided at the bottom of this example#
#
#NB: The Richards curve fit is not successful for these data due to 
#    high correlations in parameters, i.e., it is too complicated for these data.

#source("trawlfunctions.R")  #Load trawlfunctions if not already installed
#attach("C:/Documents and Settings/rmil013/My Documents/trawl/R/.RData")

haddock=matrix(scan("haddock.dat",skip=1),ncol=3,byrow=T)

#For this example, columns 2 and 3 need to be swapped because ttfit() and ccfit() 
#assume that the experimental catch numbers are in column 2 and the control
#catch numbers in column 3
haddockdat=haddock[,c(1,3,2)]

#Set up two plots on the figure (for observed propns and deviance resids)
par(mfrow=c(2,1))

#Fit the model and save as haddock.fit
haddock.fit=ttfit(haddockdat)

#See the parameters and diagnostics
haddock.fit
# Output is:

# Converged is 1 or 2
# $converged
# [1] 1

# Parameters a, b and p
#$x
#[1] -27.6544451   0.9164563   0.5728422

#Log-likelihoods of fitted model, null model, and full model
#$l
#          [,1] [,2]
#[1,] -1016.828   21
#[2,] -1056.988   23
#[3,] -1009.410    0

# Covariance matrix of parameters
# $covar
            [,1]         [,2]          [,3]
# [1,] 33.40173818 -1.145220650  0.0473210147
# [2,] -1.14522065  0.039341497 -0.0016948297
# [3,]  0.04732101 -0.001694830  0.0003009187

# Estimates of L25, L50, l75 (default) and their s.e's
# $lens
         [,1]      [,2]
# [1,] 28.97665 0.2789660
# [2,] 30.17541 0.3608154
# [3,] 31.37417 0.5631851

# Estimate of selection range and its s.e.
# $sr
# [1] 2.3975224 0.5188912

# Estimated retention probs for each lengthclass
#$r
# [1] 0.003472186 0.008636906 0.021319562 0.051655373 0.119869409 0.254034160
# [7] 0.459897645 0.680419785 0.841863323 0.930125212 0.970831701 0.988126797
#[13] 0.995217431 0.998081784 0.999231957 0.999692692 0.999877074 0.999950834
#[19] 0.999980336 0.999992136 0.999996855 0.999998742 0.999999497 0.999999799

# Deviance residuals
# $devres
#  [1] -0.09639071 -0.15176298 -0.41126355  0.02947730  0.07448783  0.47433458
#  [7] -1.09201911  0.55966541  0.88955053 -0.87087881  0.20954439 -0.27074318
# [13] -0.91090234  0.22714660 -0.42580615  0.14651504  1.96667109  1.17173377
# [19] -0.65497361  1.07441031 -0.50736883 -0.29293141  1.36598722  1.05559950

# Check of sufficient data for reliable LR statistic
# Based on expected count and total minus expected count.
# (Don't worry about a few FALSES's)
# $suff.dat
#  [1] FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
# [13]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE FALSE FALSE FALSE


#################################
#Plot the fitted selection curve#
#################################
mtext("Haddock data",cex=2,outer=T)
plot(haddockdat[,1],haddock.fit$r,type="l",xlab="Length (cm)",las=1,
      ylab="Retention probability",main="Haddock selection by 87 mm mesh") 
abline(h=c(0.25,0.5,0.75),lty=2)

#To change the length values used for the plot, use the plotlens argument
mtext("Haddock data",cex=2,outer=T)
plens=seq(20,50,0.1)
r=ttfit(haddockdat,plotlens=plens,plots=F)$r
plot(plens,r,type="l",xlab="Length (cm)",las=1,
      ylab="Retention probability",main="Haddock selection by 87 mm mesh") 
abline(h=c(0.25,0.5,0.75),lty=2)

 
