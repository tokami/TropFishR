#' @title Length converted catch curve
#'
#' @description
#'
#' @param
#'
#' @details
#'
#' @export

CatchCurve.R <- function(){

  ####   CALCULATIONS
  # t of lower length classes
  data$lowerLength <- data$midLength - (para$interval[i] / 2)
  data$catch_cor <- data$catch * para$catch_cor_fac[i]
  data$t_L1 <- (para$t0[i] - (1/para$k[i])) * log(1 - (data$lowerLength / para$L8[i]))
  # delta t
  data$dt <- NA
  for(k in 1:(length(data$dt)-1)){
    data$dt[k] <- data$t_L1[k+1]-data$t_L1[k]
  }
  # t of midlengths
  data$t_midL <- (para$t0[i] - (1/para$k[i])) * log(1 - (data$midLength / para$L8[i]))
  # ln( Catch / delta t)
  data$lnC_dt <- log(data$catch_cor - data$dt)


  #identify plot
  plot(x=data$t_midL,y=data$lnC_dt,
       xlab="Realtive age [yrs - t0]",ylab="ln(N/dt)",
       ylim=c(0,10))
  cutter <- identify(x=data$t_midL,y=data$lnC_dt,labels=rownames(data))

  #calculations + model
  data_cut <- data[cutter[1]:cutter[2],]
  lm1 <- lm(lnC_dt ~ t_midL,data=data_cut)
  sum_lm1 <- summary(lm1)
  r_lm1 <- sum_lm1$r.squared
  intercept_lm1 <- sum_lm1$coefficients[1]
  slope_lm1 <- sum_lm1$coefficients[2]
  se_slope_lm1 <- sum_lm1$coefficients[4]

  sum_lm1$r.squared  ########## NEW

  Z_lm1 <- abs(slope_lm1)
  SE_Z_lm1 <- abs(se_slope_lm1)

  #save plot
  dev.off()
  setwd(direc_data)
  pdf(paste("catchi_2_",speci,".pdf",sep=''),height=6,width=6)
  #final plot
  plot(x=data$t_midL,y=data$lnC_dt, main = speci,
       xlab="Realtive age [yrs - t0]",ylab="ln(N / dt)",
       ylim=c(0,12),xlim=c(0,10),cex=1.5,bty='n')
  par(new=T)
  points(x=data_cut$t_midL,y=data_cut$lnC_dt,
         pch=19,col='blue',cex=1.5)
  ablineclip(a = intercept_lm1,b=slope_lm1,lwd=1.7,col='blue',x1 = data$t_midL[cutter[1]],x2 = data$t_midL[cutter[2]])
  mtext(side = 3,text = paste("Z =",round(Z_lm1,2),"±",round(SE_Z_lm1,2)),col='blue')
  dev.off()



  #save Z
  results_list[[i]] <- data.frame(species = speci,
                                  Z = Z_lm1,
                                  SE = SE_Z_lm1)
}












Zti <- do.call(rbind,results_list)

#write excel table
setwd(direc_data)
WriteXLS('Zti','Zti.xls',FreezeRow = 1,BoldHeaderRow = T,row.names = F)



##############
#for two

#identify plot
plot(x=data$t_midL,y=data$lnC_dt,
     xlab="Realtive age [yrs - t0]",ylab="ln(N/dt)",
     ylim=c(0,10))
cutter2 <- identify(x=data$t_midL,y=data$lnC_dt,labels=rownames(data))

#calculations + model
data_cut2 <- data[cutter2[1]:cutter2[2],]
lm2 <- lm(lnC_dt ~ t_midL,data=data_cut2)
sum_lm2 <- summary(lm2)
r_lm2 <- sum_lm2$r.squared
intercept_lm2 <- sum_lm2$coefficients[1]
slope_lm2 <- sum_lm2$coefficients[2]
se_slope_lm2 <- sum_lm2$coefficients[4]

Z_lm2 <- abs(slope_lm2)
SE_Z_lm2 <- abs(se_slope_lm2)


#save plot
dev.off()
setwd(direc_data)
pdf(paste("catchi_4_",speci,".pdf",sep=''),height=6,width=6)
#final plot
plot(x=data$t_midL,y=data$lnC_dt, main = speci,
     xlab="Realtive age [yrs - t0]",ylab="ln(N / dt)",
     ylim=c(0,12),xlim=c(0,10),cex=1.5,bty='n')
par(new=T)
points(x=data_cut$t_midL,y=data_cut$lnC_dt,
       pch=19,col='darkgreen',cex=1.5)
ablineclip(a = intercept_lm1,b=slope_lm1,lwd=1.7,col='darkgreen',x1 = min(data_cut$t_midL),x2 = max(data_cut$t_midL))
mtext(side = 3,text = paste("Z =",round(Z_lm1,2),"±",round(SE_Z_lm1,2)),col='darkgreen')
points(x=data_cut2$t_midL,y=data_cut2$lnC_dt,
       pch=19,col='goldenrod3',cex=1.5)
ablineclip(a = intercept_lm2,b=slope_lm2,lwd=1.7,col='goldenrod3',x1 = min(data_cut2$t_midL),x2 = max(data_cut2$t_midL))
mtext(side = 3,line=-1.2,text = paste("Z =",round(Z_lm2,2),"±",round(SE_Z_lm2,2)),col='goldenrod3')
dev.off()

