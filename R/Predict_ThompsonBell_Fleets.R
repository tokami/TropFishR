#' @title Thompson and Bell prediction model
#
#' @description  This is a function to calculate the total mortality (Z) from length composition data via the length converted catch curve or from age at length data with catch curve.
#'
#' @param param A list
#' @param stock_size_1 Stock size to start with
#' @param plus.group Should a plus group be created, by default (NA) not. BUt authors advise to do so
#' @param FM_change Vector containing new FM values
#' @param fleet_mat Matrix providing catch or FM values separated by fleets
#' @param fleet_unit Either 'Catch' or 'FM' indicating in which unit the fleet_mat is provided
#' @param fleet_FM_change Matrix containing new FM values for each fishery
#' @param fleet_plot_name Name of fishery, which should be used for plotting Yield per recurit against changes of FM
#' @param unit.time Time unit "year" or "month", Default "year"
#'
#' @examples
#'
#' # load data
#' data(data_Predict_ThompsonBell)
#'
#' # construct fleet matrix
#' fleet_mat=matrix(ncol = 2,nrow = 12)
#' fleet_mat[,1] <- c(0.72,0.96,0.84,0.48,0.6,0.48,1.08,0.48,0.084,0.12,0.24,0.24)
#' fleet_mat[,2] <- c(0.48,0.36,0.48,0.96,1.32,0.72,0.48,0.72,1.116,1.68,2.52,2.28)
#' colnames(fleet_mat) <- c("artisanal", "industrial")
#'
#' # construct matrix with changes of fishing mortalities by fleet
#' fleet_FM_change = matrix(ncol = 2,nrow = 8)
#' fleet_FM_change[,1] <- rep(1,8)
#' fleet_FM_change[,2] <- c(0,0.4,0.8,1.0,1.2,1.5,2,3)
#' colnames(fleet_FM_change) <- c("artisanal", "industrial")
#'
#' # run the model
#' output <- Predict_ThompsonBell_Fleets(data_Predict_ThompsonBell,seq(0,3,0.2),
#'    fleet_mat,"FM",fleet_FM_change,fleet_plot_name = "industrial",
#'    'month',plus.group = 12)
#'
#' # investigate results
#' output
#'
#' @details better to treat last group always as a plus group..... For variable parameter system vectors are reuqired for constant parameter systems matrices or data.frames have to be inserted. or vectors The length converted linearised catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @references example 1 : Kuwait (Garcia and van zalinge 1982)
#'   Millar, R. B., & Holst, R. (1997). Estimation of gillnet and hook selectivity using log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
#' @export

Predict_ThompsonBell_Fleets <- function(param, FM_change,fleet_mat, fleet_unit,
                                        fleet_FM_change,fleet_plot_name = NA,
                                        unit.time = 'year', stock_size_1 = NA,
                                        plus.group = NA){

  res <- param
  meanWeight <- res$meanWeight
  meanValue <- res$meanValue

  #mortalities
  FM <- res$FM
  if(!is.null(res$M)){
    nM <- res$M
    Z <- FM + nM
  }else{
    Z <- res$Z
    nM <- mean(Z - FM,na.rm=T)
  }


  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                        Age data                          #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if('age' %in% names(res)) classes <- as.character(res$age)
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                       Length data                        #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if('midLengths' %in% names(res)) classes <- as.character(res$midLengths)

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
  classes.num <- as.numeric(classes.num[,1])

  if(is.na(plus.group)){
    if(length(classes.num[1:length(classes.num)]) > length(FM)) stop("Please provide fishing mortality values for all age/length groups!")
  }
  if(!is.na(plus.group)){
    if(length(classes.num[1:plus.group]) > length(FM)) stop("Please provide fishing mortality values for all age/length groups!")
    if(length(classes.num[1:plus.group]) > length(FM)){
      FM <- FM[1:plus.group]
      Z <- Z[1:plus.group]
    }
  }

  #prediction by fleet
  res2 <- Predict_ThompsonBell1(res, stock_size_1 = stock_size_1,
                                plus.group = plus.group, unit.time = unit.time)

  #correct fleet_mat if provided in Catch not fishing mortalities
  if(fleet_unit == 'Catch'){
    total.FM <- sum(res2$FM,na.rm = TRUE)
    fleet_mat.FM <- (res2$tot.C * fleet_mat ) / total.FM
  }                                                      #fleet_mat.Catch needed? ???
  if(fleet_unit == 'FM'){
    fleet_mat.FM <- fleet_mat
  }

  pred_fleet_mat_list <- list()
  for(x8 in 1:dim(fleet_mat)[2]){
    pred.fleet_mat <- as.matrix(fleet_mat[,x8]) %*% fleet_FM_change[,x8]
    colnames(pred.fleet_mat) <- fleet_FM_change[,x8]

    pred_fleet_mat_list[[x8]] <- pred.fleet_mat
  }


  fleet_FM.list <- list()
  fleet_Yield.list <- list()
  fleet_Biomass.list <- list()
  fleet_Value.list <- list()

  for(x9 in 1:dim(pred_fleet_mat_list[[1]])[2]){
    #FISHING MORTALITY
    FM_fleets.list <- lapply(pred_fleet_mat_list, function(x) x[,x9])
    FM_fleets.df <- do.call(cbind,FM_fleets.list)
    FM_fleets.df <- as.data.frame(FM_fleets.df)
    FM_fleets.df$FMtot <- rowSums(FM_fleets.df,na.rm=TRUE)

    param.loop <- res2
    param.loop$FM <- FM_fleets.df$FMtot
    param.loop$Z <- FM_fleets.df$FMtot + nM

    if(is.na(plus.group)){
      if(length(classes.num[1:length(classes.num)]) > dim(fleet_mat)[1]) stop("Please provide fishing mortality values for all age/length groups!")
    }
    if(!is.na(plus.group)){
      if(length(classes.num[1:plus.group]) > dim(fleet_mat)[1]) stop("Please provide fishing mortality values for all age/length groups!")
      if(length(classes.num[1:plus.group]) > dim(fleet_mat)[1]){
        param.loop$FM <- param.loop$FM[1:plus.group]
        param.loop$Z <- param.loop$Z[1:plus.group]
      }
    }

    res3 <- Predict_ThompsonBell1(param.loop,unit.time = unit.time,
                                  stock_size_1 = stock_size_1, plus.group = plus.group)

    res3$dt <- as.numeric(as.character(res3$dt))
    res3$V <- as.numeric(as.character(res3$V))
    FM_fleets.df$Ctot <- as.numeric(as.character(res3$C))   #or tot.C
    FM_fleets.df$w <- as.numeric(as.character(res3$meanWeight))
    Ftots <-  colSums(FM_fleets.df[,1:3],na.rm=TRUE)

    #YIELD
    Y_fleets.list <- list()
    for(x10 in 1:length(FM_fleets.list)){
      Y_fleets.list[[x10]] <- ((FM_fleets.df$Ctot * FM_fleets.list[[x10]]) /
                                 FM_fleets.df$FMtot ) * FM_fleets.df$w
    }
    Y_fleets.df <- do.call(cbind,Y_fleets.list)
    Y_fleets.df <- as.data.frame(Y_fleets.df)
    Y_fleets.df$Ytot <- rowSums(Y_fleets.df,na.rm=TRUE)

    Ytots <- colSums(Y_fleets.df,na.rm=TRUE)

    #BIOMASS
    B_fleets.list <- list()
    for(x11 in 1:length(FM_fleets.list)){
      B_fleets.list[[x11]] <- Y_fleets.list[[x11]] /
        (FM_fleets.list[[x11]] * res3$dt)
    }
    B_fleets.df <- do.call(cbind,B_fleets.list)
    B_fleets.df <- as.data.frame(B_fleets.df)
    B_fleets.df$Btot <- rowSums(B_fleets.df,na.rm=TRUE)

    Btots <- colSums((B_fleets.df * res3$dt),na.rm=TRUE) /
      sum(res3$dt,na.rm=TRUE)

    #VALUE
    V_fleets.list <- list()
    for(x12 in 1:length(FM_fleets.list)){
      V_fleets.list[[x12]] <- Y_fleets.list[[x12]] * res3$V
    }
    V_fleets.df <- do.call(cbind,V_fleets.list)
    V_fleets.df <- as.data.frame(V_fleets.df)
    V_fleets.df$Vtot <- rowSums(V_fleets.df,na.rm=TRUE)

    Vtots <- colSums(V_fleets.df,na.rm=TRUE)

    fleet_FM.list[[x9]] <- Ftots
    fleet_Yield.list[[x9]] <- Ytots
    fleet_Biomass.list[[x9]] <- Btots
    fleet_Value.list[[x9]] <- Vtots
  }

  fleet_FM.df <- do.call(rbind,fleet_FM.list)
  fleet_Yield.df <- do.call(rbind,fleet_Yield.list)
  fleet_Biomass.df <- do.call(rbind,fleet_Biomass.list)
  fleet_Value.df <- do.call(rbind,fleet_Value.list)

  fleet_FM.df <- as.data.frame(cbind(fleet_FM.df,fleet_FM_change))
  fleet_Yield.df <- as.data.frame(cbind(fleet_Yield.df,fleet_FM_change))
  fleet_Biomass.df <- as.data.frame(cbind(fleet_Biomass.df,fleet_FM_change))
  fleet_Value.df <- as.data.frame(cbind(fleet_Value.df,fleet_FM_change))


  #dimensions of plotting area
  max_val <- round(max(fleet_Value.df$Vtot,na.rm=TRUE),digits=0)
  dim_val <- 10 ^ (nchar(max_val)-1)
  max_yiel <- round(max(fleet_Yield.df$Ytot,na.rm=TRUE),digits=0)
  dim_yiel <- 10 ^ (nchar(max_yiel)-1)
  max_bio <- round(max(fleet_Biomass.df$Btot,na.rm=TRUE),digits=0)
  dim_bio <- 10 ^ (nchar(max_bio)-1)

  #some info for plotting
  num_fleets <- dim(fleet_mat)[2]
  if(!is.na(fleet_plot_name)){
    plot_fleet <- which(colnames(fleet_mat) == fleet_plot_name)
  }else if(is.na(fleet_plot_name)){
    rgi.list <- list()
    for(x13 in 1:dim(fleet_FM_change)[2]){
      stri <- 1
      sopi <- dim(fleet_FM_change)[1] * x13

      rgi.list[[x13]] <- range(fleet_FM_change[stri:sopi],na.rm=TRUE)

      stri <- 1 + dim(fleet_FM_change)[1]
    }

    rgi.upd <- lapply(rgi.list,function(x) x[2]-x[1])
    rgi.upd <- do.call(rbind,rgi.upd)
    plot_fleet <- which(rgi.upd == max(rgi.upd,na.rm=TRUE))
  }



  #4 panel plot
  op <- par(new=FALSE,mar = c(4, 4, 4, 4) + 0.1, mfrow=c(2,2))
  plot(fleet_FM_change[,plot_fleet],fleet_Value.df$Vtot, type ='o',main='Value',
       col ='darkorange', ylim = c(0,ceiling(max_val/dim_val)*dim_val),
       lwd=1.6, xlab = "F-factor X", ylab = "Value",lty=1)
  abline(v=1,col='grey',lty=2)
  for(j in 1:dim(fleet_FM_change)[2]){
    points(fleet_FM_change[,plot_fleet],fleet_Value.df[,j],
           col ='darkorange',lwd=1.6,lty = (j+1),type='o')
  }
  plot(fleet_FM_change[,plot_fleet], fleet_Yield.df$Ytot,
       ylim = c(0,ceiling(max_yiel/dim_yiel)*dim_yiel),
       lty = 1, col='dodgerblue',lwd=1.6,main='Yield',
       type='o',ylab='Yield',xlab='F-factor X')
  abline(v=1,col='grey',lty=2)
  for(j in 1:dim(fleet_FM_change)[2]){
    points(fleet_FM_change[,plot_fleet],fleet_Yield.df[,j],
           lty = (j+1), col='dodgerblue',lwd=1.6,type='o')
  }
  plot(fleet_FM_change[,plot_fleet],fleet_Biomass.df$Btot,
       ylim = c(0,ceiling(max_bio/dim_bio)*dim_bio),
       lty = 1, col = 'darkgreen',lwd=1.6,main='Biomass',
       type='o',ylab='Biomass',xlab='F-factor X')
  abline(v=1,col='grey',lty=2)
  for(j in 1:dim(fleet_FM_change)[2]){
    points(fleet_FM_change[,plot_fleet],fleet_Biomass.df[,j],
           lty = (j+1), col = 'darkgreen',lwd=1.6,type='o')
  }
  plot(1,axes=F,xlab='',ylab='',type='n')
  legend("center", c("total",colnames(fleet_mat)),bty = "n",
         lty = 1:(dim(fleet_FM_change)[2]+1), seg.len = 0.2,
         col = 'black',lwd=2,xpd=TRUE,y.intersp = 0.15,
         text.width=0.3,x.intersp=0.3)
  par(op)

  res4 <- list(fleet_FM = fleet_FM.df,
               fleet_Yield = fleet_Yield.df,
               fleet_Biomass = fleet_Biomass.df,
               fleet_Value = fleet_Value.df)
  ret <- c(res,res2,res4)
  return(ret)
} ## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?

