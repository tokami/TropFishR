#' @title Thompson and Bell prediction model
#
#' @description  This is a function to calculate the total mortality (Z) from length composition data via the length converted catch curve or from age at length data with catch curve.
#'
#' @param param A list
#' @param mean_weight Average weight of fish
#' @param datatype Type of data which is used for analysis, either 'length' or 'age', for length frequency or age composition data, respectively
#' @param FM Fishing mortality
#' @param Z Total mortality
#' @param value Information about the value/price of fish per kilo or gramm
#' @param Tr Age of recruitment
#' @param stock_size_1 Stock size to start with
#' @param plus.group Should a plus group be created, by default (NA) not. BUt authors advise to do so
#' @param FM_change Vector containing new FM values
#' @param fleet_mat Matrix providing catch or FM values separated by fleets
#' @param fleet_unit Either 'Catch' or 'FM' indicating in which unit the fleet_mat is provided
#' @param fleet_FM_change Matrix containing new FM values for each fishery
#' @param fleet_plot_name Name of fishery, which should be used for plotting Yield per recurit against changes of FM
#'
#' @param FM reference F-at-age-array
#'
#' @examples
#' #data("ex.Predict_ThompsonBell")
#' #with(ex.Predict_ThompsonBell,Predict_ThompsonBell(classes = age,
#' #  mean_weight = meanWeight,value = valueGramm,FM = FM,Z = Z,
#' #  datatype = 'age',Tr = 1,stock_size_1 = NA,plus.group = 12))
#'
#' @details better to treat last group always as a plus group..... For variable parameter system vectors are reuqired for constant parameter systems matrices or data.frames have to be inserted. or vectors The length converted linearised catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @references example 1 : Kuwait (Garcia and van zalinge 1982)
#'   Millar, R. B., & Holst, R. (1997). Estimation of gillnet and hook selectivity using log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
#' @export


data(data_Predict_ThompsonBell)
param = data_Predict_ThompsonBell
unit.time = "year"


#param can be age or midLenghts



Predict_ThompsonBell <- function(param, select.param = NA, unit.time = "year",
                                 stock_size_1 = NA, plus.group = NA,


   Tr,stock_size_1 = NA,
   FM_change = NA, fleet_mat = NA, fleet_unit = NA,
   fleet_FM_change = NA, fleet_plot_name = NA,
   s_list = NA, Linf = NA, K = NA, t0 = 0, L50 = NA,
   L75 = NA){



  res <- param
  meanWeight <- res$meanWeight
  meanValue <- res$meanValue

  #mortalities
  FM <- res$FM
  Z <- res$Z
  #natural Mortality
  nM <- mean(Z - FM,na.rm=T)

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                        Age data                          #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if('age' %in% names(res)){

    classes <- as.character(res$age)
    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    #function with all calculations
    calc.age <- function(classes.num, FM, Z, meanWeight,
                         meanValue, unit.time, stock_size_1){

      # delta t
      dt <- rep(NA,length(classes.num))
      for(x0 in 1:(length(classes.num)-1)){
        dt[x0] <- classes.num[x0+1] - classes.num[x0]
      }
      if(unit.time == 'month'){
        dt <- dt * 1/12
      }

      #population per age group
      pop <- rep(NA,length(classes.num))
      if(!is.na(stock_size_1)){
        pop[1] <- stock_size_1
      }else pop[1] <- 1000  #if not provided results are just negative

      for(x1 in 2:length(classes.num)){
        pop[x1] <- pop[x1-1] * exp(-Z[x1] * dt[x1])
        if(x1 == length(pop)){                                ####CORRECT?????
          pop[x1] <- pop[x1-1] * exp(-Z[x1-1] * dt[x1-1])
        }
      }

      #number of deaths per time step month or year
      dead <- NA
      for(x2 in 1:length(classes.num)){
        dead[x2] <- pop[x2] - pop[x2+1]
      }

      #catch in numbers
      C <- dead * FM / Z

      #yield in kg or g
      Y <- C * meanWeight

      #mean biomass in kg or g
      B <- Y / (FM * dt)

      #value expressed in money units
      value.Y <- Y * meanValue

      #total catch, yield, value and average biomass
      tot.C <- sum(C, na.rm=TRUE)
      tot.Y <- sum(Y, na.rm=TRUE)
      tot.value <- sum(value.Y, na.rm=TRUE)
      meanB <- sum((B * dt), na.rm=TRUE) / sum(dt, na.rm=TRUE)   ### more complicated biomass concept if dt is not constant, see Chapter 5

      calc.res <- list(dt = dt, pop = pop, dead = dead, C = C,
                       Y = Y, B = B,value.Y = value.Y,
                       tot.C = tot.C, tot.Y = tot.Y, meanB = meanB,
                       tot.value = tot.value
        )

      return(calc.res)
    }


    res2 <- calc.age(classes.num,  FM, Z, meanWeight, meanValue,
                              unit.time, stock_size_1)


    plus.group=11

    ########### GO ON FROM HERE! DOES IT MAKE SENSE TO DO IT WITHOUT DATAFRAMES?


    #with plus group
    if(!is.na(plus.group)){

      df <- do.call(cbind,res2[1:7])
      df <- df[1:plus.group,]
      classes <- classes[1:plus.group]

      #new class
      classes[length(classes)] <-
        paste(classes[length(classes)],"+",sep='')

      #num deaths
      df[plus.group, "dead"] <- df[plus.group, "pop"]

      #catch
      new.C <- (FM[plus.group] / Z[plus.group]) * df[plus.group, "pop"]
      catch.plus.dif <- new.C - df[plus.group, "C"]
      df[plus.group, "C"] <- new.C

      #yield
      df[plus.group, "Y"] <- meanWeight[plus.group] * catch.plus.dif

      #value
      df[plus.group, "value.Y"] <-
        df[plus.group, "Y"] * meanValue[plus.group]

      #biomass       ####not sure....omitted in manual
      df[plus.group, "B"] <-
        df[plus.group, "Y"] / (FM[plus.group] * df[plus.group, "dt"])
    }
    # does not cut all vectors from beginning!!! careful if you return "res" in the end!

    #prediction based on f_change
    if(!is.na(plus.group)){
      df.TBpred <- df.TBnew
    }else df.TBpred <- df.TB

    pred_mat <- as.matrix(df.TBpred$FM) %*% FM_change
    colnames(pred_mat) <- FM_change

    df.TBpred.loop <- df.TBpred
    pred_res_list <- list()
    for(x7 in 1:length(FM_change)){
      df.TBpred.loop$FM <- pred_mat[,x7]
      df.TBpred.loop$Z <- pred_mat[,x7] + nM
      res <- calc_TB.age(df.TB = df.TBpred.loop, unit.time = unit.time, stock_size_1 = stock_size_1)
      pred_res_list[[x7]] <- res$totals
    }

    pred_res_df <- do.call(rbind, pred_res_list)
    pred_res_df$Xfact <- FM_change
    rownames(pred_res_df) <- FM_change


    #save x axis positions
    max_val <- round(max(pred_res_df$total.value,na.rm=TRUE),digits=0)
    dim_val <- 10 ^ (nchar(max_val)-1)

    par(oma = c(1, 1, 1.5, 1),new=FALSE,mar = c(5, 4, 4, 4) + 0.3)
    plot(pred_res_df$Xfact,pred_res_df$total.value, type ='l',
         col ='darkorange', ylim = c(0,ceiling(max_val/dim_val)*dim_val),
         lwd=1.6, xlab = "F-factor X", ylab = "Value")
    lines(pred_res_df$Xfact,pred_res_df$mean.biomass,
          col = 'darkgreen',lwd=1.6)    # draw lines with small intervals: seq(0,max(),0.05) but y as to be dependent of x (formula of calculaiton of y)
    lines(pred_res_df$Xfact,pred_res_df$total.yield,
          col='dodgerblue',lwd=1.6)
    par(oma = c(0, 0, 0, 0), new = TRUE)
    legend("top", c("value", "yield", "biomass"), xpd = TRUE,
           horiz = TRUE, inset = c(0, -0.1), bty = "n",lty = 1,seg.len = 0.7,
           col = c('darkorange','dodgerblue','darkgreen'), cex = 0.8,lwd=2,
           text.width=0.3,x.intersp=0.3)
    plot1 <- recordPlot()

    #prediction by fleet
    if(length(fleet_mat) != 1){
      df.TBfleet <- df.TBpred

      #correct fleet_mat if provided in Catch not fishing mortalities
      if(fleet_unit == 'Catch'){
        total.catch <- sum(df.TBfleet$num.caught,na.rm = TRUE)
        total.FM <- sum(df.TBfleet$FM,na.rm = TRUE)
        fleet_mat.FM <- (total.catch * fleet_mat ) / total.FM
      } #fleet_mat.Catch needed? ??? ? ? ?? ? ? ?? ??????????????????????????
      if(fleet_unit == 'FM'){
        fleet_mat.FM <- fleet_mat
      }

      #if no fleet_FM_change matrix is provided FM_change is used, fleet_FM_change
      #matrix provides a detailed manual changeable matrix with specific FM changes
      #per fleet, FM_change is just a standard assumed vector varying FM for both
      #fleets simultaneously
      if(length(fleet_FM_change) == 1){
        pred_fleet_mat_list <- list()
        for(x8 in 1:dim(fleet_mat)[2]){
          pred.fleet_mat <- as.matrix(fleet_mat[,x8]) %*% FM_change
          colnames(pred.fleet_mat) <- FM_change

          df.TBfleet.loop <- df.TBfleet
          pred.fleet_res_list <- list()
          for(x9 in 1:length(FM_change)){
            df.TBfleet.loop$FM <- pred.fleet_mat[,x9]
            df.TBfleet.loop$Z <- pred.fleet_mat[,x9] + nM
            res <- calc_TB.age(df.TB = df.TBfleet.loop, unit.time = unit.time, stock_size_1 = stock_size_1)
            pred.fleet_res_list[[x9]] <- res$totals
          }

          pred.fleet_res_df <- do.call(rbind, pred.fleet_res_list)
          pred.fleet_res_df$Xfact <- FM_change
          rownames(pred.fleet_res_df) <- FM_change

          pred_fleet_mat_list[[x8]] <- pred.fleet_res_df
        }

        #dimensions of plotting area
        max_val <- round(max(pred_res_df$total.value,na.rm=TRUE),digits=0)
        dim_val <- 10 ^ (nchar(max_val)-1)
        max_yiel <- round(max(pred_res_df$total.yield,na.rm=TRUE),digits=0)
        dim_yiel <- 10 ^ (nchar(max_yiel)-1)
        max_bio <- round(max(pred_res_df$mean.biomass,na.rm=TRUE),digits=0)
        dim_bio <- 10 ^ (nchar(max_bio)-1)

        #4 panel plot
        par(new=FALSE,mar = c(4, 4, 4, 4) + 0.1, mfrow=c(2,2))
        plot(pred_res_df$Xfact,pred_res_df$total.value, type ='o',main='Value',
             col ='darkorange', ylim = c(0,ceiling(max_val/dim_val)*dim_val),
             lwd=1.6, xlab = "F-factor X", ylab = "Value",lty=1)
        abline(v=1,col='grey',lty=2)
        for(j in 1:length(pred_fleet_mat_list)){
          lines(pred_fleet_mat_list[[j]]$Xfact,pred_fleet_mat_list[[j]]$total.value,
                col ='darkorange',lwd=1.6,lty = (j+1),type='o')
        }
        plot(pred_res_df$Xfact,pred_res_df$total.yield,
             ylim = c(0,ceiling(max_yiel/dim_yiel)*dim_yiel),
             lty = 1, col='dodgerblue',lwd=1.6,main='Yield',
             type='o',ylab='Yield',xlab='F-factor X')
        abline(v=1,col='grey',lty=2)
        for(j in 1:length(pred_fleet_mat_list)){
          points(pred_fleet_mat_list[[j]]$Xfact,pred_fleet_mat_list[[j]]$total.yield,
                 lty = (j+1), col='dodgerblue',lwd=1.6,type='o')
        }
        plot(pred_res_df$Xfact,pred_res_df$mean.biomass,
             ylim = c(0,ceiling(max_bio/dim_bio)*dim_bio),
             lty = 1, col = 'darkgreen',lwd=1.6,main='Biomass',
             type='o',ylab='Biomass',xlab='F-factor X')
        abline(v=1,col='grey',lty=2)
        for(j in 1:length(pred_fleet_mat_list)){
          points(pred_fleet_mat_list[[j]]$Xfact,pred_fleet_mat_list[[j]]$mean.biomass,
                 lty = (j+1), col = 'darkgreen',lwd=1.6,type='o')
        }
        plot(1,axes=F,xlab='',ylab='',type='n')
        legend("center", c("total",colnames(fleet_mat)),bty = "n",
               lty = 1:(length(pred_fleet_mat_list)+1), seg.len = 0.2,
               col = 'black',lwd=2,xpd=TRUE,y.intersp = 0.15,
               text.width=0.3,x.intersp=0.3)
        plot2 <- recordPlot()
      }


      #if fleet_FM_change matrix is provided
      if(length(fleet_FM_change) != 1){

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

          df.TBfleet.loop <- df.TBfleet[,1:4]
          df.TBfleet.loop$FM <- FM_fleets.df$FMtot
          df.TBfleet.loop$Z <- FM_fleets.df$FMtot + nM

          res <- calc_TB.age(df.TB = df.TBfleet.loop, unit.time = unit.time, stock_size_1 = stock_size_1)
          FM_fleets.df$Ctot <- res$dfs$num.caught
          FM_fleets.df$w <- res$dfs$mean_weight
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
              (FM_fleets.list[[x11]] * res$dfs$dt)
          }
          B_fleets.df <- do.call(cbind,B_fleets.list)
          B_fleets.df <- as.data.frame(B_fleets.df)
          B_fleets.df$Btot <- rowSums(B_fleets.df,na.rm=TRUE)

          Btots <- colSums((B_fleets.df * res$dfs$dt),na.rm=TRUE) /
            sum(res$dfs$dt,na.rm=TRUE)

          #VALUE
          V_fleets.list <- list()
          for(x12 in 1:length(FM_fleets.list)){
            V_fleets.list[[x12]] <- Y_fleets.list[[x12]] * res$dfs$value
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
        par(new=FALSE,mar = c(4, 4, 4, 4) + 0.1, mfrow=c(2,2))
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
        plot2 <- recordPlot()
      }
    }

    #if different Lcs are provided
    if(!is.na(Lc_mat)){

      Lt <- Linf * (1- exp(-K * (df.TBpred$classes.num - t0)))


      if(s_list$selecType == 'trawl_ogive'){
        sel <- 1 / (1 + exp(- (Lt - L50)/
                              ((2*(L75-L50))/(log(0.75/(1-0.75))-
                                                log(0.25/(1-0.25))))))

        #         #alternative:
        #         sel <- 1 / (1 + exp(S1.TS - S2.TS * Lt))


      }else if(s_list$selecType == 'gillnet'){
        if (s_list$selecDist == "lognormal") {
          sel <- (1/df.TBpred$classes.num) * exp(s_list$select_p1 +
                                                   log(s_list$mesh_size/s_list$mesh_size1) -
                                                   (s_list$select_p2^2/2) -
                                                   (((log(df.TBpred$classes.num) -
                                                        s_list$select_p1 -
                                                        log(s_list$mesh_size/s_list$mesh_size1))^2)/
                                                      (2 *s_list$select_p2^2)))
        }
        if (s_list$selecDist == "normal_fixed") {
          sel <- exp(-((df.TBpred$classes.num - s_list$mesh_size *
                          s_list$select_p1)^2/
                         (2 * s_list$select_p2^2)))
        }
      }



    }
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                       Length data                        #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if('midLengths' %in% names(res)){}
} ## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?

