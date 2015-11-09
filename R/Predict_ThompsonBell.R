#' @title Thompson and Bell prediction model
#
#' @description  This is a function to calculate the total mortality (Z) from length composition data via the length converted catch curve or from age at length data with catch curve.
#'
#' @param classes Midpoints of the length class as vector (length frequency data) or ages as vector (age composition data).
#' @param catch Catch as vector, or a matrix with catches of subsequent years if the catch curve with constat time intervals should be applied.
#' @param datatype Type of data which is used for analysis, either 'length' or 'age', for length frequency or age composition data, respectively
#' @param Linf Infinite length for investigated species in cm [cm].
#' @param K Growth coefficent for investigated species per year [1/year].
#' @param t0 Theoretical time zero, at which individuals of this species hatch (default: 0).
#' @param Winf Infinite weight in grams!
#'
#' @param FM reference F-at-age-array
#'
#' @examples
#' data("ex.Predict_ThompsonBell")
#' with(ex.Predict_ThompsonBell,Predict_ThompsonBell(classes = age,
#'   mean_weight = meanWeight,value = valueGramm,FM = FM,Z = Z,
#'   datatype = 'age',Tr = 1,stock_size_1 = NA,plus.group = 12))
#'
#' @details better to treat last group always as a plus group..... For variable parameter system vectors are reuqired for constant parameter systems matrices or data.frames have to be inserted. or vectors The length converted linearised catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @references
#' example 1 : Kuwait (Garcia and van zalinge 1982)
#'
#' @export


data("ex.Predict_ThompsonBell")
classes = ex.Predict_ThompsonBell$age
mean_weight = ex.Predict_ThompsonBell$meanWeight
value = ex.Predict_ThompsonBell$valueGramm
FM = ex.Predict_ThompsonBell$FM
Z = ex.Predict_ThompsonBell$Z
datatype = 'age'
Tr = 1
stock_size_1 = NA
unit.time = 'month'
plus.group = 12
FM_change = seq(0,3,0.2)
fleet_mat=matrix(ncol = 2,nrow = 12)
fleet_mat[,1] <- c(0.72,0.96,0.84,0.48,0.6,0.48,1.08,0.48,0.084,0.12,0.24,0.24)
fleet_mat[,2] <- c(0.48,0.36,0.48,0.96,1.32,0.72,0.48,0.72,1.116,1.68,2.52,2.28)
colnames(fleet_mat) <- c("artisanal", "industrial")
fleet_unit = "FM"   #or "Catch"
fleet_FM_change = matrix(ncol = 2,nrow = 8)
fleet_FM_change[,1] <- rep(1,8)
fleet_FM_change[,2] <- c(0,0.4,0.8,1.0,1.2,1.5,2,3)
colnames(fleet_FM_change) <- c("artisanal", "industrial")



Predict_ThompsonBell <- function(classes, mean_weight, FM, Z, value = NA, Tr,
                                 datatype, stock_size_1 = NA, plus.group = NA,
                                 FM_change = NA, fleet_mat = NA,fleet_FM_change = NA){

  df.TB <- cbind(classes,mean_weight,value,FM,Z)
  df.TB <- as.data.frame(df.TB)
  df.TB$classes <- as.character(df.TB$classes)

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(df.TB$classes, split="\\+"))
  df.TB$classes.num <- as.numeric(classes.num[,1])

  #natural Mortality
  nM <- mean(df.TB$Z - df.TB$FM,na.rm=T)

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                        Age data                          #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if(datatype == 'age'){


    calc_TB.age <- function(df.TB,unit.time,stock_size_1){
      # delta t
      df.TB$dt <- NA
      for(x0 in 1:(length(df.TB$dt)-1)){
        df.TB$dt[x0] <- df.TB$classes.num[x0+1] - df.TB$classes.num[x0]
      }
      if(unit.time == 'month'){
        df.TB$dt <- df.TB$dt * 1/12
      }

      #population per age group
      df.TB$population <- NA
      if(!is.na(stock_size_1)){
        df.TB$population[1] <- stock_size_1
      }else df.TB$population[1] <- 1000  #if not provided results are just negative

      for(x1 in 2:length(df.TB$population)){
        df.TB$population[x1] <- df.TB$population[x1-1] * exp(-df.TB$Z[x1]*df.TB$dt[x1])
        if(x1 == length(df.TB$population)){           ####CORRECT?????
          df.TB$population[x1] <- df.TB$population[x1-1] * exp(-df.TB$Z[x1-1]*
                                                                 df.TB$dt[x1-1])
        }
      }

      #number of deaths per time step month or year
      df.TB$num.deaths <- NA
      for(x2 in 1:length(df.TB$num.deaths)){
        df.TB$num.deaths[x2] <- df.TB$population[x2] - df.TB$population[x2+1]
      }

      #catch in numbers
      df.TB$num.caught <- df.TB$num.deaths * df.TB$FM / df.TB$Z

      #yield in kg or g
      df.TB$yield <- df.TB$num.caught * df.TB$mean_weight

      #mean biomass in kg or g
      df.TB$biomass <- df.TB$yield / (df.TB$FM * df.TB$dt)

      #value expressed in money units
      df.TB$value_yield <- df.TB$yield * df.TB$value

      #total yield, value and average biomass
      total.yield <- sum(df.TB$yield,na.rm=TRUE)
      total.value <- sum(df.TB$value_yield,na.rm=TRUE)
      mean.biomass <- sum((df.TB$biomass * df.TB$dt),na.rm=TRUE) /
        sum(df.TB$dt,na.rm=TRUE)   ### more complicated biomass concept if dt is not constant, see Chapter 5

      res_funct <- data.frame(total.yield = total.yield,
                 mean.biomass = mean.biomass,
                 total.value = total.value)
      return(list(totals = res_funct,dfs = df.TB))

    }


    TB.results <- calc_TB.age(df.TB = df.TB, unit.time = unit.time, stock_size_1 = stock_size_1)
    df.TB <- TB.results$dfs

    #with plus group
    if(!is.na(plus.group)){
      #new dataframe
      df.TBnew <- df.TB[1:plus.group,]
      #new class
      df.TBnew$classes[length(df.TBnew$classes)] <-
        paste(df.TBnew$classes[length(df.TBnew$classes)],"+",sep='')
      #num deaths
      df.TBnew$num.deaths[length(df.TBnew$num.deaths)] <-
        df.TBnew$population[plus.group]
      #catch
      df.TBnew$num.caught[plus.group] <- (df.TBnew$FM[plus.group]/
                                            df.TBnew$Z[plus.group]) *
        df.TBnew$population[plus.group]
      catch.plus.dif <- df.TBnew$num.caught[plus.group] - df.TB$num.caught[plus.group]
      #yield
      df.TBnew$yield[plus.group] <- df.TBnew$mean_weight[plus.group] * catch.plus.dif
      #value
      df.TBnew$value_yield[plus.group] <-
        df.TBnew$yield[plus.group] * df.TBnew$value[plus.group]
      #biomass       ####not sure....omitted in manual
      df.TBnew$biomass[plus.group] <-
        df.TBnew$yield[plus.group] / (df.TBnew$FM[plus.group] * df.TBnew$dt[plus.group])
    }

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
           horiz = TRUE, inset = c(0, -0.1), bty = "n",lty = 1,seg.len = 0.5,
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
      }
      if(fleet_unit == 'FM'){
        fleet_mat.FM <- fleet_mat
      }

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
        plot(pred_res_df$Xfact,pred_res_df$total.value, type ='l',main='Value',
             col ='darkorange', ylim = c(0,ceiling(max_val/dim_val)*dim_val),
             lwd=1.6, xlab = "F-factor X", ylab = "Value",lty=1)
        abline(v=1,col='grey',lty=2)
        for(j in 1:length(pred_fleet_mat_list)){
          lines(pred_fleet_mat_list[[j]]$Xfact,pred_fleet_mat_list[[j]]$total.value,
                col ='darkorange',lwd=1.6,lty = (j+1))
        }
        plot(pred_res_df$Xfact,pred_res_df$total.yield,
             ylim = c(0,ceiling(max_yiel/dim_yiel)*dim_yiel),
             lty = 1, col='dodgerblue',lwd=1.6,main='Yield',
             type='l',ylab='Yield',xlab='F-factor X')
        abline(v=1,col='grey',lty=2)
        for(j in 1:length(pred_fleet_mat_list)){
          lines(pred_fleet_mat_list[[j]]$Xfact,pred_fleet_mat_list[[j]]$total.yield,
                lty = (j+1), col='dodgerblue',lwd=1.6)
        }
        plot(pred_res_df$Xfact,pred_res_df$mean.biomass,
             ylim = c(0,ceiling(max_bio/dim_bio)*dim_bio),
             lty = 1, col = 'darkgreen',lwd=1.6,main='Biomass',
             type='l',ylab='Biomass',xlab='F-factor X')
        abline(v=1,col='grey',lty=2)
        for(j in 1:length(pred_fleet_mat_list)){
          lines(pred_fleet_mat_list[[j]]$Xfact,pred_fleet_mat_list[[j]]$mean.biomass,
                lty = (j+1), col = 'darkgreen',lwd=1.6)
        }
        plot(1,axes=F,xlab='',ylab='',type='n')
        legend("center", c("total",colnames(fleet_mat)),bty = "n",
               lty = 1:(length(pred_fleet_mat_list)+1), seg.len = 0.2,
               col = 'black',lwd=2,xpd=TRUE,y.intersp = 0.15,
               text.width=0.3,x.intersp=0.3)
        plot2 <- recordPlot()
      }


      if(length(fleet_FM_change) != 1){

        pred_fleet_mat_list <- list()
        for(x8 in 1:dim(fleet_mat)[2]){
          pred.fleet_mat <- as.matrix(fleet_mat[,x8]) %*% fleet_FM_change[,x8]
          colnames(pred.fleet_mat) <- fleet_FM_change[,x8]

          df.TBfleet.loop <- df.TBfleet
          pred.fleet_res_list <- list()
          for(x9 in 1:length(fleet_FM_change[,x8])){
            df.TBfleet.loop$FM <- pred.fleet_mat[,x9]
            df.TBfleet.loop$Z <- pred.fleet_mat[,x9] + nM
            res <- calc_TB.age(df.TB = df.TBfleet.loop, unit.time = unit.time, stock_size_1 = stock_size_1)
            pred.fleet_res_list[[x9]] <- res$totals
          }

          pred.fleet_res_df <- do.call(rbind, pred.fleet_res_list)
          pred.fleet_res_df$Xfact <- fleet_FM_change[,x8]
          #rownames(pred.fleet_res_df) <- fleet_FM_change[,x8]

          pred_fleet_mat_list[[x8]] <- pred.fleet_res_df
        }

        pred_fleet_mat_list




      }


    }
  }


  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                       Length data                        #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if(datatype == 'length'){
    }
} ## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?

