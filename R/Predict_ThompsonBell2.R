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

data("data_Predict_ThompsonBell")
param = data_Predict_ThompsonBell
unit.time = "year"
stock_size_1 = NA
plus.group = NA
FM_change = seq(1,10,1)


#param can be age or midLenghts



Predict_ThompsonBell2 <- function(param, FM_change,

                                 select.param = NA, unit.time = "year",
                                 stock_size_1 = NA, plus.group = NA,


   Tr,stock_size_1 = NA,
   FM_change = NA,
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

    #prediction based on f_change
    pred_mat <- as.matrix(FM) %*% FM_change

    pred_res_list <- list()
    for(x7 in 1:length(FM_change)){
      param$Z <- pred_mat[,x7] + nM
      param$FM <- pred_mat[,x7]
      res <- Predict_ThompsonBell1(param)
      pred_res_list[[x7]] <- res$totals
    }

    pred_res_df <- do.call(rbind, pred_res_list)
    pred_res_df$Xfact <- FM_change

    #save x axis positions
    maxY <- round(max(pred_res_df$tot.Y,na.rm=TRUE),digits=0)
    maxV <- round(max(pred_res_df$tot.V,na.rm=TRUE),digits=0)
    maxB <- round(max(pred_res_df$meanB,na.rm=TRUE),digits=0)
    #####HERE

    maxis <- c(maxY,maxV,maxB)
    which(maxis == max(maxis))
    max_val <- round(max(pred_res_df$tot.V,na.rm=TRUE),digits=0)
    dim_val <- 10 ^ (nchar(max_val)-1)

    par(oma = c(1, 1, 1.5, 1),new=FALSE,mar = c(5, 4, 4, 4) + 0.3)
    plot(pred_res_df$Xfact,pred_res_df$tot.V, type ='l',
         col ='darkorange', ylim = c(0,ceiling(max_val/dim_val)*dim_val),
         lwd=1.6, xlab = "F-factor X", ylab = "Value")
    lines(pred_res_df$Xfact,pred_res_df$meanB,
          col = 'darkgreen',lwd=1.6)    # draw lines with small intervals: seq(0,max(),0.05) but y as to be dependent of x (formula of calculaiton of y)
    lines(pred_res_df$Xfact,pred_res_df$tot.Y,
          col='dodgerblue',lwd=1.6)
    par(oma = c(0, 0, 0, 0), new = TRUE)
    legend("top", c("value", "yield", "biomass"), xpd = TRUE,
           horiz = TRUE, inset = c(0, -0.1), bty = "n",lty = 1,seg.len = 0.7,
           col = c('darkorange','dodgerblue','darkgreen'), cex = 0.8,lwd=2,
           text.width=0.3,x.intersp=0.3)

    ret <- c(res,res2)
    return(ret)

  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                       Length data                        #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if('midLengths' %in% names(res)){}
} ## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?

