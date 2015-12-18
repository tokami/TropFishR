#' @title Prediction models: Thompson and Bell model
#
#' @description  This is a function ...
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{age} or \strong{midLengths}: midpoints of the length class as vector (length-frequency
#'   data) or ages as vector (age composition data),
#'   \item \strong{Linf}: infinite length for investigated species in cm [cm],
#'   \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'   \item \strong{t0}: theoretical time zero, at which individuals of this species hatch,
#'   \item \strong{catch}: catch as vector, or a matrix with catches of subsequent years if
#'   the catch curve with constat time intervals should be applied;
#' }
#' @param FM_change Vector with ascending Fishing mortalities
#'
#' @keywords function prediction
#'
#' @examples
#' # yield per recruit with age structured data
#' # load data
#' data(shrimps)
#'
#' # run model
#' predict_mod(shrimps,seq(0.1,3,0.1))
#'
#' # yield per recruit with length structured data
#' # load data
#' data(hake)
#'
#' # run model
#' predict_mod(hake,seq(0.1,3,0.1))
#'
#' @details better to treat last group always as a plus group..... For variable parameter system vectors are reuqired for constant parameter systems matrices or data.frames have to be inserted. or vectors The length converted linearised catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{tplusdt_2} or \strong{t_midL}: relative,
#'   \item \strong{lnC_dt}: rearranged,
#'   \item \strong{reg_int}: the,
#'   \item \strong{Z}: the,
#'   \item \strong{se}: the;}
#' in case of calc_ogive, additionally:
#' \itemize{
#'   \item \strong{intercept}: intercep,
#'   \item \strong{Sobs}: observed,
#'   \item \strong{ln_1_S_1}: dependent,
#'   \item \strong{Sest}: estimated,
#'   \item \strong{t50}: age,
#'   \item \strong{t75}: age,
#'   \item \strong{L50}: length,
#'   \item \strong{L75}: length;
#' }
#'
#' @references example 1 : Kuwait (Garcia and van zalinge 1982)
#'   Millar, R. B., & Holst, R. (1997). Estimation of gillnet and hook selectivity using log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
#' @export

predict_mod <- function(param, FM_change){

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

#   FM <- res$FM
#   Z <- res$Z
#   #natural Mortality
#   nM <- mean(Z - FM,na.rm=T)

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

  #prediction based on f_change
  pred_mat <- as.matrix(FM) %*% FM_change

  pred_res_list <- list()
  for(x7 in 1:length(FM_change)){
    param$Z <- pred_mat[,x7] + nM
    param$FM <- pred_mat[,x7]
    res <- stock_sim(param)
    pred_res_list[[x7]] <- res$totals
  }

  pred_res_df <- do.call(rbind, pred_res_list)
  pred_res_df$Xfact <- FM_change

  #save x axis positions
  max_val <- round(max(pred_res_df$tot.V,na.rm=TRUE),digits=0)
  dim_val <- 10 ^ (nchar(max_val)-1)
  max_yiel <- round(max(pred_res_df$tot.Y,na.rm=TRUE),digits=0)
  dim_yiel <- 10 ^ (nchar(max_yiel)-1)
  max_bio <- round(max(pred_res_df$meanB,na.rm=TRUE),digits=0)
  dim_bio <- 10 ^ (nchar(max_bio)-1)

#   max_val <- round(max(pred_res_df$tot.V,na.rm=TRUE),digits=0)
#   dim_val <- 10 ^ (nchar(max_val)-1)

  op <- par(oma = c(1, 1, 1.5, 1),new=FALSE,mar = c(5, 4, 4, 6) + 0.3)
  plot(pred_res_df$Xfact,pred_res_df$tot.V, type ='o',ylab='Value',xlab='F-factor X',
       col ='darkorange', ylim = c(0,ceiling(max_val/dim_val)*dim_val),
       lwd=1.6)
  par(new=TRUE)
  plot(pred_res_df$Xfact,pred_res_df$tot.Y,type ='o',ylab='',xlab='',
       col='dodgerblue',lwd=1.6,axes=FALSE,
       ylim = c(0,ceiling(max_yiel/dim_yiel)*dim_yiel))
  axis(4,at=pretty(c(0,pred_res_df$tot.Y)))
  mtext("Yield", side=4, line=2)
  par(new=TRUE)
  plot(pred_res_df$Xfact,pred_res_df$meanB,type='o',axes=FALSE,ylab='',xlab='',
       col = 'darkgreen',lwd=1.6,
       ylim = c(0,ceiling(max_bio/dim_bio)*dim_bio))    # draw lines with small intervals: seq(0,max(),0.05) but y as to be dependent of x (formula of calculaiton of y)
  axis(4,at=pretty(c(0,pred_res_df$meanB)),line = 3)
  mtext("Biomass", side=4, line=5)

  par(oma = c(0, 0, 0, 0), new = TRUE)
  legend("top", c("value", "yield", "biomass"), xpd = TRUE,
         horiz = TRUE, inset = c(0, -0.1), bty = "n",lty = 1,seg.len = 0.7,
         col = c('darkorange','dodgerblue','darkgreen'), cex = 0.8,lwd=2,
         text.width=0.3,x.intersp=0.3)
  par(op)

  res2 <- pred_res_df
  ret <- c(res,res2)
  return(ret)
} ## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?

