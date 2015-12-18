#' @title Prediction models: YPR and Thompson and Bell model
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
#' @param FM_change vector with ascending Fishing mortalities
#' @param Lc_tc_change vector with ascending Lc values
#' @param type indicating which model should applied: \code{"ypr"} for Beverton and Holt's
#'   yield per recruit model and \code{"ThompBell"} for the Thompson and Bell model
#'
#' @keywords function prediction
#'
#' @examples ...
#'
#' @details better to treat last group always as a plus group.....
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{tplusdt_2} or \strong{t_midL}: relative,
#'   \item \strong{lnC_dt}: rearranged,
#'   \item \strong{reg_int}: the,
#'   \item \strong{Z}: the,
#'   \item \strong{se}: the;}
#' }
#'
#' @references
#' example 1 : Kuwait (Garcia and van Zalinge 1982)
#'
#' Millar, R. B., & Holst, R. (1997). Estimation of gillnet and hook selectivity using
#' log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3):471-477
#'
#' #@export

predict_mod <- function(param, FM_change, Lc_tc_change, type){


  # Beverton and Holt's ypr
#------------

  res <- param
  M <- res$M
  K <- res$K
  t0 <- ifelse(!is.null(res$t0),res$t0,0)
  a <- ifelse(!is.null(res$a),res$a,NA)
  b <- ifelse(!is.null(res$b),res$b,NA)

  if(length(FM) == 1 & is.na(FM[1])){
    FM <- seq(0,10,0.1)
    print(noquote("No fishing mortality (FM) was provided, a default range of 0 to 10 is used."))
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                        Age data                          #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if("Winf" %in% names(res)){
    Winf <- res$Winf
    tr <- res$tr
    tc <- tc_Lc
    list_tc_runs <- list()
    for(i in 1:length(tc)){
      tci <- tc[i]

      #Biomass per Recruit
      S <- exp(-K * (tci - t0))
      B_R <- exp(-M*(tci-tr)) * Winf *
        ((1/(FM+M)) - ((3*S)/((M+FM)+K)) +
           ((3*(S^2))/((M+FM)+(2*K))) - ((S^3)/((M+FM) + (3*K))))

      #virgin biomass
      Bv_R <- B_R[which(FM == 0)]
      #biomass of exploited part of the cohort (biomass of fish older than tc)

      #biomass in percetage of virgin biomass
      B_R.percent <- round((B_R / Bv_R ) * 100, digits = 1)

      #Yield per Recruit
      Y_R <- B_R * FM

      #relative yield per recruit - mostly done with length frequency data (exclusively?)
      Y_R.rel <- Y_R * (exp(-M *(tr - t0))) / Winf


      #mean age in annual yield
      Ty <- (1 / (M+FM)) + tci

      #mean length in the annual yield
      #Ly <- Linf * (1 - (((M+FM)*S)/((M+FM)+K)))

      #mean weight in annual yield
      Wy <- (M+FM) * Winf *
        ((1/(FM+M)) - ((3*S)/((M+FM)+K)) +
           ((3*(S^2))/((M+FM)+(2*K))) - ((S^3)/((M+FM) + (3*K))))


      results.PBH <- data.frame(FM = FM,
                                Y_R = Y_R,
                                Y_R.rel = Y_R.rel,
                                B_R = B_R,
                                B_R.percent = B_R.percent,
                                Ty = Ty,
                                Wy = Wy)


      list_tc_runs[[i]] <- results.PBH

    }
    names(list_tc_runs) <- tc
    ret <- list(res,
                FM = FM,
                tc = tc,
                list_tc_runs = list_tc_runs)
    class(ret) <- "ypr_mod"
    # plot results
    plot(ret)
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                       Length data                        #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if("Linf" %in% names(res)){
    Linf <- res$Linf
    Lr <- res$Lr
    Lc <- tc_Lc
    list_Lc_runs <- list()
    for(i in 1:length(Lc)){
      Lci <- Lc[i]

      E <- FM/(FM + M)
      # transform Linf in Winf #### CHECK!!
      Winf <- a * (Linf ^ b)

      #yiel per recurit for length data   # also possbile inputing option: F/K
      S <- 1 - (Lci/Linf)  # == U  ##(U <- 1 - (Lci/Linf))
      A <- ((Linf - Lci)/(Linf-Lr))^(M/K)
      Y_R <- FM * A * Winf * ((1/(M+FM)) - (3*S)/((M+FM)+K) +
                                (3*S^2)/((M+FM)+2*K) -
                                (S^3)/((M+FM)+3*K))

      #biomass per recruit for length data?
      B_R <- Y_R / FM

      #virgin biomass
      Bv_R <- B_R[which(FM == 0)]
      #biomass of exploited part of the cohort (biomass of fish older than tc)

      #biomass in percetage of virgin biomass
      B_R.percent <- round((B_R / Bv_R ) * 100, digits = 1)

      #relative yield per recruit - mostly done with length frequency data (exclusively?)
      m <- K/(M+FM)
      Y_R.rel <- E * S^(M/K) * (1 - ((3*S)/(1+m)) +
                                  ((3*S^2)/(1+2*m)) - ((S^3)/(1+3*m)))

      #mean length in the annual yield
      Ly <- Linf * (1 - (((M+FM)*S)/((M+FM)+K)))

      #mean weight in annual yield
      # Wy <- (M+FM) * Winf *
      #    ((1/(FM+M)) - ((3*S)/((M+FM)+K)) +
      #       ((3*(S^2))/((M+FM)+(2*K))) - ((S^3)/((M+FM) + (3*K))))


      results.PBH <- data.frame(FM = FM,
                                Ly = Ly,
                                E = E,
                                Y_R = Y_R,
                                Y_R.rel = Y_R.rel,
                                B_R = B_R,
                                B_R.pecent = B_R.percent)


      list_Lc_runs[[i]] <- results.PBH

    }
    names(list_Lc_runs) <- Lc
    ret <- list(res,
                FM = FM,
                Lc = Lc,
                list_Lc_runs = list_Lc_runs
    )
    class(ret) <- "ypr_mod"

    # plot results
    plot(ret)
  }
#------------


  # Thompson and Bell model
#------------
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
#------------


}





## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?

