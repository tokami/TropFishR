#' @title ELEFAN
#'
#' @description Electronic LEngth Frequency ANalysis for estimating growth parameter of exploited fish populations.
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row) and sampling date (column);
#' }
#' @param Linf_fix numeric; if used the K-Scan method is applied with a fixed Linf value (i.e. varying K only).
#' @param Linf_range numeric vector with potential Linf values. Default is the last length class plus/minus 5 cm
#' @param K_range K values for which the score of growth functions should be calculated
#'    (by default: exp(seq(log(0.1),log(10),length.out = 100)))
#' @param C growth oscillation amplitude (default: 0)
#' @param WP winter point (default: 0)
#' @param MA number indicating over how many length classes the moving average should be performed (defalut: 5, for
#'    more information see \link{lfqRestructure}).
#' @param addl.sqrt Passed to \link{lfqRestructure}. Applied an additional square-root transformation of positive values according to Brey et al. (1988).
#'    (default: FALSE, for more information see \link{lfqRestructure}).
#'
#' @examples
#' data(trout)
#'
#' # K-Scan
#' ELEFAN(trout, Linf_fix = 16)
#'
#' # Surface response analysis
#' ELEFAN(trout, K_range = seq(0.1,2,0.1), Linf_range = seq(14,17,1))
#'
#' @details This functions allows to perform the K-Scan and Response surface analysis to estimate growth parameters.
#'    It combines the step of restructuring length-frequency data (\link{lfqRestructure}) followed by the fitting of VBGF
#'    curves through the restructured data (\link{lfqFitCurves}). K-Scan is a method used to search for the K
#'    parameter with the best fit while keeping the Linf fixed. In contrast, with response surface analysis
#'    both parameters are estimated and the fits are displayed in a heatmap.
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline axis grid image mtext par plot text title
#'
#' @references
#' Brey, T., Soriano, M., and Pauly, D. 1988. Electronic length frequency analysis: a revised and expanded
#' user's guide to ELEFAN 0, 1 and 2.
#'
#' Pauly, D. 1981. The relationship between gill surface area and growth performance in fish:
#' a generalization of von Bertalanffy's theory of growth. \emph{Meeresforsch}. 28:205-211
#'
#' Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the objective extraction of
#' growth parameters from length-frequency data. \emph{Meeresforschung}, 28(4):205-211
#'
#' Pauly, D., 1985. On improving operation and use of ELEFAN programs. Part I: Avoiding
#' "drift" of K towards low values. \emph{ICLARM Conf. Proc.}, 13-14
#'
#' Pauly, D., 1987. A review of the ELEFAN system for analysis of length-frequency data in
#' fish and aquatic invertebrates. \emph{ICLARM Conf. Proc.}, (13):7-34
#'
#' Pauly, D. and G. R. Morgan (Eds.), 1987. Length-based methods in fisheries research.
#' (No. 13). WorldFish
#'
#' Pauly, D. and G. Gaschuetz. 1979. A simple method for fitting oscillating length growth data, with a
#' program for pocket calculators. I.C.E.S. CM 1979/6:24. Demersal Fish Cttee, 26 p.
#'
#' Pauly, D. 1984. Fish population dynamics in tropical waters: a manual for use with programmable
#' calculators (Vol. 8). WorldFish.
#'
#' Quenouille, M. H., 1956. Notes on bias in estimation. \emph{Biometrika}, 43:353-360
#'
#' Somers, I. F., 1988. On a seasonally oscillating growth function. ICLARM Fishbyte 6(1): 8-11.
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2): 407 p.
#'
#' Tukey, J., 1958. Bias and confidence in not quite large samples.
#' \emph{Annals of Mathematical Statistics}, 29: 614
#'
#' Tukey, J., 1986. The future of processes of data analysis. In L. V. Jones (Eds.),
#' The Collected Works of John W. Tukey-philosophy and principles of data analysis:
#' 1965-1986 (Vol. 4, pp. 517-549). Monterey, CA, USA: Wadsworth & Brooks/Cole
#'
#' @export

ELEFAN <- function(param, Linf_fix = NA, Linf_range = NA,
                   K_range = exp(seq(log(0.1),log(10),length.out=100)),
                   C = 0, WP = 0,
                   MA = 5, addl.sqrt = FALSE){

  res <- param
  classes <- res$midLengths
  catch <- res$catch
  dates <- res$dates
  interval <- (classes[2]-classes[1])/2

  n_samples <- dim(catch)[2]
  n_classes <- length(classes)

  # if(is.na(tmax)) tmax <- n_classes - 1   # minus 1 because one cohort is start cohort, the others are cohorts which are added to the start cohort and together they should not exceed n_classes
  # n_cohorts <- tmax # + 1 cohort for each year of sampling period (e.g. if firs and last sampling 1.5 years = + 1 cohort) # a combination of: length of sampling period (per year one new cohort), and tmax (surviving cohorts from last years)
  # n_cohorts <- 3

  if(is.na(Linf_fix) & is.na(Linf_range[1])) Linf_range <- seq(classes[n_classes]-5,classes[n_classes]+5,1) ### OLD: c(classes[n_classes]-5,classes[n_classes]+5)

  # ELEFAN 0
  res <- lfqRestructure(res, MA = MA, addl.sqrt = addl.sqrt)
  catch_aAF_F <- res$rcounts
  peaks_mat <- res$peaks_mat
  ASP <- res$ASP

  if(!is.na(Linf_fix)){
    Linfs <- Linf_fix
  }else Linfs <-  Linf_range ## OLD: seq(Linf_range[1],Linf_range[2],Linf_step)
  Ks <- K_range #OLD: seq(K_range[1],K_range[2],K_step)

  ESP_starting_point_L <- array(NA,dim=c(length(Ks),2,length(Linfs)))
  ESP_list_L <- matrix(NA,nrow=length(Ks),ncol=length(Linfs))
  for(li in 1:length(Linfs)){

    ESP_starting_point_K <- matrix(NA,nrow=length(Ks),ncol=2)
    ESP_list_K <- rep(NA,length(Ks))
    for(ki in 1:length(Ks)){

      # ELEFAN 1
      paras <- list(Linf=Linfs[li],K = Ks[ki], C=C, WP=WP)
      ele1_res <- lfqFitCurves(lfq = res, par = paras)

      max_ESP <- ele1_res$max_ESP
      startingSample <- ele1_res$startingSample
      startingLength <- ele1_res$startingLength

      # export max_ESP, startingSample, startingClass
      ESP_list_K[ki] <- max_ESP
      ESP_starting_point_K[ki,] <- c(startingSample,startingLength)
    }
    ESP_list_L[,li] <- ESP_list_K
    ESP_starting_point_L[,,li] <- ESP_starting_point_K
  }

  dimnames(ESP_list_L) <- list(Ks,Linfs)
  score_mat <- round((10^(ESP_list_L/ASP)) /10,digits = 3)
  rownames(score_mat) <- round(as.numeric(rownames(score_mat)),digits = 2)

  colnames(ESP_starting_point_L) <- c("startingSample","startingLength")

  # Graphs
  if(is.na(Linf_fix)){
    plot_dat <- reshape2::melt(score_mat)
    image(x = Linfs,
          y = Ks,
          z = t(score_mat), col=colorRampPalette(c("yellow","red"), space="Lab")(5),
          main = 'Response surface analysis', ylab = 'K', xlab='Linf')
    #grid (NULL,NULL, lty = 6, col = "cornsilk2")
    text(x=plot_dat$Var2,y=plot_dat$Var1,round(as.numeric(plot_dat$value),digits = 2),cex = 0.6)
  }else{
    if(all(Ks %in% exp(seq(log(0.1),log(10),length.out=100)))){
      K_labels <- c(seq(0.1,1,0.1),2:10)
      K_plot <- log10(Ks)
      K_ats <- log10(K_labels)
    }else{
      K_labels <- Ks
      K_plot <- Ks
      K_ats <- Ks
    }
    phis <- round(log10(K_labels) + 2 * log10(Linfs), digits = 2)

    par(mar = c(12,5,4,2))
    plot(K_plot,score_mat,type = 'l', ylim=c(0,max(score_mat, na.rm = TRUE)*1.4),
         ylab = "Score function", xlab = "Growth constant K (/year)", col = "red", lwd=2,xaxt='n')
    axis(1,at = K_ats,labels = K_labels)
    axis(1,at = K_ats,labels = phis,
         line = 5.5)
    mtext(text = expression(paste("Growth performance index (",phi,"')")),side = 1,line = 8.5)
    grid(nx = 0, NULL, lty = 6, col = "gray40")
    abline(v = K_ats, lty = 6, col = "gray40")
    title("K-Scan", line = 2)
  }

  if(!is.na(Linf_fix)){
    Rn_max <- max(score_mat, na.rm = TRUE)
    K_max <- as.numeric(dimnames(score_mat)[[1]][which(score_mat[,1] == Rn_max)])
    startingPoints <- ESP_starting_point_L[which(score_mat == Rn_max),,1]
  }


  ret <- c(res,list(restructuredFrequencies = catch_aAF_F,
                    score_mat = score_mat,
                    ESP_starting_point_L=ESP_starting_point_L,
                    C = C,
                    WP = WP))
  if(!is.na(Linf_fix)){
    ret$Rn_max <- Rn_max
    ret$Linf_fix <- Linf_fix
    ret$K_opt <- K_max
    ret$startingPoints <- startingPoints
  }
  class(ret) <- "lfq"
  return(ret)
}

