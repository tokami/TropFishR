#' @title Thompson and Bell prediction model
#
#' @description  This is a function to calculate the total mortality (Z) from length composition data via the length converted catch curve or from age at length data with catch curve.
#'
#' @param param A list
#' @param stock_size_1 Stock size to start with
#' @param plus.group Should a plus group be created, by default (NA) not. BUt authors advise to do so
#' @param FM_change Vector containing new FM values
#' @param Linf growth parameter 1
#' @param K growth parameter 2
#' @param t0 growth parameter 3
#' @param s_list List wit all selectivity information
#' @param Lc_change vector with ascending Lc values
#'
#' @param FM reference F-at-age-array
#'
#' @examples
#'
#' data(data_Predict_ThompsonBell)
#'
#' select.list <- list(selecType = 'knife_edge',  #or 'gillnet' or 'trawl_ogive'
#'    Lc = 34,tc = 5,selecDist = 'lognormal',    #or 'normal_fixed'
#'    mesh_size = 8.1,mesh_size1 = 9.1,select_p1 = 21.1,select_p2 = 23.8)
#'
#' Predict_ThompsonBell_FM_Lc(data_Predict_ThompsonBell,FM_change = seq(0,3,0.2),
#'   unit.time = 'month', Linf=50, K=0.3, t0=0.01, s_list=select.list,
#'   Lc_change=seq(24,44,2))
#'
#'
#' @details better to treat last group always as a plus group..... For variable parameter system vectors are reuqired for constant parameter systems matrices or data.frames have to be inserted. or vectors The length converted linearised catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @references example 1 : Kuwait (Garcia and van zalinge 1982)
#'   Millar, R. B., & Holst, R. (1997). Estimation of gillnet and hook selectivity using log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
#' @export

Predict_ThompsonBell_FM_Lc <- function(param,
                                       FM_change,
                                       stock_size_1 = NA,
                                       unit.time = 'year',
                                       plus.group = NA,
                                       Linf,
                                       K,
                                       t0 = 0,
                                       s_list,
                                       Lc_change){

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

  # instead of s_list the outcome of one of the other select functions?
  #as a option or put values per hand

  Lt <- Linf * (1- exp(-K * (classes.num - t0)))

  sel <- select_ogive(s_list,classes.num,Lt)

  sel.list <- list()
  for(x19 in 1:length(Lc_change)){
    sel.list[[x19]] <- select_ogive(s_list,classes.num, Lt, Lc_change[x19])
  }
  Lc_mat <- do.call(cbind,sel.list)
  colnames(Lc_mat) <- Lc_change

  Lc_mat_FM <- Lc_mat * max(FM,na.rm=TRUE)

  #list with FM_Lc_matrices per FM_change
  FM_Lc_com_mat.list <- list()
  for(x20 in 1:length(colnames(Lc_mat_FM))){
    FM_Lc_com_mat.list[[x20]] <- as.matrix(Lc_mat_FM[,x20]) %*% FM_change
    colnames(FM_Lc_com_mat.list[[x20]]) <- FM_change
  }

  param.loop <- res

  pred.FM_Lc_com_res_loopC_list <- list()
  pred.FM_Lc_com_res_loopY_list <- list()
  pred.FM_Lc_com_res_loopB_list <- list()
  pred.FM_Lc_com_res_loopV_list <- list()

  for(x21 in 1:length(FM_Lc_com_mat.list)){  #loop for length of list == Lc changes
    mati <- FM_Lc_com_mat.list[[x21]]

    pred.FM_Lc_com_res_loop1_list <- list()
    for(x22 in 1:dim(mati)[2]){

      param.loop$FM <- mati[,x22]
      param.loop$Z <- mati[,x22] + nM
      res2 <- Predict_ThompsonBell1(param.loop, unit.time,
                                    stock_size_1,plus.group=plus.group)
      pred.FM_Lc_com_res_loop1_list[[x22]] <- res2$totals
    }
    prev_mat <- do.call(rbind, pred.FM_Lc_com_res_loop1_list)
    prev_matC <- prev_mat[,'tot.C']
    prev_matY <- prev_mat[,'tot.Y']
    prev_matB <- prev_mat[,'meanB']
    prev_matV <- prev_mat[,'tot.V']


    pred.FM_Lc_com_res_loopC_list[[x21]] <- prev_matC
    pred.FM_Lc_com_res_loopY_list[[x21]] <- prev_matY
    pred.FM_Lc_com_res_loopB_list[[x21]] <- prev_matB
    pred.FM_Lc_com_res_loopV_list[[x21]] <- prev_matV
  }

  #for catch
  mat_FM_Lc_com.C <- do.call(rbind, pred.FM_Lc_com_res_loopC_list)
  rownames(mat_FM_Lc_com.C) <- Lc_change
  colnames(mat_FM_Lc_com.C) <- FM_change

  #for yield
  mat_FM_Lc_com.Y <- do.call(rbind, pred.FM_Lc_com_res_loopY_list)
  rownames(mat_FM_Lc_com.Y) <- Lc_change
  colnames(mat_FM_Lc_com.Y) <- FM_change

  #for biomass
  mat_FM_Lc_com.B <- do.call(rbind, pred.FM_Lc_com_res_loopB_list)
  rownames(mat_FM_Lc_com.B) <- Lc_change
  colnames(mat_FM_Lc_com.B) <- FM_change

  #for value
  mat_FM_Lc_com.V <- do.call(rbind, pred.FM_Lc_com_res_loopV_list)
  rownames(mat_FM_Lc_com.V) <- Lc_change
  colnames(mat_FM_Lc_com.V) <- FM_change

  # transvers matrices for plotting (the opposite arrangement from book)
  mat_FM_Lc_com.C <- t(mat_FM_Lc_com.C)
  mat_FM_Lc_com.Y <- t(mat_FM_Lc_com.Y)
  mat_FM_Lc_com.B <- t(mat_FM_Lc_com.B)
  mat_FM_Lc_com.V <- t(mat_FM_Lc_com.V)

  # colours for plot
  pal <- colorRampPalette(c(
    rgb(1,0.5,0.5), rgb(1,1,0.5), rgb(0.5,1,1), rgb(0.5,0.5,1)
  ))


  #plot
  image(x = FM_change,
        y = Lc_change,
        z = mat_FM_Lc_com.Y, col=pal(100),
        xlab = 'Fishing mortality', ylab = 'Lc')
  contour(x = FM_change,
          y = Lc_change,
          z = mat_FM_Lc_com.Y, add=TRUE)
  #mtext("Yield", line=0.5, side=3)

  ret <- c(res,
           list(Lt=Lt,sel=sel,
             mat_FM_Lc_com.C=mat_FM_Lc_com.C,
           mat_FM_Lc_com.Y=mat_FM_Lc_com.Y,
           mat_FM_Lc_com.V=mat_FM_Lc_com.V,
           mat_FM_Lc_com.B=mat_FM_Lc_com.B))
  return(ret)

} ## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?
