#' @title Prediction models: YPR and Thompson & Bell model
#
#' @description  This is a function predicting catch, yield, biomass and economic values for different
#'    fishing mortality scenarions (in combination with gear changes). It allows to choose if the
#'    calculation should be based upon Beverton and Holt's yield per recruit model or the
#'    more sophisticated Thompson and Bell model.
#'
#' @param param a list consisting of following parameters (not all are required):
#' \itemize{
#'   \item \strong{Linf} or \strong{Winf}: infinite length for investigated species in cm [cm],
#'   \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'   \item \strong{t0}: theoretical time zero, at which individuals of this species hatch,
#'   \item \strong{M}: natural mortality,
#'   \item \strong{Z}: total mortality,
#'   \item \strong{FM}: fishing mortality,
#'   \item \strong{a}: ,
#'   \item \strong{b}: ,
#'   \item \strong{Lr} or \strong{tr}: length or age of recuritment;}
#' additional list objects in for the Thompson and Bell model:
#'  \itemize{
#'   \item \strong{age} or \strong{midLengths}: midpoints of the length class as vector (length-frequency
#'   data) or ages as vector (age composition data),
#'   \item \strong{meanWeight}: vector with mean weight per length group or age class,
#'   \item \strong{meanValue}: vector with mean value per length group or age class,
#' }
#' @param FM_change vector with ascending fishing mortalities
#' @param Lc_tc_change vector with ascending lengths or ages at first capture (Lc/tc)
#' @param type indicating which model should be applied: \code{"ypr"} for Beverton and Holt's
#'   yield per recruit model and \code{"ThompBell"} for the Thompson and Bell model
#' @param s_list
#' @param stock_size_1
#' @param unit.time
#' @param plus.group
#'
#' @keywords function prediction
#'
#' @examples
#' #______________________________________
#' #       Yiel Per Recruit (YPR)
#' #______________________________________
#' # age structured data
#' # Nemipterus marginatus
#' threadfin <- list(Winf = 286,K = 0.37, t0 = -0.2, M = 1.1, tr = 0.4)
#'
#' # run model
#' predict_mod(threadfin, FM_change = seq(0,6,0.1),
#'    Lc_tc_change = seq(0.2,1,0.2), type = 'ypr')  #where it is maximal  = MSY
#'
#' # Leiognathus spendens (Pauly 1980)
#' ponyfish <- list(Winf = 64, K = 1, t0 = -0.2, M = 1.8, tr = 0.2)
#'
#' # run model
#' predict_mod(ponyfish, Lc_tc_change = c(0.2,0.3,1.0), type = 'ypr')
#'
#' #______________________________________
#' # length structured data
#' # Xiphias gladius (Berkeley and Houde 1980)
#' swordfish <- list(Linf = 309, K = 0.0949, M = 0.18, a=0.0003, b=3, Lr = 90)  ## T_Lr , a, b ??? assumed
#'
#' # run model
#' predict_mod(swordfish, Lc_tc_change = c(100,118,150,180), type = 'ypr')
#'
#' ####test: E <- seq(0,0.9,0.1) FM <- E * M / (1 - E)
#'
#'
#' #______________________________________
#' #      Thompson and Bell model
#' #______________________________________
#' # with age structured data
#' # load data
#' data(shrimps)
#'
#' # run model
#' predict_mod(shrimps, FM_change = seq(0.1,3,0.1), type = 'ThompBell')
#'
#' # create list with selectivity information
#' select.list <- list(selecType = 'knife_edge',  #or 'gillnet' or 'trawl_ogive'
#'    Lc = 34, tc = 5, selecDist = 'lognormal',    #or 'normal_fixed'
#'    mesh_size = 8.1, mesh_size1 = 9.1, select_p1 = 21.1, select_p2 = 23.8)
#'
#' # add additional parameters to data list
#' shrimps <- c(shrimps,  list(Linf = 50, K = 0.3, t0 = 0.01))
#'
#' # run model
#' predict_mod(shrimps,FM_change = seq(0,3,0.2), Lc_tc_change = seq(24,44,2),
#'    type = 'ThompBell', s_list = select.list, unit.time = 'month')
#'
#' #______________________________________
#' # with length structured data
#' # load data
#' data(hake)
#'
#' # run model
#' predict_mod(hake,FM_change = seq(0.1,3,0.1), type = 'ThompBell')
#'
#'
#' @details better to treat last group always as a plus group.....
#'    The Thompson and Bell model incorporates an iteration step simulating the stock by means
#'    of the \code{\link{stock_sim}} function. In case changes in gear characteristics -
#'    here measured in terms of Lc or tc, the length or age at first capture, respectively -
#'    should be explored, a list with selectivity information about the gear has to be provided and
#'    the prediction models make use of the selectivity \code{\link{select_ogive}} function.
#'
#' @return A list with the input parameters and dependent on the model type following
#'    list objects:
#' \itemize{
#'   \item \code{type = 'ypr'}
#'   \itemize{
#'      \item \strong{FM}: fishing mortalities,
#'      \item \strong{Lc} or \strong{tc}: lengths or ages at first capture,
#'      \item \strong{list_Lc_runs} or \strong{list_tc_runs}: a list with the dataframes for each Lc or tc value:
#'      \itemize{
#'      \item FM_change
#'      \item Y_R
#'      \item Y_R.rel
#'      \item B_R
#'      \item B_R.percent
#'      \item Ty
#'      \item Wy
#'      }
#'   }
#'   \item \code{type = 'ThomBell'}
#'   \itemize{
#'      \item \strong{dt}: delta t,
#'      \item \strong{N}: population number,
#'      \item \strong{dead}: deaths due to natural reasons,
#'      \item \strong{C}: catch,
#'      \item \strong{Y}: yield,
#'      \item \strong{B}: biomass,
#'      \item \strong{V}: value,
#'      \item \strong{totals}: summed up values (total catch, total yield, total value, average biomass),
#'      \item \strong{tot.C}: total catches for different x factors,
#'      \item \strong{tot.Y}: total yield values for different x factors,
#'      \item \strong{tot.V}: total values for different x factors,
#'      \item \strong{meanB}: average biomasses for different x factors,
#'      \item \strong{Xfact}: fishing mortality changes;
#'   }
#'   \item \code{type = 'ThomBell'} and \code{Lc_tc_change} provided
#'   \itemize{
#'      \item \strong{FM_change}: fishing mortality changes,
#'      \item \strong{Lc_tc_change}: changes in length or age at first capture,
#'      \item \strong{Lt}: lengths at age,
#'      \item \strong{sel}: probability of capture,
#'      \item \strong{mat_FM_Lc_com.C}: catch matrix for all fishing mortality and Lc/tc combinations,
#'      \item \strong{mat_FM_Lc_com.Y}: yield matrix for all fishing mortality and Lc/tc combinations,
#'      \item \strong{mat_FM_Lc_com.V}: value matrix for all fishing mortality and Lc/tc combinations,
#'      \item \strong{mat_FM_Lc_com.B}: biomass matrix for all fishing mortality and Lc/tc combinations;
#'   }
#'
#' }
#'
#' @references
#' example 1 : Kuwait (Garcia and van Zalinge 1982)
#'
#' Beverton and Holt 1966
#'
#' Boerema, L.K., and J.A. Gulland, 1973. Stock assessment of the Peruvian anchovy
#' (Engraulis ringens) and management of the fishery. Journal of the Fisheries Board of
#' Canada, 30(12):2226-2235
#'
#' Millar, R.B., and R. Holst, 1997. Estimation of gillnet and hook selectivity using
#' log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3):471-477
#'
#' Sparre, P., and Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

predict_mod <- function(param, FM_change = NA, Lc_tc_change = NULL, type,  s_list,
                        stock_size_1 = NA, unit.time = 'year', curr.E = NA,
                        plus.group = NA){
  res <- param

  # Beverton and Holt's ypr
#------------
  if(type == "ypr"){
    M <- res$M
    K <- res$K
    t0 <- ifelse(!is.null(res$t0),res$t0,0)
    a <- ifelse(!is.null(res$a),res$a,NA)
    b <- ifelse(!is.null(res$b),res$b,NA)

    if(length(FM_change) == 1 & is.na(FM_change[1])){
      FM_change <- seq(0,10,0.1)
      print(noquote("No fishing mortality (FM_change) was provided, a default range of 0 to 10 is used."))
    }

    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #                        Age data                          #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if("Winf" %in% names(res)){
      Winf <- res$Winf
      tr <- res$tr
      tc <- Lc_tc_change
      list_tc_runs <- list()
      for(i in 1:length(tc)){
        tci <- tc[i]

        #Biomass per Recruit
        S <- exp(-K * (tci - t0))
        B_R <- exp(-M*(tci-tr)) * Winf *
          ((1/(FM_change+M)) - ((3*S)/((M+FM_change)+K)) +
             ((3*(S^2))/((M+FM_change)+(2*K))) - ((S^3)/((M+FM_change) + (3*K))))

        #virgin biomass
        Bv_R <- B_R[which(FM_change == 0)]
        #biomass of exploited part of the cohort (biomass of fish older than tc)

        #biomass in percetage of virgin biomass
        B_R.percent <- round((B_R / Bv_R ) * 100, digits = 1)

        #Yield per Recruit
        Y_R <- B_R * FM_change

        #relative yield per recruit - mostly done with length frequency data (exclusively?)
        Y_R.rel <- Y_R * (exp(-M *(tr - t0))) / Winf


        #mean age in annual yield
        Ty <- (1 / (M+FM_change)) + tci

        #mean length in the annual yield
        #Ly <- Linf * (1 - (((M+FM_change)*S)/((M+FM_change)+K)))

        #mean weight in annual yield
        Wy <- (M+FM_change) * Winf *
          ((1/(FM_change+M)) - ((3*S)/((M+FM_change)+K)) +
             ((3*(S^2))/((M+FM_change)+(2*K))) - ((S^3)/((M+FM_change) + (3*K))))


        results.PBH <- data.frame(FM_change = FM_change,
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
                  FM = FM_change,
                  tc = tc,
                  list_tc_runs = list_tc_runs)
      class(ret) <- "predict_mod"
      # plot results
      plot(ret)

      return(ret)
    }

    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #                       Length data                        #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if("Linf" %in% names(res)){
      Linf <- res$Linf
      Lr <- res$Lr
      Lc <- Lc_tc_change

      # yield functions
      # ___________________________________________________
      S <- 1 - (Lci/Linf)  # == U  ##(U <- 1 - (Lci/Linf))

      ypr <- function(x){
        A <- ((Linf - Lci)/(Linf - Lr)) ^ (M/K)
        G <- x + M   # G == Z
        y <- x * A * Winf * ((1/G) - (3*S)/(G + K) + (3*S^2)/(G + 2*K) - (S^3)/(G + 3*K))
        return(y)
      }

      #according to Sparre & Venema and Gayanilo 2006:
      ypr.rel <- function(x){
        m <- (1-x)/(M/K)    ## == K/Z
        y <- x * S^(M/K) * (1 - ((3*S)/(1+m)) + ((3*S^2)/(1+2*m)) - ((S^3)/(1+3*m)))
        return(y)
      }

      derivative <- function(x){
        C <- ((K*(1-x))/M)
        B <- (S^(M/K)) * (1 - ((3*S)/(1+C)) + ((3*S^2)/(1+2*C)) - ((S^3)/(1+3*C)))
        D <- (-((3*K*S^3)/(M*((3*K*(1-x)/M)+1)^2)) + ((6*K*S^2)/(M*((2*K*(1-x)/M)+1)^2)) -
                ((3*K*S)/(M*((K*(1-x)/M)+1)^2)))
        des <- x * (S^(M/K)) * D  + B
        return(des)
      }
      # ___________________________________________________

      list_Lc_runs <- list()
      list_Es <- list()
      for(i in 1:length(Lc)){
        Lci <- Lc[i]

        Z <- (M + FM_change)
        E <- FM_change/Z
        # transform Linf in Winf
        Winf <- a * (Linf ^ b)

        #yiel per recurit for length data   # also possbile inputing option: F/K
        Y_R <- ypr(FM_change)

        #biomass per recruit for length data?
        B_R <- Y_R / FM_change

        #virgin biomass
        Bv_R <- B_R[(which(FM_change == 0)+1)]  ### CHECK: if == 0 than 0 NaN because yield and FM_change is 0
        #biomass of exploited part of the cohort (biomass of fish older than tc)

        #biomass in percetage of virgin biomass
        B_R.percent <- round((B_R / Bv_R ) * 100, digits = 1)

        #relative yield per recruit - mostly done with length frequency data (exclusively?)
        Y_R.rel <- ypr.rel(E)
        #according to Gayanilo 1997 Fisat description (wrong???):
#         Y_R.rel <- E * S^m * (1 - ((3*S)/(1+m)) +
#                                     ((3*S^2)/(1+2*m)) - ((S^3)/(1+3*m)))

        #mean length in the annual yield
        Ly <- Linf * (1 - ((Z*S)/(Z+K)))

        #mean weight in annual yield
        # Wy <- (M+FM_change) * Winf *
        #    ((1/(FM_change+M)) - ((3*S)/((M+FM_change)+K)) +
        #       ((3*(S^2))/((M+FM_change)+(2*K))) - ((S^3)/((M+FM_change) + (3*K))))

        results.PBH <- data.frame(FM = FM_change,
                                  Ly = Ly,
                                  E = E,
                                  Y_R = Y_R,
                                  Y_R.rel = Y_R.rel,
                                  B_R = B_R,
                                  B_R.percent = B_R.percent)
        list_Lc_runs[[i]] <- results.PBH


        # First derivative of relative yield per recruit model
        deri <- derivative(E)

        # reference points
        N01 <- which.min(abs(deri - (deri[1] * 0.1)))
        N05 <- which.min(abs(B_R.percent - 50))    ##which.min(abs(deri - (deri[1] * 0.5)))
        Nmax <- which.min(abs(deri))

        df_loop_Es <- data.frame(Lc = Lci,
          F01 = FM_change[N01],
          F05 = FM_change[N05],
          Fmax = FM_change[Nmax],
          E01 = E[N01],
          E05 = E[N05],
          Emax = E[Nmax]
          )
        list_Es[[i]] <- df_loop_Es
      }

      df_Es <- do.call(rbind,list_Es)

      # current exploitation rate
      curr.F = (M * curr.E)/(1-curr.E)
      df_currents <- data.frame(curr.E = curr.E,
                                curr.F = curr.F,
                                curr.YPR = ypr((curr.F)),
                                curr.YPR.rel = ypr.rel(curr.E))

      names(list_Lc_runs) <- Lc
      ret <- list(res,
                  FM = FM_change,
                  Lc = Lc,
                  list_Lc_runs = list_Lc_runs,
                  df_Es = df_Es,
                  currents = df_currents)
      class(ret) <- "predict_mod"

      # plot results
      plot(ret)

      return(ret)
    }
  }

#------------


  # Thompson and Bell model
#------------
  if(type == "ThompBell"){
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

    Linf <- res$Linf
    K <- res$K
    t0 <- ifelse("t0" %in% names(res),res$t0,0)

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

    #HHHHHHHHHHHH#
    #     ALL    #
    #HHHHHHHHHHHH#

    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    # Only FM change provided without Lc_tc change
    if(is.null(Lc_tc_change)){
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

      res2 <- pred_res_df
      ret <- c(res,res2)
      class(ret) <- "predict_mod"

      # plot results
      plot(ret)

      return(ret)
    }


    # FM and Lc_tc change provided
    if(!is.null(Lc_tc_change)){
      # instead of s_list the outcome of one of the other select functions?
      #as a option or put values per hand

      Lt <- Linf * (1- exp(-K * (classes.num - t0)))

      sel <- select_ogive(s_list,Lt = Lt) #classes.num

      sel.list <- list()
      for(x19 in 1:length(Lc_tc_change)){
        sel.list[[x19]] <- select_ogive(s_list, Lt = Lt, Lc = Lc_tc_change[x19]) #classes.num
      }
      Lc_mat <- do.call(cbind,sel.list)
      colnames(Lc_mat) <- Lc_tc_change

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
          res2 <- stock_sim(param.loop, unit.time,
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
      rownames(mat_FM_Lc_com.C) <- Lc_tc_change
      colnames(mat_FM_Lc_com.C) <- FM_change

      #for yield
      mat_FM_Lc_com.Y <- do.call(rbind, pred.FM_Lc_com_res_loopY_list)
      rownames(mat_FM_Lc_com.Y) <- Lc_tc_change
      colnames(mat_FM_Lc_com.Y) <- FM_change

      #for biomass
      mat_FM_Lc_com.B <- do.call(rbind, pred.FM_Lc_com_res_loopB_list)
      rownames(mat_FM_Lc_com.B) <- Lc_tc_change
      colnames(mat_FM_Lc_com.B) <- FM_change

      #for value
      mat_FM_Lc_com.V <- do.call(rbind, pred.FM_Lc_com_res_loopV_list)
      rownames(mat_FM_Lc_com.V) <- Lc_tc_change
      colnames(mat_FM_Lc_com.V) <- FM_change

      # transvers matrices for plotting (the opposite arrangement from book)
      mat_FM_Lc_com.C <- t(mat_FM_Lc_com.C)
      mat_FM_Lc_com.Y <- t(mat_FM_Lc_com.Y)
      mat_FM_Lc_com.B <- t(mat_FM_Lc_com.B)
      mat_FM_Lc_com.V <- t(mat_FM_Lc_com.V)

      ret <- c(res,
               list(FM_change = FM_change,
                    Lc_tc_change = Lc_tc_change,
                 Lt=Lt,sel=sel,
                    mat_FM_Lc_com.C=mat_FM_Lc_com.C,
                    mat_FM_Lc_com.Y=mat_FM_Lc_com.Y,
                    mat_FM_Lc_com.V=mat_FM_Lc_com.V,
                    mat_FM_Lc_com.B=mat_FM_Lc_com.B))
      class(ret) <- "predict_mod"

      # plot results
      plot(ret)

      return(ret)
    }
  }
#------------

}





## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?

