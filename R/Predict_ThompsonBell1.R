#' @title Thompson and Bell prediction model
#
#' @description  This is a function to calculate the total mortality (Z) from
#'    length composition data via the length converted catch curve or from age at
#'    length data with catch curve.
#'
#' @param param A list containing all information
#' @param unit.time Indicates if the age groups are per month (\code{"month"}) or per year (\code{"year"}). Default: \code{"year"}
#' @param stock_size_1 Stock size of smallest age/length group
#' @param plus.group Indicates age/length group, which should be turned into a plus group (i.e. all groups above are comprised in one group)
#'
#' @examples
#'
#' # age-based yield per recruit model
#' # load data
#' data("data_Predict_ThompsonBell")
#'
#' # run model
#' # option 1: without plus group
#' Predict_ThompsonBell1(data_Predict_ThompsonBell)
#'
#' # option 2: with plus group
#' Predict_ThompsonBell1(data_Predict_ThompsonBell,plus.group=11)
#'
#' # length-based yield per recruit model
#' # load data
#' data(hake)
#'
#' # run model
#' output <- Predict_ThompsonBell1(param = hake, stock_size_1 = 98919.3)
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

Predict_ThompsonBell1 <- function(param, unit.time = "year",
                                 stock_size_1 = NA, plus.group = NA){


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
  #ifelse(is.null(res$M),nM <- mean(Z - FM,na.rm=T),nM <- res$M)
  #ifelse(is.null(res$Z),Z <- FM + nM,Z <- res$Z)


  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                        Age data                          #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if('age' %in% names(res)){

    classes <- as.character(res$age)
    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    # delta t
    dt <- rep(NA,length(classes.num))
    for(x0 in 1:(length(classes.num)-1)){
      dt[x0] <- classes.num[x0+1] - classes.num[x0]
    }
    if(unit.time == 'month'){
      dt <- dt * 1/12
    }

    #population per age group
    N <- rep(NA,length(classes.num))
    if(!is.na(stock_size_1)){
      N[1] <- stock_size_1
    }else N[1] <- 1000  #if not provided results are just negative

    for(x1 in 2:length(classes.num)){
      N[x1] <- N[x1-1] * exp(-Z[x1] * dt[x1])
      if(x1 == length(N)){                                ####CORRECT?????
        N[x1] <- N[x1-1] * exp(-Z[x1-1] * dt[x1-1])
      }
    }

    #number of deaths per time step month or year
    dead <- NA
    for(x2 in 1:length(classes.num)){
      dead[x2] <- N[x2] - N[x2+1]
    }

    #catch in numbers
    C <- dead * FM / Z

    #yield in kg or g
    Y <- C * meanWeight

    #mean biomass in kg or g
    B <- Y / (FM * dt)

    #value expressed in money units
    V <- Y * meanValue

    #total catch, yield, value and average biomass
    tot.C <- sum(C, na.rm=TRUE)
    tot.Y <- sum(Y, na.rm=TRUE)
    tot.V <- sum(V, na.rm=TRUE)
    meanB <- sum((B * dt), na.rm=TRUE) / sum(dt, na.rm=TRUE)   ### more complicated biomass concept if dt is not constant, see Chapter 5
    totals <- data.frame(tot.C=tot.C,
                         tot.Y=tot.Y,
                         tot.V=tot.V,
                         meanB=meanB)

    res2 <- list(dt = dt, N = N, dead = dead, C = C,
                 Y = Y, B = B, V = V,
                 totals = totals)


    #with plus group
    if(!is.na(plus.group)){

      df <- do.call(cbind,res2[1:7])
      df <- df[1:plus.group,]
      classes <- classes[1:plus.group]

      #new class
      classes[length(classes)] <-
        paste(classes[length(classes)],"+",sep='')

      #num deaths
      df[plus.group, "dead"] <- df[plus.group, "N"]

      #catch
      new.C <- (FM[plus.group] / Z[plus.group]) * df[plus.group, "N"]
      catch.plus.dif <- new.C - df[plus.group, "C"]
      df[plus.group, "C"] <- new.C

      #yield
      df[plus.group, "Y"] <- meanWeight[plus.group] * catch.plus.dif

      #value
      df[plus.group, "V"] <-
        df[plus.group, "Y"] * meanValue[plus.group]

      #biomass       ####not sure....omitted in manual
      df[plus.group, "B"] <-
        df[plus.group, "Y"] / (FM[plus.group] * df[plus.group, "dt"])


      res2 <- as.list(as.data.frame(df))



      df2 <- do.call(cbind,res)
      df2 <- df2[1:plus.group,]
      res <- as.list(as.data.frame(df2))
      res$age <- classes

      #total catch, yield, value and average biomass
      tot.C <- sum(res2$C, na.rm=TRUE)
      tot.Y <- sum(res2$Y, na.rm=TRUE)
      tot.V <- sum(res2$V, na.rm=TRUE)
      meanB <- sum((res2$B * res2$dt), na.rm=TRUE) / sum(dt, na.rm=TRUE)   ### more complicated biomass concept if dt is not constant, see Chapter 5
      totals <- data.frame(tot.C=tot.C,
                           tot.Y=tot.Y,
                           tot.V=tot.V,
                           meanB=meanB)

      res2 <- c(res2,totals)
    }

    ret <- c(res,res2)
    return(ret)
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                       Length data                        #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if('midLengths' %in% names(res)){

    Linf <- res$Linf
    K <- res$K
    a <- res$a
    b <- res$b

    classes <- as.character(res$midLengths)
    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    #calculate size class interval
    int <- classes.num[2] - classes.num[1]

    # t of lower and upper length classes
    lowL <- classes.num - (int / 2)
    upL <- classes.num + (int/2)

    # H
    H <- ((Linf - lowL)/(Linf - upL)) ^ (nM/(2*K))

    #population
    N <- rep(NA,length(classes.num))
    if(is.na(stock_size_1)) stock_size_1 <- 1000
    N[1] <- stock_size_1
    for(x1 in 2:length(classes.num)){
      N[x1] <- N[x1-1] * ((1/H[x1-1]) - (FM[x1-1]/Z[x1-1])) /
        (H[x1-1] - (FM[x1-1]/Z[x1-1]))
    }

    #number of deaths per time step month or year
    dead <- NA
    for(x2 in 1:length(classes.num)){
      dead[x2] <- N[x2] - N[x2+1]
    }

    #catch
    C <- rep(NA,length(classes.num))
    for(x3 in 1:(length(classes.num)-1)){
      C[x3] <- (N[x3] - N[x3+1]) * (FM[x3]/Z[x3])
    }

    #average weight
    W <- a * ((lowL + upL)/2 ) ^ b

    #yield
    Y <- C * W

    #Value
    V <- Y * meanValue

    #biomass
    B <- rep(NA,length(classes.num))
    for(x4 in 1:(length(classes.num)-1)){
    B[x4] <- ((N[x4] - N[x4+1]) / Z[x4] ) * W[x4]
    }

    #last length group
    C[length(C)] <- (N[length(N)] - 0) * FM[length(FM)]/Z[length(Z)]
    W[length(W)] <- a * ((lowL[length(lowL)] + Linf)/2 ) ^ b
    Y[length(Y)] <- C[length(C)] * W[length(W)]
    B[length(B)] <- (N[length(N)] - 0) / Z[length(Z)] * W[length(W)]
    V[length(V)] <- Y[length(Y)] * meanValue[length(meanValue)]


    #total catch, yield, value and average biomass
    tot.C <- sum(C, na.rm=TRUE)
    tot.Y <- sum(Y, na.rm=TRUE)
    tot.V <- sum(V, na.rm=TRUE)
    meanB <- sum((B), na.rm=TRUE)
    totals <- data.frame(tot.C=tot.C,
                         tot.Y=tot.Y,
                         tot.V=tot.V,
                         meanB=meanB)

    res2 <- list(N = N, dead = dead, C = C,
                 Y = Y, B = B, V = V,
                 totals = totals)

    ret <- c(res,res2)
    return(ret)

  }
} ## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?

