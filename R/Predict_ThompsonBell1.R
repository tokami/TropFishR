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
with(data_Predict_ThompsonBell,Predict_ThompsonBell1(param,plus.group=11))
param = data_Predict_ThompsonBell
unit.time = "year"
stock_size_1 = NA
plus.group = NA

#test
plus.group=11


#param can be age or midLenghts



Predict_ThompsonBell1 <- function(param, select.param = NA, unit.time = "year",
                                 stock_size_1 = NA, plus.group = NA){


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

    res2 <- list(dt = dt, pop = pop, dead = dead, C = C,
                 Y = Y, B = B,value.Y = value.Y,
                 tot.C = tot.C, tot.Y = tot.Y, meanB = meanB,
                 tot.value = tot.value)


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


      res2 <- as.list(as.data.frame(df))



      df2 <- do.call(cbind,res)
      df2 <- df2[1:plus.group,]
      res <- as.list(as.data.frame(df2))
      res$age <- classes
    }
    # does not cut all vectors from beginning!!! careful if you return "res" in the end!

    return(c(res,res2))
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                       Length data                        #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if('midLengths' %in% names(res)){}
} ## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?

