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



Predict_ThompsonBell <- function(classes, mean_weight, FM, Z, value = NA, Tr,
                                 datatype, stock_size_1 = NA, plus.group = NA){

  df.TB <- cbind(classes,mean_weight,value,FM,Z)
  df.TB <- as.data.frame(df.TB)
  df.TB$classes <- as.character(df.TB$classes)

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(df.TB$classes, split="\\+"))
  df.TB$classes.num <- as.numeric(classes.num[,1])

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                        Age data                          #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if(datatype == 'age'){

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
    df.TB$num.caught <- NA
    for(x3 in 1:length(df.TB$num.caught)){
      df.TB$num.caught[x3] <- df.TB$num.deaths[x3] * (df.TB$FM[x3] / df.TB$Z[x3])
    }

    #yield in kg or g
    df.TB$yield <- NA
    for(x4 in 1:length(df.TB$yield)){
      df.TB$yield[x4] <- df.TB$num.caught[x4] * df.TB$mean_weight[x4]
    }

    #mean biomass in kg or g
    df.TB$biomass <- NA
    for(x5 in 1:length(df.TB$biomass)){
      df.TB$biomass[x5] <- df.TB$yield[x5] / (df.TB$FM[x5] * df.TB$dt[x5])
    }

    #value expressed in money units
    df.TB$value_yield <- NA
    for(x6 in 1:length(df.TB$value_yield)){
      df.TB$value_yield[x6] <- df.TB$yield[x6] * df.TB$value[x6]
    }

    #total yield, value and average biomass
    total.yield <- sum(df.TB$yield,na.rm=TRUE)
    total.value <- sum(df.TB$value_yield,na.rm=TRUE)
    mean.biomass <- sum((df.TB$biomass * df.TB$dt),na.rm=TRUE) /
      sum(df.TB$dt,na.rm=TRUE)   ### more complicated biomass concept if dt not constant, see Chapter 5


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


  }



  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                       Length data                        #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if(datatype == 'length'){
    }
} ## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?

