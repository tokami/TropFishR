#' @title Date - Year conversion
#'
#' @description Convert dates to numeric years with decimal as fraction of a year
#'
#' @param date a date (class 'Date')
#'
#' @examples
#' date2yeardec(Sys.Date())
#'
#' @return a scalar (class 'numeric')
#' @export

date2yeardec <- function(date){
  # adapted from lubridate::decimal_date.default
  if (any(!inherits(date, c("POSIXt", "POSIXct", "POSIXlt", "Date")))){
    stop("date(s) not in POSIXt or Date format")
  }
  Y <- as.numeric(format(date, format="%Y"))
  start <- as.POSIXct(paste0(Y,  "/01/01"), tz="UTC")
  end   <- as.POSIXct(paste0(Y+1,"/01/01"), tz="UTC")
  sofar <- as.numeric(difftime(date, start, units = "secs"))
  total <- as.numeric(difftime(end, start, units = "secs"))
  res <- Y + sofar/total
  return(res)
}

#' @title Convert FiSAT's starting point to ta value
#'
#' @description Starting points returned or chosen within FiSAT are not supported
#'    in TropFishR. Instead \code{ta} takes on the job of anchoring VBGF growth curves on
#'    a temporal axis. This function allows to convert FiSAT's starting points to \code{ta} values
#'
#' @param lfq list with dates, midLengths, and catch
#' @param par list with growth parameters 'Linf' and 'K' of VBGF
#' @param startingLength starting length as returned by FiSAT, indicating the length within
#'    the starting sample cut by a growth curve
#' @param startingSample starting sample as returned by FiSAT, indicating the sample which is
#'    cut by a growth curve
#'
#' @keywords function lfq startingPoints ta
#'
#' @examples
#' data(synLFQ5)
#' lfqNEW <- startingPoint2tanchor(synLFQ5, par = list(Linf = 92, K = 0.37),
#'    startingLength = 31, startingSample = 4)
#' lfqRest <- lfqRestructure(lfqNEW, MA = 11)
#' plot(lfqRest,par=list(Linf=lfqRest$Linf,K=lfqRest$K,ta=lfqRest$ta))
#'
#' @return list with input elements and estimated ta value
#'
#' @export

startingPoint2tanchor <- function(lfq, par, startingLength, startingSample){

  res <- lfq
  if(("Linf" %in% names(par)) == FALSE) stop(noquote("Please add Linf to param!"))
  if(("K" %in% names(par)) == FALSE) stop(noquote("Please add K to param!"))
  if(("dates" %in% names(res)) == FALSE) stop(noquote("Please add a dates to param!"))

  if(length(res$dates) < startingSample) stop(noquote(paste("Only",
                                                            length(res$dates),
                                                            "sample dates in param$dates, but startingSample is",
                                                            startingSample)))

  tx <- res$dates[startingSample]
  tx <- date2yeardec(tx)
  agex <- VBGF(L = startingLength, pars = list(Linf = par$Linf, K = par$K, t0 = 0))
  ta <- (tx-agex) %% floor(tx-agex)

  ret <- c(res, list(Linf = par$Linf, K = par$K, ta = ta))
  return(ret)
}

#' @title Year - Date conversion
#'
#' @description Convert numeric years to dates
#'
#' @param yeardec numeric year
#'
#' @examples
#'  yeardec2date(2014.14)
#'
#' @return date in format "\%Y-\%m-\%d" (class 'Date').
#'
#' @export

yeardec2date <- function(yeardec){
  # adapted from lubridate::date_decimal
  start <- as.POSIXct(paste0(trunc(yeardec),  "/01/01"), tz="UTC")
  end   <- as.POSIXct(paste0(trunc(yeardec)+1,"/01/01"), tz="UTC")
  res <- as.Date(start + (difftime(end, start, units="secs") * (yeardec - trunc(yeardec))))
  return(res)
}

#' @title Create lfq data from length measurements
#'
#' @description Convert raw length measurements to length frequency data (lfq class).
#'
#' @param data data with at least two columns, one with the length measurements, one
#'    with the sampling date
#' @param Lname name of the length column
#' @param Dname name of the date column
#' @param Fname optional; name of column with frequency, in case each length was measured more than
#'              one time
#' @param bin_size size of the bins in cm (Default: 2)
#' @param species character; to store species name in lfq list
#' @param stock character; to store stock ID or name in lfq list
#' @param comment optional character; to store comments conerning the lfq list
#' @param Lmin minimum length for the midLengths vector (default: 0)
#' @param length_unit unit of length measurements, either "cm" (default), "mm" or "m"
#' @param plus_group logical; should a plus group be created? If yes you will be
#'    asked to insert the length for the plus group in the console (default: FALSE).
#'    Instead of inserting the length of the plus group via the console, the value
#'    can be incorporated in a vector, e.g. plus_group = c(TRUE, 30).
#' @param aggregate_dates logical; indicating whether dates should be lumped in monthly
#'    sampling times (assuming sampling always aound the 15th of
#'    each month; default is FALSE). More exact lumping can only done manually and
#'    then sampling dates provided in data.
#' @param plot logical; should a graph of lfq data be displayed? (Default: FALSE)
#'
#' @keywords function lfq length-frequency
#'
#' @examples
#' # create random data
#' set.seed(1)
#' data <- data.frame(length.mm. = sample(c(rpois(300, lambda = 60),
#'            rpois(200, lambda = 100), rpois(100, lambda = 150)),
#'            size = 1000, replace = TRUE),
#'            dates = seq.Date(as.Date("2015-10-02"),as.Date("2016-08-28"),
#'            length.out = 1000))
#' # create lfq data
#' lfq_dat <- lfqCreate(data,Lname = "length.mm.", Dname = "dates", aggregate_dates = TRUE,
#'    length_unit = "mm", bin_size = 0.5, plot=TRUE, plus_group=c(TRUE,15.75))
#'
#' @return A list of "lfq" class with
#'    \itemize{
#'    \item \strong{dates} dates of sampling times (class Date),
#'    \item \strong{midLengths} midpoints of the length classes,
#'    \item \strong{catch} matrix with catches/counts per length class (row) and
#'   sampling date (column).
#'    }
#'
#' @export

lfqCreate <- function(data, Lname, Dname, Fname = NA, bin_size = 1,
                      species = NA, stock = NA, comment = "",
                      Lmin = 0,
                      length_unit = "cm", plus_group = FALSE,
                      aggregate_dates = FALSE,
                      plot = FALSE){

    data$length <- get(Lname, data)
    data$date <- get(Dname, data)
    if(!is.na(Fname)){
        data$freq <- get(Fname, data)
    }else{
        data$freq <- rep(1,nrow(data))
    }
    if(class(data$date) != "Date") stop(noquote("Please provide the date as 'Date' class (e.g. as.Date())."))

    # convert length if necessary
    if(length_unit == "m") data$length <- data$length * 100
    if(length_unit == "mm") data$length <- data$length / 10

    # delete any length of 0
    data$length[which(data$length == 0)] <- NA

    # show histogram
    ## hist(data2$length, breaks = 50)

    # no NAs allowed in length nor date column
    data <- data[!is.na(data$length),]
    data <- data[!is.na(data$date),]

    # order according to date
    data <- data[order(data$date),]

    # create monthly grouped dates
    if(aggregate_dates){
      data$samplings <- as.Date(paste(format(data$date, "%Y-%m"),"15",sep="-"))
    }else data$samplings <- data$date

    # rearrange data into LFQ data
    bin.breaks <- seq(Lmin, max(data$length) + bin_size, by=bin_size)
    midLengths <- bin.breaks + bin_size/2

    data2 <- aggregate(list(freq=data$freq),
                       by=list(date=data$samplings, length=data$length), sum)

    # order according to date
    data2 <- data2[order(data2$date),]


    listi <- vector("list",length(unique(data2$date)))
    LF_dat <- data.frame(bin = bin.breaks)
    for(i in 1:length(unique(data2$date))){

        sampli <- unique(data2$date)[i]
        lengthi <- as.numeric(data2$length[data2$date == sampli])
        freqi <- as.numeric(data2$freq[data2$date == sampli])

        bin.breaks2 <- rep(NA, length(bin.breaks))
        for(ii in 1:length(bin.breaks)){
            if(ii == length(bin.breaks)){
                bin.breaks2[ii] <- length(which(lengthi >= bin.breaks[ii]))
            }else{
                bin.breaks2[ii] <- length(which(lengthi >= bin.breaks[ii] & lengthi < bin.breaks[ii+1]))
            }
        }

        bin.breaks3 <- rep(bin.breaks, bin.breaks2)
        dati <- aggregate(list(freq=freqi), by=list(bin=bin.breaks3), sum)

        listi[[i]] <- merge(LF_dat, dati, by.x = "bin", all.x =TRUE)[,2]
    }
    catch_mat <- do.call(cbind,listi)
    catch_mat[is.na(catch_mat)] <- 0

    # plus group
    if(plus_group[1]){
      if(length(plus_group) == 1){
        if(is.vector(catch_mat)){
          print(data.frame(midLengths = midLengths, frequency = catch_mat))
        }else print(data.frame(midLengths = midLengths, frequency = rowSums(catch_mat)))
        writeLines("Check the table above and insert the length of the plus group (Esc to cancel).")
        pg = -1
        while(pg > max(midLengths) | pg < min(midLengths)){
          pg <- readline(paste0("Enter a length group between ", min(midLengths)," and ",
                                max(midLengths),":"))
          pg = as.numeric(as.character(pg))
          if(!(pg %in% midLengths)){
            writeLines(paste0(pg, " is not an element of midLengths (see table)."))
            pg = -1
            #pg <- ifelse(grepl("\\D",pg),-1,as.integer(pg))
            if(is.na(pg)){break}  # breaks when hit enter
          }
        }
      }else if(length(plus_group) == 2){
        pg = as.numeric(as.character(plus_group[2]))
      }

      midLengths <- midLengths[1:which(midLengths == pg)]
      if(is.vector(catch_mat)){
        addplus <- sum(catch_mat[(which(midLengths == pg):length(catch_mat))])
        catch_mat <- catch_mat[1:which(midLengths == pg)]
        catch_mat[which(midLengths == pg)] <-
          catch_mat[which(midLengths == pg)] + addplus
      }else{
        addplus <- colSums(catch_mat[(which(midLengths == pg):nrow(catch_mat)),])
        catch_mat <- catch_mat[1:which(midLengths == pg),]
        catch_mat[which(midLengths == pg),] <-
          catch_mat[which(midLengths == pg),] + addplus
      }
    }



    res <- list(species = species,
                stock = stock,
                dates = unique(data$samplings),
                midLengths = midLengths,
                catch = catch_mat,
                comment = comment)
    class(res) <- "lfq"
    if(plot) plot(res, Fname = "catch")
    return(res)
}

#' @title Modify lfq data for further analysis
#'
#' @description Modify length-freqeuncy (LFQ) data. Allows to summarise catch matrix
#'    of LFQ data to one column per year. This is required for e.g. \code{\link{catchCurve}}.
#'    Allows to change bin size of LFQ data. Allows to ad plus group to catch matrix.
#'
#' @param lfq lfq object with dates, midLengths, and catch
#' @param par growth parameters as resulting from e.g. \code{\link{ELEFAN}}
#' @param bin_size Bin size for length frequencies (in cm)
#' @param aggregate Factor to aggregate catch per year (\code{"year"}),
#'    per quarter (\code{"quarter"}), or per month (\code{"month"}). By default data
#'    is not aggregated (\code{NA}).
#' @param vectorise_catch logical; indicating if the catch matrix should be summarised to
#'    yearly vectors (default: FALSE).
#' @param plus_group logical or numeric; should a plus group be created? If yes you will be
#'    asked to insert the length for the plus group in the console (default: FALSE).
#'    Instead of inserting the length of the plus group via the console, the value
#'    can be inserted, e.g. plus_group = 85.5.
#' @param minDate minimum date to subset lfq data
#' @param maxDate maximum date to subset lfq data
#' @param years numeric with year(s) to subset lfq data
#' @param Lmin minimum length to subset lfq data
#' @param Lmax maximum length to subset lfq data
#' @param lfq2 optional second lfq object which will be merged with lfq. This might be interesting for
#'    fleet specific lfq objects. Default: NA. Be aware that catches are combined without weighting!
#'
#' @keywords function lfq length-frequency
#'
#' @examples
#' data(synLFQ4)
#'
#' ## summarise catch matrix per year
#' lfq_sum <- lfqModify(synLFQ4, vectorise_catch = TRUE)
#'
#' ## change bin size
#' lfq_bin <- lfqModify(synLFQ4, bin_size = 4)
#'
#' ## add plus_group
#' lfq_plus <- lfqModify(synLFQ4, plus_group = 85.5)
#'
#' @return lfq object with rearranged catch matrix (yearly sums) and growth parameters
#'    if provided.
#'
#' @export

lfqModify <- function(lfq, par = NULL,
                      bin_size = NA,
                      aggregate = NA,  ## either: year, quarter, month  (year will substitute vectorise_catch)
                      vectorise_catch = FALSE,
                      plus_group = FALSE,
                      minDate = NA,
                      maxDate = NA,
                      years = NA,
                      Lmin = NA,
                      Lmax = NA,
                      lfq2 = NA){

    if(class(lfq) != "lfq") stop("Your lfq data set has to have class 'lfq'!")
    dates <- lfq$dates
    midLengths <- lfq$midLengths
    catch <- lfq$catch

    ## linf for plus_group definition
    if(!is.null(par)){
        linf <- par$Linf
    }else{
        if("Linf" %in% names(lfq)){
            linf <- lfq$Linf
        }else if("par" %in% names(lfq)){
            linf <- lfq$par$Linf
        }else{
            linf <- NA
        }
    }

    ## replace NAs in catch
    catch[which(is.na(catch))] <- 0

    ## select beyond certain date
    if(!is.na(minDate)){
        catch <- lfq$catch[,which(dates >= minDate)]
        dates <- lfq$dates[which(dates >= minDate)]
    }

    ## select before certain date
    if(!is.na(maxDate)){
        catch <- catch[,which(dates <= maxDate)]
        dates <- dates[which(dates <= maxDate)]
    }

    ## select certain years
    if(!is.na(years[1])){
        catch <- catch[,which(format(dates,"%Y") %in% years)]
        dates <- dates[which(format(dates,"%Y") %in% years)]
    }

    ## select above certain length
    if(!is.na(Lmin)){
        catch <- catch[which(midLengths >= Lmin),]
        midLengths <- midLengths[which(midLengths >= Lmin)]
    }

    ## select below certain length
    if(!is.na(Lmax)){
        catch <- catch[which(midLengths <= Lmax),]
        midLengths <- midLengths[which(midLengths <= Lmax)]
    }


    ## merge two lfq data sets (ADD weighing factor)
    if(!any(is.na(lfq2))){
        if(class(lfq2) != "lfq") stop("Your lfq2 data set has to have class 'lfq'!")

        ## extract variables
        dates2 <- lfq2$dates
        midLengths2 <- lfq2$midLengths
        catch2 <- lfq2$catch

        ## error messages
        if(diff(midLengths)[1] != diff(midLengths2)[1]) stop("The bin sizes do not fit eachother")
        if(any(!dates2 %in% dates)) warning("At least one sampling date of lfq2 does not match with the dates \nin lfq, not matching dates will be added!")

        mergi <- merge(data.frame(dates=dates,x=dates),
                       data.frame(dates=dates2,y=dates2),
                       by="dates",all=TRUE)
        mergi2 <- merge(data.frame(midLengths=midLengths,x=midLengths),
                        data.frame(midLengths=midLengths2,y=midLengths2),
                        by="midLengths",all=TRUE)
        indY <- which(is.na(mergi2$y) & mergi2$midLengths > max(midLengths2))
        matY <- matrix(0, nrow=length(indY), ncol=ncol(catch2))
        catch2 <- rbind(catch2,matY)
        indY <- which(is.na(mergi2$y) & mergi2$midLengths < min(midLengths2))
        matY <- matrix(0, nrow=length(indY), ncol=ncol(catch2))
        catch2 <- rbind(matY,catch2)
        ind <- which(is.na(mergi2$x) & mergi2$midLengths > max(midLengths))
        mat <- matrix(0, nrow=length(ind), ncol=ncol(catch))
        catch <- rbind(catch,mat)
        ind <- which(is.na(mergi2$x) & mergi2$midLengths < min(midLengths))
        mat <- matrix(0, nrow=length(ind), ncol=ncol(catch))
        catch <- rbind(mat,catch)


        ## both catch matrices should have same sampling dates
        designMat <- matrix(0, ncol=length(mergi$dates), nrow=length(mergi2$midLengths))
        temp <- designMat
        ind = 1
        for(i in which(!is.na(mergi$x))){
            temp[,i] <- catch[,ind]
            ind <- ind + 1
        }
        catch <- temp
        temp <- designMat
        ind = 1
        for(i in which(!is.na(mergi$y))){
                temp[,i] <- catch2[,ind]
            ind <- ind + 1
        }
        catch2 <- temp

        ## combine lfq data sets
        for(i in 1:dim(designMat)[2]){
            temp1 <- data.frame(midLengths = mergi2$midLengths,
                                catch1 = catch[,i])
            temp2 <- data.frame(midLengths = mergi2$midLengths,
                                catch2 = catch2[,i])
            temp3 <- merge(temp1, temp2, by="midLengths", all=TRUE)
            designMat[,i] <- rowSums(temp3[,c(2,3)])
        }

        ## reassign combined data to vectors
        dates <- mergi$dates
        midLengths <- mergi2$midLengths
        catch <- designMat

    }


    if(!is.na(bin_size)){

        ## error and warning messages
        if(bin_size < midLengths[2]-midLengths[1]) stop("The specified bin_size is smaller than the bin size \nin your data. This is not possible!")

        ## rearrange data into LFQ data
        bin.breaks <- seq(0, max(midLengths) + bin_size, by=bin_size)
        midLengthsNEW <- bin.breaks + bin_size/2
        listi <- vector("list",length(unique(dates)))
        LF_dat <- data.frame(bin = bin.breaks)
        for(i in 1:length(unique(dates))){
            sampli <- unique(dates)[i]
            lengthi <- as.numeric(midLengths)

            if(length(unique(dates)) > 1){
                freqi <- as.numeric(catch[,dates == sampli])
            }else{
                freqi <- as.numeric(catch[dates == sampli])
            }

            bin.breaks2 <- rep(NA, length(bin.breaks))
            for(ii in 1:length(bin.breaks)){
                if(ii == length(bin.breaks)){
                    bin.breaks2[ii] <- length(which(lengthi >= bin.breaks[ii]))
                }else{
                    bin.breaks2[ii] <- length(which(lengthi >= bin.breaks[ii] & lengthi < bin.breaks[ii+1]))
                }
            }
            bin.breaks3 <- rep(bin.breaks, bin.breaks2)
            dati <- aggregate(list(freq=freqi), by=list(bin=bin.breaks3), sum)

            listi[[i]] <- merge(LF_dat, dati, by.x = "bin", all.x =TRUE)[,2]
        }
        catch_mat <- do.call(cbind,listi)
        catch_mat[is.na(catch_mat)] <- 0
        catch <- catch_mat
        midLengths <- midLengthsNEW


        if(any(catch != 0)){

            ## get rid of 0 bins at both ends
            lowRow <- 0
            resi <- TRUE
            while(resi == TRUE){
              lowRow <- lowRow + 1
              resi <- rowSums(catch, na.rm = TRUE)[lowRow] == 0
            }

            upRow <- nrow(catch)
            resi <- TRUE
            while(resi == TRUE){
              resi <- rowSums(catch, na.rm = TRUE)[upRow] == 0
              upRow <- upRow - 1
            }
            upRow <- upRow + 1

            catch <- catch[lowRow:upRow,]
            midLengths <- midLengths[lowRow:upRow]

        }

        ## correct if catch was numeric already
        if(class(lfq$catch) == "numeric"){
            catch <- as.numeric(catch)
        }

      }
      if(vectorise_catch & !is.matrix(catch)){
        stop(paste0("Catch is ", class(catch), ". To vectorise catch, it has to be a matrix."))
      }
    if(vectorise_catch) aggregate = "year"
    if(!is.na(aggregate) & is.matrix(catch)){
        if(aggregate == "year"){
            ## sum numbers per year
            c_sum <- by(t(catch),format(dates,"%Y"), FUN = colSums)
            # rearrange in data frame
            c_list <- lapply(as.list(c_sum), c)
            c_dat <- as.data.frame(c_list)

              if(any(c_dat != 0)){
                  # get rid of 0 bins at both ends
                  lowRow <- 0
                  resi <- TRUE
                  while(resi == TRUE){
                    lowRow <- lowRow + 1
                    resi <- rowSums(c_dat, na.rm = TRUE)[lowRow] == 0
                  }

                  upRow <- nrow(c_dat)
                  resi <- TRUE
                  while(resi == TRUE){
                    resi <- rowSums(c_dat, na.rm = TRUE)[upRow] == 0
                    upRow <- upRow - 1
                  }
                  upRow <- upRow + 1

                  catch <- c_dat[lowRow:upRow,]
                  midLengths <- midLengths[lowRow:upRow]
              }else{
                  catch <- c_dat
              }

            # override old dates
            dates <- unique(as.Date(paste0(format(dates,"%Y"),"-01-01")))
        }else if(aggregate == "quarter"){
            months <- format(dates, "%m")
            seasons <- rep(NA,length(months))
            seasons[months == "01"] <- 2
            seasons[months == "02"] <- 2
            seasons[months == "03"] <- 2
            seasons[months == "04"] <- 5
            seasons[months == "05"] <- 5
            seasons[months == "06"] <- 5
            seasons[months == "07"] <- 8
            seasons[months == "08"] <- 8
            seasons[months == "09"] <- 8
            seasons[months == "10"] <- 11
            seasons[months == "11"] <- 11
            seasons[months == "12"] <- 11
            dateFac <- as.Date(paste0(format(dates,"%Y"),"-",seasons,"-15"))
            ## sum numbers per year
            c_sum <- by(t(catch),dateFac, FUN = colSums)
            # rearrange in data frame
            c_list <- lapply(as.list(c_sum), c)
            c_dat <- as.data.frame(c_list)

              if(any(c_dat != 0)){
                  # get rid of 0 bins at both ends
                  lowRow <- 0
                  resi <- TRUE
                  while(resi == TRUE){
                    lowRow <- lowRow + 1
                    resi <- rowSums(c_dat, na.rm = TRUE)[lowRow] == 0
                  }

                  upRow <- nrow(c_dat)
                  resi <- TRUE
                  while(resi == TRUE){
                    resi <- rowSums(c_dat, na.rm = TRUE)[upRow] == 0
                    upRow <- upRow - 1
                  }
                  upRow <- upRow + 1

                  catch <- c_dat[lowRow:upRow,]
                  midLengths <- midLengths[lowRow:upRow]
              }else{
                  catch <- c_dat
              }

            # override old dates
            dates <- unique(dateFac)
        }else if(aggregate == "month"){
            ## sum numbers per year
            c_sum <- by(t(catch),format(dates,"%Y-%m"), FUN = colSums)
            # rearrange in data frame
            c_list <- lapply(as.list(c_sum), c)
            c_dat <- as.data.frame(c_list)

              if(any(c_dat != 0)){
                  # get rid of 0 bins at both ends
                  lowRow <- 0
                  resi <- TRUE
                  while(resi == TRUE){
                    lowRow <- lowRow + 1
                    resi <- rowSums(c_dat, na.rm = TRUE)[lowRow] == 0
                  }

                  upRow <- nrow(c_dat)
                  resi <- TRUE
                  while(resi == TRUE){
                    resi <- rowSums(c_dat, na.rm = TRUE)[upRow] == 0
                    upRow <- upRow - 1
                  }
                  upRow <- upRow + 1

                  catch <- c_dat[lowRow:upRow,]
                  midLengths <- midLengths[lowRow:upRow]
              }else{
                  catch <- c_dat
              }

            # override old dates
            dates <- unique(as.Date(paste0(format(dates,"%Y-%m"),"-15")))
        }else{
            stop('aggregate has to be either NA, "year", "quarter", or "month"')
        }
    }

      # plus group
      if(isTRUE(plus_group) | is.numeric(plus_group) | plus_group == "Linf"){
        if(isTRUE(plus_group)){
          if(is.vector(catch)){
            print(data.frame(midLengths = midLengths, frequency = catch))
          }else if(length(unique(format(lfq$dates, "%Y"))) == 1){
            print(data.frame(midLengths = midLengths, frequency = rowSums(catch)))
          }else{
            # sum numbers per year
            c_sum <- by(t(catch),format(dates,"%Y"), FUN = colSums)

            # rearrange in data frame
            c_list <- lapply(as.list(c_sum), c)
            c_dat <- as.data.frame(c_list)

            tmp <- data.frame(midLengths = midLengths)
            tmp <- cbind(tmp, c_dat)
            print(tmp)
          }

          writeLines(paste0("Linf = ",round(linf,2),
                            ". Check the table above and insert the length of the plus group (Esc to cancel)."))
          pg = -1
          while(pg > max(midLengths) | pg < min(midLengths)){
            pg <- readline(paste0("Enter a length group between ", min(midLengths)," and ",
                                  max(midLengths),":"))
            pg = as.numeric(as.character(pg))
            if(!(pg %in% midLengths)){
              writeLines(paste0(pg, " is not an element of midLengths (see table)."))
              pg = -1
              #pg <- ifelse(grepl("\\D",pg),-1,as.integer(pg))
              if(is.na(pg)){break}  # breaks when hit enter
            }
          }
        }else if(is.numeric(plus_group)){
            pg = as.numeric(as.character(plus_group))
        }else if(plus_group == "Linf"){
            interval <- midLengths[2] - midLengths[1]
            upperLength <- midLengths + (interval / 2)
            if(!is.na(linf)){
                pg <- midLengths[which.min(abs(upperLength - floor(linf)))]
            }else{
                writeLines("Please provide Linf in par or lfq!")
            }
        }
        if(!(pg %in% midLengths)){
          stop(paste0(pg, " is not an element of midLengths. Set 'plus_group' TRUE and pick a length class \n or check the vector 'midLengths' in your data."))
        }
        midLengths <- midLengths[1:which(midLengths == pg)]
        if(is.vector(catch)){
          if(which(midLengths == pg) < (length(catch)-1)){
            addplus <- sum(catch[((which(midLengths == pg)+1):length(catch))])
          }else if(which(midLengths == pg) == (length(catch)-1)){
            addplus <- catch[(which(midLengths == pg)+1)]
          }else if(which(midLengths == pg) == (length(catch))){
            addplus <- 0
          }
          catch <- catch[1:which(midLengths == pg)]
          catch[which(midLengths == pg)] <-
            catch[which(midLengths == pg)] + addplus
        }else{
          if(which(midLengths == pg) < (nrow(catch)-1)){
            addplus <- colSums(catch[((which(midLengths == pg)+1):nrow(catch)),])
          }else if(which(midLengths == pg) == (nrow(catch)-1)){
            addplus <- catch[(which(midLengths == pg)+1),]
          }else if(which(midLengths == pg) == (nrow(catch))){
            addplus <- 0
          }
          catch <- catch[1:which(midLengths == pg),]
          catch[which(midLengths == pg),] <-
            catch[which(midLengths == pg),] + addplus
        }
      }


        ## combine results
        if(is.vector(catch)){
            catches <- as.vector(catch)
        }else catches <- as.matrix(catch)

    res <- list()
    if("species" %in% names(lfq)) res$species <- lfq$species
    if("stock" %in% names(lfq)) res$stock <- lfq$stock
    res$dates = dates
    res$midLengths = midLengths
    res$catch = catches
    if("comment" %in% names(lfq)) res$comment <- lfq$comment

    ## add growth parameter if known
    if("par" %in% names(lfq)){
        if(class(lfq$par) == "list"){
            res$par <- lfq$par
        }else{
            res$par <- as.list(lfq$par)
        }
    }
    if(!is.null(par)){
        if(class(par) == "list"){
            res$par <- par
        }else{
            res$par <- as.list(par)
        }
    }

    ## add all objects to which in lfq to new lfq list
    idx <- names(lfq)[which(!(names(lfq) %in% names(res)))]
    tmpList <- lfq[which(names(lfq) %in% idx)]
    res <- c(res, tmpList)

    if(class(res) != "lfq"){
        class(res) <- "lfq"
    }
    ## if(!is.na(bin_size)){class(res) <- "lfq"}

    return(res)
}

#' @title Restructuring of length frequency data
#'
#' @description First step of the Electronic LEngth Frequency ANalysis (ELEFAN),
#' which is restructuring length-frequency data (lfq).
#' This is done according to a certain protocol, described by many authors (see
#' Details or References for more information).
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes
#'   \item \strong{dates} dates of sampling times (class Date)
#'   \item \strong{catch} matrix with catches/counts per length class (row) and
#'   sampling date (column)
#' }
#' @param MA number indicating over how many length classes the moving average
#' should be performed (default: 5)
#' @param addl.sqrt additional squareroot transformation of positive values
#' according to Brey et al. (1988) (default: FALSE).
#' Particularly useful if many observations have a low frequency (<10)
#'
#' @examples
#' # data and plot of catch frequencies
#' data(synLFQ4)
#' plot(synLFQ4, Fname="catch")
#'
#' # restructuring and calculation of ASP
#' synLFQ4 <- lfqRestructure(synLFQ4, MA=11)
#' synLFQ4$ASP
#'
#' # plot of restructured scores and fit of soVBGF growth curves
#' plot(synLFQ4)
#' lfqFitCurves(synLFQ4,
#'  par=list(Linf=80, K=0.5, ta=0.25, C=0.75, ts=0),
#'  draw=TRUE
#' )$fASP
#'
#'
#' @details This function is used prior to fitting of growth curves (e.g. in
#' \code{\link{ELEFAN}}, \code{\link{ELEFAN_SA}} functions). It restructures a length
#' frequency data set according to a list of steps to emphasise cohorts in the data.
#' The steps can be found in various publications, see e.g. Brey et al. (1988) or
#'  Pauly and David (1981). Here, the most recent steps documented in Gayanilo (1997)
#'  are followed.
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{rcounts}: restructured frequencies,
#'   \item \strong{peaks_mat}: matrix with uniquely numbered positive peaks,
#'   \item \strong{ASP}: available sum of peaks, sum of posititve peaks which
#'   could be potential be hit by growth curves. This is calculated as the sum of
#'   maximum values from each run of posive restructured scores,
#'   \item \strong{MA}: moving average used for restructuring.
#' }
#'
#'
#' @references
#' Brey, T., Soriano, M., and Pauly, D. 1988. Electronic length frequency analysis:
#' a revised and expanded user's guide to ELEFAN 0, 1 and 2.
#'
#' Gayanilo, Felimon C. FAO-ICLARM stock assessment tools: reference manual.
#' No. 8. Food & Agriculture Org., 1997.
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
#' Pauly, D. and G. Gaschuetz. 1979. A simple method for fitting oscillating length
#' growth data, with a program for pocket calculators. I.C.E.S. CM 1979/6:24.
#' Demersal Fish Cttee, 26 p.
#'
#' Pauly, D. 1984. Fish population dynamics in tropical waters: a manual for use
#' with programmable calculators (Vol. 8). WorldFish.
#'
#' Quenouille, M. H., 1956. Notes on bias in estimation. \emph{Biometrika}, 43:353-360
#'
#' Somers, I. F., 1988. On a seasonally oscillating growth function.
#' ICLARM Fishbyte 6(1): 8-11.
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

lfqRestructure <- function(param, MA=5, addl.sqrt=FALSE){

  lfq <- param

  # replace NAs in catch
  lfq$catch[which(is.na(lfq$catch))] <- 0

  if(MA%%2 == 0) stop("MA must be an odd integer")

  # Steps refer to Gayanilo (1997) FAO-ICLARM stock assessment tools: reference manual
  rcounts <- 0*lfq$catch
  for(i in seq(ncol(lfq$catch))){
    pm <- (MA-1)/2 # plus minus

    # positions of first and last non-zero valules
    val_first <- min(which(lfq$catch[,i] != 0))
    val_last <- max(which(lfq$catch[,i] != 0))
    val_pos <- seq(val_first, val_last)
    val_string <- lfq$catch[val_pos,i]

    # number of values
    n <- length(val_string)

    AF <- NaN*val_string
    nz <- NaN*val_string

    if(n > 1){
        temp <- seq(val_string)
    }else{
        temp <- 1
    }
    ## Steps A & B - Computation of the moving average
    for(j in temp){
      idx <- (j-pm):(j+pm)
      idx <- idx[which(idx %in% temp)]
      idxn <- idx[-which(idx==j)] # neighbors only
      nz[j] <- sum(val_string[idxn] == 0) + (MA-length(idx)) # number of adjacent zeros
      MA.j <- sum(val_string[idx])/MA
      AF[j] <- val_string[j]/MA.j
    }
    # intermediate step to remove Inf or NA
    AF <- replace(AF, which(AF==Inf | is.na(AF)), 0)

    # Calculate mean quotient
    mprime <- mean(AF, na.rm=TRUE)

    ## Step C Divide by mean quotient and subtract 1.0
    Fs <- AF / mprime - 1 # restructured frequencies

    ## Steps D & E - Identify isolated peaks; Adjust for zero frequency
    posFs <- which(Fs > 0)
    if(length(posFs)>0) {Fs[posFs] <- (Fs * 0.5^nz)[posFs]}
    # replace ultimate length bin with zero if negative
    if(sign(Fs[length(Fs)]) == -1){Fs[length(Fs)] <- 0}
    # divide penultimate length bin by 2 if negative
    if(length(sign(Fs[length(Fs)-1])) > 0 &&
       sign(Fs[length(Fs)-1]) == -1){Fs[length(Fs)-1] <- Fs[length(Fs)-1]*0.5}

    ## Step F - Adjust for Fi
    SPV <- sum(Fs[which(Fs > 0)]) # sum of positive values
    SNV <- sum(Fs[which(Fs < 0)]) # sum of negative values
    # set -1 to 0
    minus1 <- which((1+Fs) < 1e-8 | is.na(Fs))
    if(length(minus1)>0) {Fs[minus1] <- 0}
    # adjust negative numbers
    isneg <- which(Fs < 0)
    Fs[isneg] <- Fs[isneg] * (SPV/-SNV)

    # optional square-root adjustment to emphasize larger length bins with lower counts
    if(addl.sqrt){
      posFs <- which(Fs > 0)
      if(length(posFs)>0) {Fs[posFs] <- Fs[posFs] / sqrt(1+2/lfq$catch[posFs,i])} #Fs[posFs] / sqrt(1+2/Fs[posFs])}
    }

    rcounts[val_pos,i] <- Fs
  }

    lfq$rcounts <- rcounts

    # create peak matrix
    prep_mat <- lfq$rcounts
    prep_mat <- ifelse(prep_mat > 0,1,0)
    peaks_mat <- NA*prep_mat
    for(i in seq(ncol(peaks_mat))){
      vec_peaki <- prep_mat[,i]
      runs <- rle(vec_peaki)
      rle_val <- runs$values
      rle_val[which(rle_val == 1)] <- 1:length(rle_val[which(rle_val == 1)])
      peaks_mat[,i] <- rep(rle_val, runs$lengths)
    }
    maxn.peaks <- max(peaks_mat, na.rm=TRUE)
    peaks_mat <- peaks_mat + (prep_mat * maxn.peaks * col(peaks_mat))
    lfq$peaks_mat <- peaks_mat

    # ASP calc
    sampASP <- NaN*seq(ncol(rcounts))

    for(i in seq(ncol(rcounts))){
        ## lfq.i <- lfq[i,]
        tmp <- rle(sign(rcounts[,i]))
        start.idx <- c(1, cumsum(tmp$lengths[-length(tmp$lengths)])+1)
        end.idx <- cumsum(tmp$lengths)
        posrun <- which(tmp$values == 1)
        peakval <- NaN*posrun
        if(length(posrun) > 0){
            for(p in seq(length(posrun))){
                peakval[p] <- max(rcounts[start.idx[posrun[p]]:end.idx[posrun[p]], i ])
            }
            sampASP[i] <- sum(peakval)
        }else{
            sampASP[i] <- 0
        }

    }

    ASP <- sum(sampASP)
    lfq$ASP <- ASP
    lfq$MA <- MA

    class(lfq) <- "lfq"
    return(lfq)
}



#' @title Calculate age-based VBGF parameters from time-based estimates
#'
#' @description Translates ta to t0 and ts to age-adjusted ts. In addition to
#' VBGF parameters calculated via e.g. ELEFAN, the length at recruitment
#' (\code{L0}) must also be provided.
#'
#' @param par list. VBGF parameters.
#' @param L0 numeric. Length at age zero (i.e. recruitment)
#' (default: \code{L0 = 0}).
#' @param plot logical. Plot graphical representation of both
#' time and age based growth curves (default: \code{plot = TRUE}).
#'
#' @return a list containging age-based VBGF parameters.
#' @export
#'
#' @examples
#' data(synLFQ4)
#' lfq <- synLFQ4
#' lfq$par <- list(
#'   Linf = 80, K = 0.5, ta = 0.25, C = 0.75, ts = 0.5
#' )
#' ta2t0(par = lfq$par, L0=10)
#'
ta2t0 <- function(par = NULL, L0 = 0, plot = TRUE){
  if(is.null(par)){
    stop("Error: must provide 'par'")
  }

  if(L0 > par$Linf){stop("Error: 'L0' must be lower than 'par$Linf'")}

  # add seasonalized pars if needed
  seasonalized <- !is.null(par$C)
  if(!seasonalized){par$ts <- 0; par$C <- 0}

  t <- seq(0, 3, 0.01)
  Lt <- VBGF(par = par, t = t)
  trecr <- t[which.min(sqrt((Lt - L0)^2))]
  t0 <- par$ta - trecr
  if(seasonalized){
    ts <- (par$ts - trecr)%%1
  } else {
    ts <- 0
  }
  par_age <- par
  par_age$ta <- NULL
  par_age$t0 <- t0
  par_age$ts <- ts
  age <- seq(t0, 3, 0.01)
  Lage <- VBGF(par = par_age, t = age)

  if(plot){
    op <- par(mfcol = c(1,2), mar = c(4,4,1,1), mgp = c(2,0.5,0), cex = 1)
    # plot by time
    plot(t, Lt, t="l", xlab = "time", ylim = c(0, max(Lt)*1.1), yaxs="i", ylab = bquote(L[t]))
    ylim <- par()$usr[3:4]
    grid(); abline(v = trecr, lty=3); box()
    lines(t, Lt)
    points(trecr, L0, pch = 1, col = 2)
    points(trecr + t0, 0, pch = 4, col = 4)
    tmp <- data.frame(var = c(names(par), "trecr", "L0"), val = c(unlist(par), trecr, L0))
    tmp$pch = NA; tmp$col <- NA
    is_t0 <- which(tmp$var=="ta")
    tmp$pch[is_t0] <- 4
    tmp$col[is_t0] <- 4
    is_L0 <- which(tmp$var=="L0")
    tmp$pch[is_L0] <- 1
    tmp$col[is_L0] <- 2
    legend("bottomright",
      legend = c(paste(tmp$var, "=", round(tmp$val,2))),
      bty = "n", cex = 0.7, pt.cex = 1,
      pch = tmp$pch, col = tmp$col
    )
    # plot by age
    plot(age, Lage, t="l", ylim = ylim, yaxs = "i", ylab = bquote(L[age]))
    grid(); abline(v = 0, lty=3); box()
    lines(age, Lage)
    points(0, L0, pch = 1, col = 2)
    points(t0, 0, pch = 4, col = 4)
    tmp <- data.frame(var = c(names(par_age), "L0"), val = c(unlist(par_age), L0))
    tmp$pch = NA; tmp$col <- NA
    is_t0 <- which(tmp$var=="t0")
    tmp$pch[is_t0] <- 4
    tmp$col[is_t0] <- 4
    is_L0 <- which(tmp$var=="L0")
    tmp$pch[is_L0] <- 1
    tmp$col[is_L0] <- 2
    legend("bottomright",
      legend = c(paste(tmp$var, "=", round(tmp$val,2))),
      bty = "n", cex = 0.7, pt.cex = 1,
      pch = tmp$pch, col = tmp$col
    )
    par(op)
  }
  return(par_age)
}



#' Convert length-frequency data into aggregate numbers by relative age and
#' (optionally) cohort
#'
#' @param lfq an object of class lfq
#' @param method method of conversion (LCCC, GOTCHA, SLICC)
#' @param agemax integer. Maximum age plus group
#' @param n.cohort integer. Number of pseudocohort to derive from lfq.
#'   Only used GOTCHA method (Defaults to number of length bins for
#'   consistency with LCCC).
#' @param n.cohort.per.yr integer. Number of pseudocohorts per year.
#'   Overrides `n.cohort`. Can be applied to GOTCHA and SLICC. Default
#' @param use.ndt
#'
#' @return list.
#' @export
#'
#' @examples
#'
#' data("synLFQ4")
#' lfq <- synLFQ4
#' lfq$par <- list(Linf = 80, K = 0.5, C = 0.75, ts = 0.5, ta = 0.25)
#' lfq <- lfqModify(lfq, year = 2006, bin_size = 2)
#' lfq <- lfqRestructure(lfq)
#' plot(lfq)
#'
#' # LCCC
#' res <- catchCurvePrep(lfq = lfq, method = "LCCC")
#' plot(log(n) ~ rel.age, res$tab)
#' incl <- which(res$tab$rel.age > 1)
#' points(log(n) ~ rel.age, res$tab[incl,], pch = 20 )
#' fit <- lm(log(n) ~ rel.age, res$tab[incl,])
#' abline(fit)
#' -coef(fit)[2] # true: Z = 1.0
#'
#' # GOTCHA (default settings)
#' res <- catchCurvePrep(lfq = lfq, method = "GOTCHA")
#' plot(log(n) ~ rel.age, res$tab)
#' incl <- which(res$tab$rel.age > 2 & res$tab$rel.age < 6)
#' points(log(n) ~ rel.age, res$tab[incl,], pch = 20 )
#' fit <- lm(log(n) ~ rel.age, res$tab[incl,])
#' abline(fit)
#' -coef(fit)[2] # true: Z = 1.0
#'
#' # GOTCHA (1 cohort per year)
#' res <- catchCurvePrep(lfq = lfq, method = "GOTCHA", n.cohort.per.yr = 1)
#' plot(log(n) ~ rel.age, res$tab)
#' incl <- which(res$tab$rel.age > 1.5 & res$tab$rel.age < 6)
#' points(log(n) ~ rel.age, res$tab[incl,], pch = 20 )
#' fit <- lm(log(n) ~ rel.age, res$tab[incl,])
#' abline(fit)
#' -coef(fit)[2] # true: Z = 1.0
#'
#' # SLICC
#' res <- catchCurvePrep(lfq = lfq, method = "SLICC")
#' plot(log(n) ~ rel.age, res$tab, col = res$tab$cohort)
#' incl <- which(res$tab$rel.age > 1 & res$tab$rel.age < 6)
#' points(log(n) ~ rel.age, res$tab[incl,], pch = 20,
#'   col = res$tab$cohort[incl])
#' fit0 <- lm(log(n) ~ rel.age + cohort, res$tab[incl,])
#' fit1 <- lm(log(n) ~ rel.age, res$tab[incl,])
#' fit <- get(c("fit0", "fit1")[which.min(AIC(fit0, fit1)$AIC)])
#' abline(fit)
#' -coef(fit)[2] # true: Z = 1.0
#'
#'
#'
catchCurvePrep <- function(
  lfq,
  method = "LCCC",
  agemax = NULL,
  n.cohort = NULL, # only applies to GOTCHA1
  n.cohort.per.yr = NULL, # this should override n.cohort in GOTCHA
  use.ndt = NULL
  ){

  # replace zeros in catch mat
  lfq$catch[which(is.na(lfq$catch))] <- 0

  if(is.null(agemax)) agemax <- ceiling(VBGF(pars = lfq$par, L = lfq$par$Linf*0.95))

  if(is.null(use.ndt)){
    if(method == "LCCC"){
      use.ndt <- TRUE
    }else{
      use.ndt <- FALSE
    }
  }

  if(method == "LCCC"){
    lfqx <- lfq
    lfqx$par$C <- 0 # remove seasonality
    if(is.null(n.cohort.per.yr)){n.cohort.per.yr <- 36}
    lfqx <- lfqCohort(lfqx, n.cohort.per.yr = n.cohort.per.yr, calc_dt = TRUE, agemax = agemax)
    sumtab = data.frame(rel.age = apply(lfqx$rel.age, 1, mean, na.rm = TRUE))
    if(use.ndt){
      sumtab$n <- apply(lfqx$catch/lfqx$dt, 1, sum, na.rm = TRUE)
    }else{
      sumtab$n <- apply(lfqx$catch, 1, sum, na.rm = TRUE)
    }
    sumtab$length <- apply(array(lfqx$midLengths, dim = dim(lfqx$catch)), 1, mean, na.rm = TRUE)
    sumtab <- sumtab[which(sumtab$n > 0),] # remove rows where n == 0 | NA
  }

  if(method == "GOTCHA"){
    # if n.cohort.per.yr not NULL, then use this criteria for slicing
    if(!is.null(n.cohort.per.yr)) n.cohort <- NULL
    # if NULL, then use thin slices (n=36 per year) and aggregate later to n.cohort bins
    if(is.null(n.cohort.per.yr)){
      n.cohort.per.yr <- 36
      if(is.null(n.cohort)) n.cohort <- dim(lfq$catch)[1]
    }
    lfqx <- lfq
    lfqx <- lfqCohort(lfqx, calc_dt = TRUE, n.cohort.per.yr = n.cohort.per.yr, agemax = agemax)
    df <- data.frame(
      length = rep(lfqx$midLengths, times = length(lfqx$dates)),
      rel.age = c(lfqx$rel.age),
      bday = c(lfqx$bday),
      cohort = factor(c(lfqx$cohort)),
      dt = c(lfqx$dt))
    if(use.ndt){
      df$n <- c(lfqx$catch)/c(lfqx$dt)
    }else{
      df$n <- c(lfqx$catch)
    }
    df <- df[which(df$n > 0),] # remove rows where n == 0 | NA
    # create new slice aggregates if n.cohort is defined
    if(!is.null(n.cohort)){
      breaks <- seq(min(df$bday), max(df$bday), length.out = n.cohort)
      mids <- breaks[-1] - (diff(breaks)[1]/2)
      df$bday.cat <- cut(df$bday, breaks = breaks)
      df$bday2 <- mids[df$bday.cat]
      df$rel.age2 <- max(date2yeardec(lfqx$dates)) - df$bday2
      agg <- aggregate(n ~ rel.age2, df, sum, na.rm = TRUE)
      agg2 <- aggregate(length ~ rel.age2, df, max, na.rm = TRUE)
      sumtab <- data.frame(rel.age = agg$rel.age2, n = agg$n, length = agg2$length)
    }
    if(is.null(n.cohort)){
      agg <- aggregate(n ~ bday, df, sum, na.rm = TRUE)
      agg$rel.age <- max(date2yeardec(lfqx$dates)) - agg$bday
      agg2 <- aggregate(length ~ bday, df, max, na.rm = TRUE)
      sumtab <- data.frame(rel.age = agg$rel.age, n = agg$n, length = agg2$length)
    }
  }

  if(method == "SLICC"){
    lfqx <- lfq
    if(!is.null(n.cohort)) n.cohort <- NULL # always set n.cohort to NULL
    if(is.null(n.cohort.per.yr)) n.cohort.per.yr <- 1 # set to 1 cohort per year as default
    lfqx <- lfqCohort(lfqx, calc_dt = TRUE, n.cohort.per.yr = n.cohort.per.yr, agemax = agemax)
    df <- data.frame(
      length = rep(lfqx$midLengths, times = length(lfqx$dates)),
      rel.age = c(lfqx$rel.age),
      bday = c(lfqx$bday),
      cohort = factor(c(lfqx$cohort)),
      dt = c(lfqx$dt))
    if(use.ndt){
      df$n <- c(lfqx$catch)/c(lfqx$dt)
    }else{
      df$n <- c(lfqx$catch)
    }
    df <- df[which(df$n > 0),] # remove rows where n == 0 | NA
    agg <- aggregate(n ~ cohort + rel.age, df, sum, na.rm = TRUE)
    agg2 <- aggregate(length ~ cohort + rel.age, df, mean, na.rm = TRUE)
    sumtab <- data.frame(rel.age = agg$rel.age, cohort = agg$cohort, n = agg$n, length = agg2$length)
  }

  res <- list(tab = sumtab, method = method, agemax = agemax,
    n.cohort = NULL, n.cohort.per.yr = n.cohort.per.yr,
    use.ndt = use.ndt)

  return(res)

}






#' Automated determination of values to include in catch curve
#'
#' @description Procedure described by Pauly (1990) to automate the selection
#'   of points used in a catch curve (\code{log(y)~x}). The sequential
#'   data indices that result in the highest F-ststistic in a linear regression
#'   are returned. Results should still be visualized to prevent suboptimal
#'   results form 'pathological' datasets.
#'
#' @param x vector. Age or relative age associated with catch numbers
#'   (\code{x}).
#' @param y vector. catch numbers.
#' @param minN numeric. Minimum number of values to include in the regression
#'   (Default = 3).
#'
#' @return list containing the indices and F-statistic of the highest scoring
#' subset of (sequantial) data points.
#'
#' @export
#'
#' @references
#' Pauly, D. (1990). Length-converted catch curves and the seasonal growth
#'   of fishes. Fishbyte, 8(3), 33–38.
#'
#' @examples
#' # generate example data
#' set.seed(1)
#' n <- 10
#' x <- seq(n)
#' logy <- 10 + -2*x + rnorm(n, sd = 0.5)
#' logy[1:2] <- logy[1:2] * c(0.5, 0.75)
#' y <- exp(logy)
#' plot(y ~ x, log = "y")
#'
#' # When lm can reduce to 3 points, this is chosen as the 'best'
#' (tmp <- bestCC(x = x, y = y, minN = 3))
#' plot(f~n, data = tmp$bestbyn, t = "b", log="y")
#' plot(y ~ x, log = "y")
#' points(x[tmp$best$samples], y[tmp$best$samples], pch = 16)
#'
#' # Any min number of points above this threshhold increases the best n a lot
#' (tmp <- bestCC(x = x, y = y, minN = 4))
#' plot(f~n, data = tmp$bestbyn, t = "b", log="y")
#' plot(y ~ x, log = "y")
#' points(x[tmp$best$samples], y[tmp$best$samples], pch = 16)
#'
bestCC <- function(x, y, minN = max(3,0.25*length(x))){

  if(minN < 3){stop("'minN' must be greater than or equal to 3")}

  # return unique combinations of sequential data points
  res <- vector("list", length(minN:length(x)))
  for(i in seq(minN:length(x))){
    len <- (minN:length(x))[i]
    res.i <- vector("list", length(1:(length(x)-len+1)))
    for(j in seq(1:(length(x)-len+1))){
      start <- (1:(length(x)-len+1))[j]
      res.i[[j]] <- seq(from = start, by = 1, length.out = len)
    }
    res[[i]] <- res.i
  }
  res <- do.call("c", res)
  res <- unique(res)

  # fit lm to each comb and return the indices and f-statistic
  res2 <- vector("list", length(res))
  for(i in seq(res2)){
    res2[[i]]$n <- length(res[[i]])
    res2[[i]]$samples <- res[[i]]
    df.i <- data.frame(y = y[res[[i]]], x = x[res[[i]]])
    fit <- lm(log(y) ~ x, data = df.i)
    S <- summary(fit)
    res2[[i]]$f <- 0
    if(coef(fit)[2] < 0 & S$coefficients[2,4] < 0.05){res2[[i]]$f <- S$fstatistic[1]}
  }

  bestbyn <- do.call("rbind", lapply(res2, FUN = function(x){data.frame(n=x$n,f=x$f)}))
  bestbyn <- aggregate(f ~ n, data = bestbyn, FUN = max)

  # 'best' model is that with highest f-statistic
  best <- which.max(do.call("c", lapply(res2, FUN = function(x){x$f})))

  return(list(
    best = res2[[best]],
    bestbyn = bestbyn
  ))
}

