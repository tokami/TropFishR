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

lfqCreate <- function(data, Lname, Dname, Fname = NA, bin_size = 2,
                    length_unit = "cm", plus_group = FALSE,
                    aggregate_dates = FALSE,
                    plot = FALSE){

  data$length <- get(Lname, data)
  data$date <- get(Dname, data)
  if(!is.na(Fname)){
    data$freq <- get(Fname, data)
  }
  if(class(data$date) != "Date") stop(noquote("Please provide the date as 'Date' class (e.g. as.Date())."))

  # convert length if necessary
  if(length_unit == "m") data$length <- data$length * 100
  if(length_unit == "mm") data$length <- data$length / 10

  # delete any length of 0
  data$length[which(data$length == 0)] <- NA

  # raise lengths according to frequencies
  if(!is.na(Fname)){
    data2 <- data[rep.int(seq_len(nrow(data)),times=data$freq),]
  }else{
    data2 <- data
  }

  # show histogram
  ## hist(data2$length, breaks = 50)

  # no NAs allowed in length nor date column
  data2 <- data2[!is.na(data2$length),]
  data2 <- data2[!is.na(data2$date),]

  # order according to date
  data2 <- data2[order(data2$date),]

  # create monthly grouped dates
  if(aggregate_dates){
    samplings <- as.Date(paste(format(data2$date, "%Y-%m"),"15",sep="-"))
  }else samplings <- data2$date

  # rearrange data into LFQ data
  bin.breaks <- seq(0, max(data2$length) + bin_size, by=bin_size)
  midLengths <- bin.breaks[-length(bin.breaks)] + bin_size/2
  listi <- vector("list",length(unique(samplings)))
  LF_dat <- data.frame(midLengths = midLengths)
  for(i in 1:length(unique(samplings))){
    sampli <- unique(samplings)[i]
    dati <- as.character(unique(samplings)[i])
    lengthi <- as.numeric(data2$length[samplings == sampli])
    cuti <- cut(lengthi, breaks = bin.breaks, labels = midLengths, include.lowest = TRUE)
    freq <- plyr::count(cuti)
    colnames(freq) <- c("midLengths", dati)
    listi[[i]] <- merge(LF_dat,freq, by.x = "midLengths", all.x =TRUE)[,2]
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



  res <- list(dates = unique(samplings),
                   midLengths = midLengths,
                   catch = catch_mat)
  class(res) <- "lfq"
  if(plot) plot(res, Fname = "catch")
  return(res)
}
