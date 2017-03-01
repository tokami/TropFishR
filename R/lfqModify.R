#' @title Modify lfq data for further analysis
#'
#' @description Modify length-freqeuncy (LFQ) data. Allows to summarise catch matrix
#'    of LFQ data to one column per year. This is required for e.g. \code{\link{catchCurve}}.
#'    Allows to change bin size of LFQ data. Allows to ad plus group to catch matrix.
#'
#' @param lfq lfq object with dates, midLengths, and catch
#' @param par growth parameter as resulting from e.g. \code{\link{ELEFAN}}
#' @param bin_size Bin size for length frequencies (in cm)
#' @param vectorise_catch logical; indicating if the catch matrix should be summarised to
#'    yearly vectors (default: FALSE).
#' @param plus_group logical; should a plus group be created? If yes you will be
#'    asked to insert the length for the plus group in the console (default: FALSE).
#'    Instead of inserting the length of the plus group via the console, the value
#'    can be incorporated in a vector, e.g. plus_group = c(TRUE, 30).
#'
#' @keywords function lfq length-frequency
#'
#' @examples
#' data(synLFQ4)
#'
#' ## summarise catch matrix per year
#' lfq_sum <- lfqModify(synLFQ4)
#'
#' ## change bin size
#' lfq_bin <- lfqModify(synLFQ4, bin_size = 4)
#'
#' @return lfq object with rearranged catch matrix (yearly sums) and growth parameters
#'    if provided.
#'
#' @export

lfqModify <- function(lfq, par = NULL, bin_size = NA, vectorise_catch = FALSE, plus_group = FALSE){

  dates <- lfq$dates
  midLengths <- lfq$midLengths
  catch <- lfq$catch

  if(!is.na(bin_size)){
    # rearrange data into LFQ data
    bin.breaks <- seq(0, max(midLengths) + bin_size, by=bin_size)
    midLengthsNEW <- bin.breaks[-length(bin.breaks)] + bin_size/2
    listi <- vector("list",length(unique(dates)))
    LF_dat <- data.frame(midLengths = midLengthsNEW)
    for(i in 1:length(unique(dates))){
      sampli <- unique(dates)[i]
      dati <- as.character(unique(dates)[i])
      lengthi <- rep.int(midLengths,times=as.numeric(catch[,dates == sampli]))
      cuti <- cut(lengthi, breaks = bin.breaks, labels = midLengthsNEW, include.lowest = TRUE)
      freq <- plyr::count(cuti)
      colnames(freq) <- c("midLengths", dati)
      listi[[i]] <- merge(LF_dat,freq, by.x = "midLengths", all.x =TRUE)[,2]
    }
    catch_mat <- do.call(cbind,listi)
    catch_mat[is.na(catch_mat)] <- 0
    catch <- catch_mat
    midLengths <- midLengthsNEW
  }
  if(vectorise_catch & !is.matrix(catch)){
    stop(paste0("Catch is ", class(catch), ". To vectorise catch, it has to be a matrix."))
  }
  if(vectorise_catch & is.matrix(catch)){
    # sum numbers per year
    c_sum <- by(t(catch),format(dates,"%Y"), FUN = colSums)

    # rearrange in data frame
    c_list <- lapply(as.list(c_sum), c)
    c_dat <- as.data.frame(c_list)

    # get rid of 0 bins at both ends
    lowRow <- 0
    resi <- TRUE
    while(resi == TRUE){
      lowRow <- lowRow + 1
      resi <- rowSums(c_dat)[lowRow] == 0
    }

    upRow <- nrow(c_dat)
    resi <- TRUE
    while(resi == TRUE){
      resi <- rowSums(c_dat)[upRow] == 0
      upRow <- upRow - 1
    }
    upRow <- upRow + 1

    catch <- c_dat[lowRow:upRow,]
    midLengths <- midLengths[lowRow:upRow]

    # override old dates
    dates <- unique(as.Date(paste0(format(dates,"%Y"),"-01-01")))
  }

  # plus group
  if(plus_group[1]){
    if(length(plus_group) == 1){
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
    if(is.vector(catch)){
      addplus <- sum(catch[((which(midLengths == pg)+1):length(catch))])
      catch <- catch[1:which(midLengths == pg)]
      catch[which(midLengths == pg)] <-
        catch[which(midLengths == pg)] + addplus
    }else{
      addplus <- colSums(catch[((which(midLengths == pg)+1):nrow(catch)),])
      catch <- catch[1:which(midLengths == pg),]
      catch[which(midLengths == pg),] <-
        catch[which(midLengths == pg),] + addplus
    }
  }

  # combine results
  if(is.vector(catch)){
    catches <- as.vector(catch)
  }else catches <- as.matrix(catch)
  res <- list(dates = dates,
              midLengths = midLengths,
              catch = catches)
  if(!is.na(bin_size)){class(res) <- "lfq"}

  # add growth parameter if known
  if(!is.null(par)){
    res$Linf <- par$Linf
    res$K <- par$K
    res$t_anchor <- par$t_anchor
  }

  if("Linf" %in% names(lfq)) res$Linf <- lfq$Linf
  if("K" %in% names(lfq)) res$K <- lfq$K
  if("t0" %in% names(lfq)) res$t0 <- lfq$t0
  if("t_anchor" %in% names(lfq)) res$t_anchor <- lfq$t_anchor
  if("M" %in% names(lfq)) res$M <- lfq$M
  if("Z" %in% names(lfq)) res$Z <- lfq$Z
  if("FM" %in% names(lfq)) res$FM <- lfq$FM
  if("E" %in% names(lfq)) res$E <- lfq$E
  if("FM_calc" %in% names(lfq)) res$FM_calc <- lfq$FM_calc
  if("a" %in% names(lfq)) res$a <- lfq$a
  if("b" %in% names(lfq)) res$b <- lfq$b
  if("L50" %in% names(lfq)) res$L50 <- lfq$L50
  if("L75" %in% names(lfq)) res$L75 <- lfq$L75
  if("L95" %in% names(lfq)) res$L95 <- lfq$L95
  if("s_list" %in% names(lfq)) res$s_list <- lfq$s_list

  return(res)
}
