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
#' @param plus_group logical or numeric; should a plus group be created? If yes you will be
#'    asked to insert the length for the plus group in the console (default: FALSE).
#'    Instead of inserting the length of the plus group via the console, the value
#'    can be inserted, e.g. plus_group = 85.5.
#' @param lfq2 optional second lfq object which will be merged with lfq. This might be interesting for
#'    fleet specific lfq objects. Default: NA.
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
                      vectorise_catch = FALSE,
                      plus_group = FALSE,
                      lfq2 = NA){

    if(class(lfq) != "lfq") stop("Your lfq data set has to have class 'lfq'!")
    dates <- lfq$dates
    midLengths <- lfq$midLengths
    catch <- lfq$catch

    ## replace NAs in catch
    catch[which(is.na(catch))] <- 0


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
    # rearrange data into LFQ data
    bin.breaks <- seq(0, max(midLengths) + bin_size, by=bin_size)
    midLengthsNEW <- bin.breaks[-length(bin.breaks)] + bin_size/2
    listi <- vector("list",length(unique(dates)))
    LF_dat <- data.frame(midLengths = midLengthsNEW)
    for(i in 1:length(unique(dates))){
      sampli <- unique(dates)[i]
      dati <- as.character(unique(dates)[i])
      if(length(unique(dates)) > 1){
        lengthi <- rep.int(midLengths,times=as.numeric(catch[,dates == sampli]))
      }else{
        lengthi <- rep.int(midLengths,times=as.numeric(catch[dates == sampli]))
      }
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
  if(isTRUE(plus_group) | is.numeric(plus_group)){
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
    }else if(is.numeric(plus_group)){
      pg = as.numeric(as.character(plus_group))
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
    
    res <- list(dates = dates,
                midLengths = midLengths,
                catch = catches)

    # add growth parameter if known
    if(!is.null(par)){
      res$Linf <- par$Linf
      res$K <- par$K
      res$t_anchor <- par$t_anchor
    }
    if("par" %in% names(lfq)) res$par <- lfq$par

    if("Linf" %in% names(lfq)) res$Linf <- lfq$Linf
    if("K" %in% names(lfq)) res$K <- lfq$K
    if("t0" %in% names(lfq)) res$t0 <- lfq$t0
    if("t_anchor" %in% names(lfq)) res$t_anchor <- lfq$t_anchor
    if("C" %in% names(lfq)) res$C <- lfq$C
    if("ts" %in% names(lfq)) res$ts <- lfq$ts
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


    if(class(res) != "lfq"){
        class(res) <- "lfq"
    }
    ## if(!is.na(bin_size)){class(res) <- "lfq"}
    
    return(res)
}
