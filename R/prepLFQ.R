#' @title Prepare lfq data for further analysis
#'
#' @description Rearrange catch matrix in length frequency data (lfq class) to
#'    have one column per year. This is required for e.g. \code{\link{catchCurve}}.
#'
#' @param lfq lfq object with dates, midLengths, and catch
#' @param par growth parameter as resulting from e.g. \code{\link{ELEFAN}}
#'
#' @keywords function lfq length-frequency
#'
#' @examples
#' data(synLFQ4)
#' newlfq <- prepLFQ(synLFQ4)
#'
#' @return lfq object with rearranged catch matrix (yearly sums) and growth parameters
#'    if provided.
#'
#' @export

prepLFQ <- function(lfq, par = NULL){

  dates <- lfq$dates
  midLengths <- lfq$midLengths
  catch <- lfq$catch

  # sum numbers per year
  c_sum <- by(t(catch),format(dates,"%Y"),FUN = colSums)

  # rearrange in data frame
  c_list <- lapply(as.list(c_sum), c)
  c_dat <- as.data.frame(c_list)

  # get rid of 0 bins at both ends
  lowRow <- 1
  resi <- TRUE
  while(resi == TRUE){
    resi <- rowSums(c_dat)[lowRow] == 0
    lowRow <- lowRow + 1
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

  # combine results
  res <- list(dates = dates,
              midLengths = midLengths,
              catch = as.matrix(catch))

  # add growth parameter if knowm
  if(!is.null(par)){
    res$Linf <- par$Linf
    res$K <- par$K
    res$t_anchor <- par$t_anchor
  }

  return(res)
}

