#' @title Jack-knife method for confidence intervals for growth parameters
#
#' @description The jackknife routine can be used to estimate confidence intervals for
#'   the growth parameters estimated with ELEFAN.
#'
#' @param elefanFit lfq list fitted with growth parameters estimated with ELEFAN_SA or ELEFAN_GA
#' 
#' @keywords function jackknife CI
#'
#' @examples
#' 
#' @return 
#'
#' @references
#' 
#' @export

jackknife <- function(elefanFit){

    res <- elefanFit

    if(class(res) != "lfq" | !("par" %in% names(res)))
        stop("Please provide an object fitted by one of the ELEFAN functions!")

    if("control" %in% names(res)){
        method <- "SA"
    }else method <- "GA"

    JK <- vector("list", length(res$dates))
    rns <- vector("numeric", length(res$dates))    
    for(i in 1:length(res$dates)){
      loop_data <- list(dates = res$dates[-i],
                      midLengths = res$midLengths,
                      catch = res$catch[,-i])
      if(method == "GA"){
          tmp <- invisible(capture.output(ELEFAN_GA(loop_data, MA = res$MA, addl.sqrt = res$addl.sqrt,
                           seasonalised = res$seasonalised,
                           popSize = res$popSize,
                           maxiter = res$maxiter,
                           run = res$run,
                           pmutation = res$pmutation,
                           pcrossover = res$pcrossover,
                           elitism  = res$elitism,
                           flagging.out = res$flagging.out,
                           seed = res$seed,
                           low_par = res$low_par,
                           up_par = res$up_par,
                           plot.score = FALSE)))
      }
      if(method =="SA"){
          SA_time <- res$control$SA_time
          maxit <- res$control$maxit
          nb.stop.improvement <- res$control$nb.stop.improvement
          SA_temp <- res$control$temperature
          verbose <- FALSE
          tmp <- invisible(capture.output(ELEFAN_SA(loop_data,
                           SA_time = SA_time, SA_temp = SA_temp,
                           maxit = maxit,
                           nb.stop.improvement = nb.stop.improvement,
                           verbose = verbose,
                           MA = res$MA, addl.sqrt = res$addl.sqrt,
                           init_par = res$init_par,
                           low_par = res$low_par,
                           up_par = res$up_par,
                           seasonalised = res$seasonalised,
                           flagging.out = res$flagging.out,
                           plot.score = FALSE)))
      }
        JK[[i]] <- unlist(tmp$par)
        rns[i] <- tmp$Rn_max
    }
    
    JKres <- do.call(cbind, JK)
    # mean
    JKmeans <- apply(as.matrix(JKres), MARGIN = 1, FUN = mean)
    # confidence intervals
    JKconf <- apply(as.matrix(JKres), MARGIN = 1, FUN = function(x) t.test(x)$conf.int[c(1,2)])
    JKconf <- t(JKconf)
    colnames(JKconf) <- c("lo","up")

    ret <- list(RnScore = rns,
                jkMean = JKmeans,
                jkCI = t(JKconf))

    return(ret)
}

