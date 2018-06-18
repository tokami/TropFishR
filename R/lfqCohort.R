#' Assign cohort and relative age to bins
#'
#' @param lfq an lfq object with fitted VBGF parameters in `lfq$par` slot.
#'
#' @return a list of class `lfq` with additional list elements "cohort" and "rel.age"
#' 
#'
#' @examples
#' data(synLFQ4)
#' synLFQ4 <- lfqRestructure(synLFQ4, MA=11)
#' synLFQ4$par <- list(Linf = 80, K = 0.5, t_anchor = 0.25, C = 0.75, ts = 0.5)
#'
#' plot(synLFQ4, Fname = "rcounts", ylim = c(0, max(synLFQ4$midLengths)))
#' PAR2 <- synLFQ4$par
#' PAR2$t_anchor <- (PAR2$t_anchor - 0.5) %% 1
#' tmp <- lfqFitCurves(synLFQ4, par = PAR2,
#'   draw = TRUE, col=8, lty=1
#' )
#' synLFQ4 <- lfqCohort(synLFQ4)
#' pal <- colorRampPalette(c(4,5,7,2))
#' with(synLFQ4, image(x = dates, y = midLengths, z = t(cohort),
#'   col = adjustcolor(pal(max(cohort, na.rm = TRUE)), 0.5),
#'   add=TRUE
#' ))
#'
#' ### Use all data in catch curve
#' df <- data.frame(
#'   rel.age = c(synLFQ4$rel.age[,-25]),
#'   n = c(synLFQ4$catch[,-25]),
#'   cohort = c(synLFQ4$cohort[,-25])
#' )
#' df$age <- floor(df$rel.age)
#'
#' # floor of rel.age
#' agg <- aggregate(n ~ age, data = df, FUN = sum)
#' plot(log(n) ~ age, agg)
#' fit <- lm(log(n) ~ age, subset(agg, age %in% seq(1,4)))
#' abline(fit)
#' coef(fit)[2]
#'
#' # direct rel.age
#' agg <- aggregate(n ~ rel.age, data = df, FUN = sum)
#' plot(log(n) ~ rel.age, agg)
#' fit <- lm(log(n) ~ rel.age, subset(agg, rel.age > 0.75 & rel.age < 3))
#' abline(fit)
#' -coef(fit)[2]
#'
#'
#' ### For a single cohort
#' dfsub <- subset(df, cohort == 6)# most complete cohort
#'
#' # direct rel.age
#' agg <- aggregate(n ~ rel.age, data = dfsub, FUN = sum)
#' plot(log(n) ~ rel.age, agg)
#' fit <- lm(log(n) ~ rel.age, subset(agg, rel.age > 1 & rel.age < 3))
#' abline(fit)
#' -coef(fit)[2]
#'
#' @export
#' 
lfqCohort <- function(lfq){

    if(!("par" %in% names(lfq)) & (!("Linf" %in% names(lfq)) | !("K" %in% names(lfq)))){
        stop("This method requires growth parameter!")
    }else if(!("par" %in% names(lfq))){
        PAR1 <- list(Linf = lfq$Linf, K = lfq$K, t_anchor = lfq$t_anchor)
    }else if("par" %in% names(lfq)){
        PAR1 <- lfq$par
    }
    PAR2 <- PAR1

    if(length(lfq$par$t_anchor) == 1){
        PAR2$t_anchor <- (lfq$par$t_anchor - 0.5) %% 1      
        lfq1 <- lfqFitCurves(lfq, par = PAR1, draw = FALSE)
        lfq2 <- lfqFitCurves(lfq, par = PAR2, agemax = lfq1$agemax+0.5, draw = FALSE)
    }else if(length(lfq$par$t_anchor) == 2){
        taDiff1 <- (lfq$par$t_anchor[2] - lfq$par$t_anchor[1])/2 + lfq$par$t_anchor[1]
        taDiff2 <- ((lfq$par$t_anchor[1]+1) - lfq$par$t_anchor[2])/2 + lfq$par$t_anchor[2]
        PAR2$t_anchor <- c(taDiff1,taDiff2)
        lfq1 <- lfqFitCurves(lfq, par = PAR1, draw = FALSE, spawningTimes = 2)        
        lfq2 <- lfqFitCurves(lfq, par = PAR2, spawningTimes = 2,
                             agemax = lfq1$agemax+max(c(taDiff1,taDiff2)), draw = FALSE)
        PAR2a <- PAR2
        PAR2a$t_anchor <- taDiff1
        PAR2b <- PAR2
        PAR2b$t_anchor <- taDiff2   
        lfq2a <- lfqFitCurves(lfq, par = PAR2a, draw = FALSE)  
        lfq2b <- lfqFitCurves(lfq, par = PAR2b, draw = FALSE)
    }else{
        stop("t_anchor has to have length 1 or 2! Not implemented yet for more than 2 spawning events per year!")
    }

    rel.age <- lfq$catch*NaN
    cohort <- lfq$catch*NaN
    if(length(lfq$par$t_anchor) == 1){
        for(i in seq(length(lfq1$dates))){
            t.use <- which(lfq2$Lt$t == date2yeardec(lfq1$dates[i]))
            Lt.use <-  lfq2$Lt[t.use,]
            ## upper cohort bound
            upper <- suppressWarnings(apply(
              outer(X = lfq1$midLengths, Y = Lt.use$Lt, FUN = "-"),
              MARGIN = 1,
              FUN = function(x){max(which(x < 0))}
            ))
            ## If no upper exists, match to the oldest cohort (i.e. row 1 of Lt.use)
            upper <- replace(upper, upper == -Inf, 1)
            rel.age[,i] <- Lt.use$rel.age[upper] - 0.5    
            cohort[,i] <- Lt.use$ct[upper]
        }
    }else if(length(lfq$par$t_anchor) == 2){
        for(i in seq(length(lfq1$dates))){
            t.use <- which(lfq2$Lt$t == date2yeardec(lfq1$dates[i]))
            Lt.use <-  lfq2$Lt[t.use,]
            ## upper cohort bound
            upper <- suppressWarnings(apply(
                outer(X = lfq1$midLengths, Y = Lt.use$Lt, FUN = "-"),
                MARGIN = 1,
                FUN = function(x){max(which(x < 0))}
            ))
            ## If no upper exists, match to the oldest cohort (i.e. row 1 of Lt.use)
            upper <- replace(upper, upper == -Inf, 1)
            tmpa <- lfq2a$Lt[which(lfq2a$Lt$t == date2yeardec(lfq1$dates[i])),]
            tmpb <- lfq2b$Lt[which(lfq2b$Lt$t == date2yeardec(lfq1$dates[i])),]
            idx <- as.numeric(Lt.use$Lt %in% tmpa$Lt) + as.numeric(Lt.use$Lt %in% tmpb$Lt) * 2
            rel.age[,i] <- Lt.use$rel.age[upper] - (c(taDiff1,taDiff2)[idx])[upper]
            cohort[,i] <- Lt.use$ct[upper]
        }
    }

    lfq$cohort <- cohort
    lfq$rel.age <- rel.age

    return(lfq)
}


lfqCohortMT <- function(lfq){
  PAR <- lfq$par
  PAR2 <- PAR
  PAR2$t_anchor <- (PAR$t_anchor - 0.5) %% 1

  lfq <- lfqFitCurves(lfq, par = PAR, draw = FALSE)
  lfq2 <- lfqFitCurves(lfq, par = PAR2, agemax = lfq$agemax+0.5, draw = FALSE)

  rel.age <- lfq$catch*NaN
  cohort <- lfq$catch*NaN
  for(i in seq(length(lfq$dates))){
    t.use <- which(lfq2$Lt$t == date2yeardec(lfq$dates[i]))
    Lt.use <-  lfq2$Lt[t.use,]

    # upper cohort bound
    upper <- suppressWarnings(apply(
      outer(X = lfq$midLengths, Y = Lt.use$Lt, FUN = "-"),
      MARGIN = 1,
      FUN = function(x){max(which(x < 0))}
    ))

    # If no upper exists, match to the oldest cohort (i.e. row 1 of Lt.use)
    upper <- replace(upper, upper == -Inf, 1)

    rel.age[,i] <- Lt.use$rel.age[upper]-0.5
    cohort[,i] <- Lt.use$ct[upper]
  }

  lfq$cohort <- cohort
  lfq$rel.age <- rel.age

  return(lfq)

}
