#' Predict gillnet selectivity (old Millar method)
#'
#' @param type ("norm.loc", "norm.sca", "gamma", "lognorm")
#' @param meshsizes mesh sizes
#' @param rel relative powers (one for each mesh size)
#' @param pars selection curve parameters
#' @param plotlens vector of new mesh sizes for selectivity prediction
#'
#' @return selectivities
#' @export
#'
rcurves=function(type,meshsizes,rel,pars,plotlens) {
  lens=rep(plotlens,length(meshsizes))
  relsizes=meshsizes/meshsizes[1]
  if(is.null(rel)) releff=1
   else releff=rep(rel,rep(length(plotlens),length(meshsizes)))
  msizes=rep(meshsizes,rep(length(plotlens),length(meshsizes)))
  relsizes=rep(relsizes,rep(length(plotlens),length(relsizes)))
  switch(type,
   "norm.loc"={ k=pars[1]; mean1=pars[3]; sigma=pars[2]
     rselect=exp(-((lens-k*msizes)^2)/(2*sigma^2)) },
   "norm.sca"={ k1=pars[1]; k2=pars[2]
     rselect=exp(-((lens-k1*msizes)^2)/(2*k2*msizes^2)) },
   "gamma"={ alpha=pars[1]; k=pars[2]
     rselect=(lens/((alpha-1)*k*msizes))^(alpha-1)*exp(alpha-1-lens/(k*msizes)) },
   "lognorm"={ mu1=pars[1]; sigma=pars[2]
     rselect=
      (1/lens)*exp(mu1+log(relsizes)-sigma^2/2-
                    (log(lens)-mu1-log(relsizes))^2/(2*sigma^2)) })
  rselect=releff*rselect/max(releff)
  rselect=matrix(rselect,ncol=length(meshsizes))
  return(rselect)
  }
