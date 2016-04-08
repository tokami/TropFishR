#' @name trawl_fishery_Java
#'
#' @title Data from the trawl fishery off the North coast of Java
#'
#' @description Times series of catch and effort data from the trawl fishery off the
#'  North coast of Java. This dataset can
#'  be used for the estimation of maximum sustainable yield by means of the
#'  production models (\code{\link{prod_mod}} and \code{\link{prod_mod_ts}}).
#'
#' @docType data
#' @format A dataframe consisting of: 1. \strong{year} a vector with years, 2. \strong{Y}
#'  yield [1000 tons], and 3. \strong{f} fishing effort [number of standard vessels].
#'
#' @source Dwiponggo, A., 1979. Review of the demersal resources and fisheries in the
#'  Java Sea. IPFC:RRD/II/79/Inf.12. Paper presented at the SCORRAD Meeting, 1979,
#'  Hong Kong.
#'
#' @usage data(trawl_fishery_Java)
#'
#' @keywords data dataset trawl Java
#'
#' @examples
#' data(trawl_fishery_Java)
#' str(trawl_fishery_Java)
#' summary(trawl_fishery_Java)
#'
NULL
