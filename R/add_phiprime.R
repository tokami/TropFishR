#' @title Add phi prime contours to a Linf-K plot
#'
#' @description Adds isolines of the growth performance index
#'   \eqn{\phi' = \log_{10}(K) + 2 \log_{10}(L_\infty)} to an existing plot
#'   with Linf on the x-axis and K on the y-axis (e.g. \link{LinfK_scatterhist}).
#'
#' @param gridsize default = 20
#' @param ... Additional arguments passed to \code{\link[graphics]{contour}}.
#'
#' @return plot
#'
#' @importFrom graphics par contour
#'
#' @examples
#' plot(c(5, 20), c(0.1, 2), t = "n", xlab = "Linf", ylab = "K")
#' add_phiprime()
#'
#' @export
#'
add_phiprime <- function(gridsize = 20, ...){
    usr <- par()$usr
    Linf = seq(usr[1], usr[2], length.out = gridsize)
    K = seq(usr[3], usr[4], length.out = gridsize)
    grd <- expand.grid(
        Linf = Linf,
        K = K
    )
    grd$phiL <- log10(grd$K) + 2 * log10(grd$Linf)

    M <- list(x = Linf, y = K, z = matrix(grd$phiL, nrow = gridsize, ncol = gridsize))
    contour(x = M, add = TRUE, ...)
}
