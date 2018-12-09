#' Compute vertices of Voronoi diagram tiles.
#'
#' For use with polygon functions.
#' @param sites Object. Data frame of sites to construct Voronoi diagram with variables "x" and "y".
#' @param observed.data Object. Data frame of secondary source of data: observations, cases, addresses, etc. with variables "x" and "y".
#' @param rw Numeric. Vector of corners of the rectangular window or bounding box: xmin, xmax, ymin, ymax. For deldir::deldir().
#' @export
#' @examples
#' polygon.vertices <- deldirPolygons(cholera::pumps)
#' cholera::snowMap()
#' invisible(lapply(polygon.vertices, polygon))
#'
#' polygon.vertices <- deldirPolygons(cholera::pumps, cholera::roads)
#' cholera::snowMap()
#' invisible(lapply(polygon.vertices, polygon))

deldirPolygons <- function(sites, observed.data = NULL, rw = NULL) {
  if (is.null(observed.data) & is.null(rw)) {
    x.rng <- range(sites$x)
    y.rng <- range(sites$y)
  } else if (is.null(observed.data) == FALSE & is.null(rw)) {
    x.rng <- range(observed.data$x)
    y.rng <- range(observed.data$y)
  } else if (is.null(observed.data) & is.null(rw) == FALSE) {
    x.rng <- rw[1:2]
    y.rng <- rw[3:4]
  } else stop("Use either 'observed.data' or 'rw'; not both.")

  voronoi <- deldir::deldir(sites[, c("x", "y")], rw = c(x.rng, y.rng),
    suppressMsge = TRUE)
  cell.data <- deldir::tile.list(voronoi)
  lapply(cell.data, function(dat) data.frame(x = dat$x, y = dat$y))
}
