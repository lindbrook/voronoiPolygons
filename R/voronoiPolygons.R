#' Compute vertices of Voronoi diagram tiles.
#'
#' For use with polygon functions.
#' @param sites Object. Data frame of sites to construct Voronoi diagram with variables "x" and "y".
#' @param cases Object. Data frame of secondary source of data: observations, cases, addresses, etc. with variables "x" and "y".
#' @param plot.frame Character. "cases" or "sites". Data to define extent or range of Voronoi diagram.
#' @export
#' @examples
#' pump.neighborhoods <- VoronoiPolygons(cholera::pumps, cases)
#' cholera::snowMap()
#' invisible(lapply(pump.neighborhoods, polygon))
#'
#' pump.neighborhoods <- VoronoiPolygons(cholera::pumps, cases = cholera::roads, plot.frame = "cases")
#' cholera::snowMap()
#' invisible(lapply(pump.neighborhoods, polygon))

VoronoiPolygons <- function(sites, cases = NULL, plot.frame = "sites") {
  if (plot.frame == "sites") {
    x.rng <- range(sites$x)
    y.rng <- range(sites$y)
  } else if (plot.frame == "cases") {
    x.rng <- range(cases$x)
    y.rng <- range(cases$y)
  } else stop('plot.frame must either be "cases" or "sites".')

  four.corners <- list(nw = data.frame(x = min(x.rng), y = max(y.rng)),
                       ne = data.frame(x = max(x.rng), y = max(y.rng)),
                       se = data.frame(x = max(x.rng), y = min(y.rng)),
                       sw = data.frame(x = min(x.rng), y = min(y.rng)))

  voronoi <- deldir::deldir(sites[, c("x", "y")], rw = c(x.rng, y.rng),
    suppressMsge = TRUE)

  cell.data <- voronoi$dirsgs

  coordinates <- lapply(seq_len(nrow(sites)), function(i) {
    dat <- sites[i, c("x", "y")]
    cell <- cell.data[cell.data$ind1 == i | cell.data$ind2 == i, ]
    a <- cell[, c("x1", "y1")]
    b <- cell[, c("x2", "y2")]
    names(a) <- c("x", "y")
    names(b) <- c("x", "y")
    coords <- unique(rbind(a, b))

    # test for "open" polygons
    test1 <- any(cell$thirdv1 < 0 | cell$thirdv2 < 0)
    test2 <- unlist(cell[, c("thirdv1", "thirdv2")])
    test2 <- length(unique(test2[test2 < 0])) != 1

    # include corner vertex
    if (test1 & test2) {
      corner.vertex <- vapply(four.corners, function(corner) {
        vertical <- signif(corner$x) == signif(coords$x)
        horizontal <- signif(corner$y) == signif(coords$y)
        any(vertical) & any(horizontal)
      }, logical(1L))

      coords <- rbind(four.corners[[names(which(corner.vertex))]], coords)
    }

    # center vertices relative to landmark's coordinates
    coords.centered <- data.frame(x = coords$x - dat$x, y = coords$y - dat$y)

    # transform coordinates from cartesian to polar
    # sort vertices by phi (angle)
    # returns vertices in counter-clockwise order
    idx <- order(apply(coords.centered, 1, pracma::cart2pol)[1, ])
    coords <- coords[idx, ]

    # adds first vertex to last to close polygon
    rbind(coords, coords[1, ])
  })

  coordinates
}
