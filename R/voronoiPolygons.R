#' Compute vertices of Voronoi diagram tiles.
#'
#' For use with polygon functions.
#' @param sites Object. Data frame of sites to construct Voronoi diagram with variables "x" and "y".
#' @param observed.data Object. Data frame of secondary source of data: observations, cases, addresses, etc. with variables "x" and "y".
#' @param rw Numeric. Vector of corners of the rectangular window or bounding box: xmin, xmax, ymin, ymax. For deldir::deldir().
#' @export
#' @examples
#' polygon.vertices <- VoronoiPolygons(cholera::pumps)
#' cholera::snowMap()
#' invisible(lapply(polygon.vertices, polygon))
#'
#' polygon.vertices <- VoronoiPolygons(cholera::pumps, cholera::roads)
#' cholera::snowMap()
#' invisible(lapply(polygon.vertices, polygon))

VoronoiPolygons <- function(sites, observed.data = NULL, rw = NULL) {
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

    # append first vertex to last observation to close polygon
    rbind(coords, coords[1, ])
  })

  coordinates
}

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

  four.corners <- list(nw = data.frame(x = min(x.rng), y = max(y.rng)),
                       ne = data.frame(x = max(x.rng), y = max(y.rng)),
                       se = data.frame(x = max(x.rng), y = min(y.rng)),
                       sw = data.frame(x = min(x.rng), y = min(y.rng)))

  voronoi <- deldir::deldir(sites[, c("x", "y")], rw = c(x.rng, y.rng),
    suppressMsge = TRUE)

  cell.data <- deldir::tile.list(voronoi)

  lapply(cell.data, function(dat) data.frame(x = dat$x, y = dat$y))
}
