#' Compute vertices of Voronoi diagram tiles (beta).
#'
#' For use with polygon functions.
#' @param sites Object. Data frame of sites to construct Voronoi diagram with variables "x" and "y".
#' @param rw.data Object. Data frame of secondary source of data for rectangular window or bounding box: observations, cases, addresses, etc. with variables "x" and "y".
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

VoronoiPolygons <- function(sites, rw.data = NULL, rw = NULL) {
  if (is.null(rw.data) & is.null(rw)) {
    x.rng <- range(sites$x)
    y.rng <- range(sites$y)
  } else if (is.null(rw.data) == FALSE & is.null(rw)) {
    x.rng <- range(rw.data$x)
    y.rng <- range(rw.data$y)
  } else if (is.null(rw.data) & is.null(rw) == FALSE) {
    x.rng <- rw[1:2]
    y.rng <- rw[3:4]
  } else stop("Use either 'rw.data' or 'rw'; not both.")

  four.corners <- list(nw = data.frame(x = min(x.rng), y = max(y.rng)),
                       ne = data.frame(x = max(x.rng), y = max(y.rng)),
                       se = data.frame(x = max(x.rng), y = min(y.rng)),
                       sw = data.frame(x = min(x.rng), y = min(y.rng)))

  voronoi <- deldir::deldir(sites[, c("x", "y")], rw = c(x.rng, y.rng),
    suppressMsge = TRUE)

  cell.data <- voronoi$dirsgs

  lapply(seq_len(nrow(sites)), function(i) {
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

    # corner vertex check
    if (test1 & test2) {
      distance.to.corner <- vapply(four.corners, function(corner) {
        stats::dist(rbind(dat, corner))
      }, numeric(1L))

      EW.test <- signif(coords$x) %in% signif(x.rng)
      NS.test <- signif(coords$y) %in% signif(y.rng)

      if (sum(NS.test) == 2 | sum(EW.test) == 2) {
        nms <- names(sort(distance.to.corner))[1:2]
        coords <- rbind(coords, do.call(rbind, four.corners[nms]))
      } else if (nrow(coords) == 2) {

        segment.test <- lapply(four.corners, function(corner) {
          test <- segmentIntersection(dat$x, dat$y, corner$x, corner$y,
          coords[1, "x"], coords[1, "y"], coords[2, "x"], coords[2, "y"])
          all(is.na(test))
        })

        sel <- names(segment.test)[unlist(segment.test)]
        coords <- rbind(coords, do.call(rbind, four.corners[sel]))

      } else {
        coords <- rbind(coords, four.corners[[which.min(distance.to.corner)]])
      }
    }

    # center vertices relative to landmark's coordinates
    coords.centered <- data.frame(x = coords$x - dat$x, y = coords$y - dat$y)

    # vertices in counter-clockwise order
    idx <- order(apply(coords.centered, 1, pracma::cart2pol)[1, ])
    coords[idx, ]
  })
}

segmentIntersection <- function(x1, y1, x2, y2, a1, b1, a2, b2) {
  # returns the point of intersection between two segments or NA if none.
  # http://stackoverflow.com/questions/20519431/finding-point-of-intersection-in-r
  # x1, y1, x2, y2 coordinates of first segment's endpoints.
  # a1, b1, a2, b2 coordinates of second segment's endpoints.
  denom <- (b2 - b1) * (x2 - x1) - (a2 - a1) * (y2 - y1)
  denom[abs(denom) < 1e-10] <- NA # parallel lines
  ua <- ((a2 - a1) * (y1 - b1) - (b2 - b1) * (x1 - a1)) / denom
  ub <- ((x2 - x1) * (y1 - b1) - (y2 - y1) * (x1 - a1)) / denom
  x <- x1 + ua * (x2 - x1)
  y <- y1 + ua * (y2 - y1)
  inside <- (ua >= 0) & (ua <= 1) & (ub >= 0) & (ub <= 1)
  data.frame(x = ifelse(inside, x, NA), y = ifelse(inside, y, NA))
}
