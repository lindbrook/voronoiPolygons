
<!-- README.md is generated from README.Rmd. Please edit that file -->
VoronoiPolygons: transform 'deldir' Voronoi tiles into polygons
---------------------------------------------------------------

Using the locations of sites or landmarks of interest, Voronoi tessellation partitions a space into cells or tiles that represent neighborhoods (i.e., catchment or service areas). Then, by computing the vertices of those tiles, we can leverage various functions that use polygons and do things like color code neighborhoods or count elements of interest within neighborhoods.

As an example, I use data from John Snow's map of the 1854 cholera outbreak in the Soho area London.

Coloring Neighborhoods
----------------------

To color code neighborhood tiles, we can use graphics::polygon().

``` r
polygon.vertices <- VoronoiPolygons(cholera::pumps, cases = cholera::roads,
  plot.frame = "cases")

snow.colors <- grDevices::adjustcolor(cholera::snowColors(), alpha.f = 1/3)

cholera::snowMap(add.cases = FALSE)
cholera::addNeighborhoodCases(metric = "euclidean")

invisible(lapply(seq_along(polygon.vertices), function(i) {
  polygon(polygon.vertices[[i]], col = snow.colors[[i]])
}))
```

<img src="README_files/figure-markdown_github/coloring-1.png" style="display: block; margin: auto;" />

Counting Cases
--------------

To count the number of cases within each neighborhood, we can use sp::point.in.polygon().

``` r
polygon.vertices <- VoronoiPolygons(cholera::pumps, cases = cholera::roads,
  plot.frame = "cases")

cases <- cholera::fatalities.unstacked

census <- lapply(polygon.vertices, function(tile) {
  sp::point.in.polygon(cases$x, cases$y, tile$x, tile$y)
})

neighborhood.pump <- paste0("p", cholera::pumps$id)

stats::setNames(vapply(census, sum, integer(1L)), neighborhood.pump)
#>  p1  p2  p3  p4  p5  p6  p7  p8  p9 p10 p11 p12 p13 
#>   0   1  13  23   6  61 361  16  27  62   2   2   4
```

### getting started

To install the development version of 'VoronoiPolygons' from GitHub:

``` r
# Note that you may need to install the 'devtools' package:
# install.packages("devtools")

# For 'devtools' (< 2.0.0)
devtools::install_github("lindbrook/VoronoiPolygons", build_vignettes = TRUE)

# For 'devtools' (>= 2.0.0)
devtools::install_github("lindbrook/VoronoiPolygons"))
```
