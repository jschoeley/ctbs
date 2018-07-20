# Demonstration of the centered ternary balance scheme
# Jonas Schöley
# 2018-05-10

# This script reproduces the figures found in the paper "[]". The functions
# used below are simplified versions of those implemented in the R package
# "tricolore".

# Init --------------------------------------------------------------------

library(tidyverse)
library(ggtern)

# Data --------------------------------------------------------------------

# an empty, ready-to-use map of the european continent as ggplot object
load('data/euro_basemap.RData')
# polygon outlines of European NUTS-2 regions
load('data/euro_geo_nuts2.RData')
# European regional education composition
load('data/euro_education.RData')
# European regional labor force composition
load('data/euro_sectors.RData')

# Misc --------------------------------------------------------------------

# for plotting polygons with holes correctly
gghole <- function (fort) {
  poly <- fort[fort$id %in% fort[fort$hole, ]$id, ]
  hole <- fort[!fort$id %in% fort[fort$hole, ]$id, ]
  out <- list(poly, hole)
  names(out) <- c('poly', 'hole')
  return(out)
}

# Compositional Data Analysis ---------------------------------------------

#' Compositional Pertubation
#'
#' Pertube a compositional data set by a compositional vector.
#'
#' @param P n by m matrix of compositions {p1, ..., pm}_i for i=1,...,n.
#' @param c Compositional pertubation vector {c1, ..., cm}.
#'
#' @return n by m matrix of pertubed compositions.
#'
#' @examples
#' P <- prop.table(matrix(runif(12), 4), margin = 1)
#' cP <- Pertube(P, 1/Centre(P))
#' Centre(cP)
#'
#' @references
#' Von Eynatten, H., Pawlowsky-Glahn, V., & Egozcue, J. J. (2002).
#' Understanding perturbation on the simplex: A simple method to better
#' visualize and interpret compositional data in ternary diagrams.
#' Mathematical Geology, 34(3), 249-257.
#'
#' Pawlowsky-Glahn, V., Egozcue, J. J., & Tolosana-Delgado, R. (2007). Lecture
#' Notes on Compositional Data Analysis. Retrieved from
#' https://dugi-doc.udg.edu/bitstream/handle/10256/297/CoDa-book.pdf
Pertube <- function (P, c = rep(1/3, 3)) {
  return(prop.table(t(t(P)*c), margin = 1))
}

# Ternary Geometry --------------------------------------------------------

# T(K=k^2):   Equilateral triangle subdivided into K equilateral sub-triangles.
#             Each side of T is divided into k intervals of equal length.
# (p1,p2,p3): Barycentric coordinates wrt. T(K).
# id:         One-dimensional index of sub-triangles in T(K).
#
#                  p2           id index
#                  /\               9
#                 /  \            6 7 8
#                /____\         1 2 3 4 5
#              p1      p3

#' Centroid Coordinates of Sub-Triangles in Segmented Equilateral Triangle
#'
#' Segment an equilateral triangle into k^2 equilateral sub-triangles and return
#' the barycentric centroid coordinates of each sub-triangle.
#'
#' @param k Number of rows in the segmented equilateral triangle.
#'
#' @return A matrix of barycentric centroid coordinates of regions id=1,...,k^2.
#'
#' @references
#' S. H. Derakhshan and C. V. Deutsch (2009): A Color Scale for Ternary Mixtures.
#'
#' @examples
#' TernaryMeshCentroids(1)
#' TernaryMeshCentroids(2)
#' TernaryMeshCentroids(3)
TernaryMeshCentroids <- function (k) {
  # total number of centroids and centroid id
  K = k^2; id = 1:K

  # centroid coordinates as function of K and id
  g <- floor(sqrt(K-id)); gsq <- g^2
  c1 <- (((-K + id + g*(g+2) + 1) %% 2) - 3*gsq - 3*id + 3*K + 1) / (6*k)
  c2 <- -(((-K + gsq + id + 2*g + 1) %% 2) + 3*g - 3*k + 1) / (3*k)
  c3 <- (((-K + gsq + id + 2*g + 1) %% 2) + 3*gsq + 6*g + 3*id - 3*K + 1) / (6*k)

  return(cbind(id = id, p1 = c1, p2 = c2, p3 = c3))
}

#' Vertex Coordinates of Sub-Triangles in Segmented Equilateral Triangle
#'
#' Given the barycentric centroid coordinates of the sub-triangles in an
#' equilateral triangle subdivided into k^2 equilateral sub-triangles, return
#' the barycentric vertex coordinates of each sub-triangle.
#'
#' @param C n by 4 matrix of barycentric centroid coordinates of n=k^2
#'          sub-triangles. Column order: id, p1, p2, p3 with id=1,...,k^2.
#' @param k Number of rows in the segmented equilateral triangle.
#'
#' @return Index, vertex id and barycentric vertex coordinates for each of the
#'         k^2 sub-triangles.
#'
#' @examples
#' k = 2
#' C <- TernaryMeshCentroids(k)
#' TernaryMeshVertices(C)
#'
#' @references
#' S. H. Derakhshan and C. V. Deutsch (2009): A Color Scale for Ternary Mixtures.
TernaryMeshVertices <- function (C) {
  k <- sqrt(nrow(C))
  j <- k - floor(sqrt(k^2-C[,1]))
  i <- C[,1] - (j-1)*(2*k-j+1)
  term1 <- ((-1)^(i %% 2) * 2) / (3*k)
  term2 <- ((-1)^(i %% 2)) / (3*k)

  v1 <- cbind(C[,2] - term1, C[,3] + term2, C[,4] + term2)
  v2 <- cbind(C[,2] + term2, C[,3] - term1, C[,4] + term2)
  v3 <- cbind(C[,2] + term2, C[,3] + term2, C[,4] - term1)

  V <- cbind(C[,1], rep(1:3, each = nrow(C)), rbind(v1, v2, v3))
  colnames(V) <- c('id', 'vertex', 'p1', 'p2', 'p3')

  return(V)
}

#' Coordinates and Labels for the Centered Gridlines of a Ternary Diagram
TernaryCentroidGrid <- function (centroid) {

  # centroid percent difference labels
  labels = seq(-1, 1, 0.1)
  labels = data.frame(
    L = labels[labels >= -centroid[1]][1:10],
    T = labels[labels >= -centroid[2]][1:10],
    R = labels[labels >= -centroid[3]][1:10]
  )

  # breaks of uncentered grid
  breaks = data.frame(
    L = labels$L + centroid[1],
    T = labels$T + centroid[2],
    R = labels$R + centroid[3]
  )

  # grid L
  gridL =
    data.frame(
      scale = 'L',
      centroid = ifelse(breaks$L == centroid[1], TRUE, FALSE),
      L_from = breaks$L,
      T_from = 1-breaks$L,
      R_from = 0,
      L_to = breaks$L,
      T_to = 0,
      R_to = 1-breaks$L
    )

  # grid T
  gridT =
    data.frame(
      scale = 'T',
      centroid = ifelse(breaks$T == centroid[2], TRUE, FALSE),
      L_from = 0,
      T_from = breaks$T,
      R_from = 1-breaks$T,
      L_to = 1-breaks$T,
      T_to = breaks$T,
      R_to = 0
    )

  # grid R
  gridR =
    data.frame(
      scale = 'R',
      centroid = ifelse(breaks$R == centroid[3], TRUE, FALSE),
      L_from = 1-breaks$R,
      T_from = 0,
      R_from = breaks$R,
      L_to = 0,
      T_to = 1-breaks$R,
      R_to = breaks$R
    )

  # grid line coordinates of uncentered grid
  grid = rbind(gridL, gridT, gridR)

  # grid line coordinates of centered grid
  cgrid = data.frame(
    grid[,1:2],
    prop.table(t(t(grid[,3:5])*(1/centroid)), margin = 1),
    prop.table(t(t(grid[,6:8])*(1/centroid)), margin = 1)
  )

  # breaks of centered grid
  cbreaks = data.frame(L = cgrid[cgrid$scale == 'L', 'L_from'],
                       T = cgrid[cgrid$scale == 'T', 'T_from'],
                       R = cgrid[cgrid$scale == 'R', 'R_from'])

  list(grid = grid, cgrid = cgrid,
       breaks = breaks, cbreaks = cbreaks, labels = labels)

}

# Ternary Color Scale -----------------------------------------------------

#' RGB Mixture of Ternary Composition
#'
#' Return the ternary balance scheme colors for a matrix of ternary compositions.
#'
#' @param P n by 3 matrix of ternary compositions {p1, p2, p3}_i for
#'          i=1, ..., n.
#' @param breaks Number of breaks in the discrete color scale. An integer >0.
#'               Values above 99 imply no discretization.
#' @param h_ Primary hue of the first ternary element in angular degrees [0, 360].
#' @param c_ Maximum possible chroma of mixed colors [0, 200].
#' @param l_ Lightness of mixed colors [0, 100].
#' @param contrast Lightness contrast of the color scale [0, 1).
#' @param center Ternary coordinates of the grey-point.
#' @param spread Spread of the color scale around center > 0.
#'
#' @return An n row data frame giving, for each row of the input P, the input
#' proportions (p1, p2, p3), parameters of the color mixture (h, c, l) and the
#' hexsrgb string of the mixed colors.
#'
#' @examples
#' P <- prop.table(matrix(runif(9), ncol = 3), 1)
#' ColorMap(P, breaks = 5, h_ = 80, c_ = 170, l_ = 80,
#'          contrast = 0.6, center = rep(1/3, 3), spread = 1)
ColorMap <- function (P, h_, c_, l_, contrast, center) {

  # generate primary colors starting with a hue value in [0, 360) and then
  # picking two equidistant points on the circumference of the color wheel.
  # input hue in degrees, all further calculations in radians.
  phi <- (h_*0.0174 + c(0, 2.09, 4.19)) %% 6.28

  # closing
  P <- P_raw <- prop.table(P, margin = 1)
  # centering
  P <- Pertube(P, 1/center)

  # calculate the chroma matrix C by scaling the row proportions
  # of the input matrix P by the maximum chroma parameter.
  C <- P*c_

  # the complex matrix Z represents each case (i) and group (j=1,2,3) specific
  # color in complex polar form with hue as angle and chroma as radius.
  Z <- matrix(complex(argument = phi, modulus = c(t(C))), ncol = 3, byrow = TRUE)

  # adding up the rows gives the CIE-Lab (cartesian) coordinates
  # of the convex color mixture in complex form.
  z <- rowSums(Z)
  # convert the cartesian CIE-Lab coordinates to polar CIE-Luv coordinates
  # and add lightness level.
  M <- cbind(h = (Arg(z)*57.3)%%360, c = Mod(z), l = l_)

  # decrease lightness and chroma towards the center of the color scale
  cfactor <- M[,2]*contrast/c_ + 1-contrast
  M[,3] <- cfactor*M[,3]
  M[,2] <- cfactor*M[,2]

  # convert the complex representation of the color mixture to
  # hex-srgb representation via the hcl (CIE-Luv) color space
  hexsrgb <- hcl(h = M[,1], c = M[,2], l = M[,3],
                 alpha = 1, fixup = TRUE)

  # original compositions, pertubed compositions, hcl values of mixtures and hexsrgb code
  result <- data.frame(P_raw, P, M[,1], M[,2], M[,3], hexsrgb,
                       row.names = NULL, check.rows = FALSE,
                       check.names = FALSE, stringsAsFactors = FALSE)
  colnames(result) <- c('p1', 'p2', 'p3', 'cp1', 'cp2', 'cp3', 'h', 'c', 'l', 'hexsrgb')
  return(result)
}

#' Ternary Balance Scheme Legend
#'
#' Plot a ternary balance scheme legend.
#'
#' @inheritParams ColorMap
#'
#' @examples
#' ColorKey(breaks = 5, h_ = 0, c_ = 140, l_ = 70, contrast = 0.5,
#'          center = rep(1/3, 3), spread = 1)
ColorKey <- function (h_, c_, l_, contrast, center) {

  # partition the ternary legend into k^2 equilateral sub-triangles
  # calculate ternary vertex coordinates and fill color for each sub-triangle.
  k = 100
  C <- TernaryMeshCentroids(k)
  V <- TernaryMeshVertices(C)
  rgbs <- ColorMap(P = C[,-1], h_, c_, l_, contrast, center)[['hexsrgb']]
  legend_surface <- data.frame(V, rgb = rep(rgbs, 3))

  # plot the legend
  legend <-
    # basic legend
    ggtern(legend_surface, aes_string(x = 'p1', y = 'p2', z = 'p3')) +
    geom_polygon(aes_string(group = 'id', fill = 'rgb', color = 'rgb'), lwd = 1) +
    geom_mask() +
    # rgb color input
    scale_color_identity(guide = FALSE) +
    scale_fill_identity(guide = FALSE) +
    # theme
    theme_classic() +
    theme(tern.axis.title.L = element_text(hjust = 0.2, vjust = 1, angle = -60),
          tern.axis.title.R = element_text(hjust = 0.8, vjust = 0.6, angle = 60))

  return(legend)
}

# Ternary diagram of regional education levels in Europe ------------------

# perform color coding
euro_education_colors <-
  ColorMap(as.matrix(euro_education[,2:4]),
           h_ = 126, c_ = 200, l_ = 90, contrast = 0.7, center = c(1/3, 1/3, 1/3))
euro_education_colors$id <- euro_education$id

# generate corresponding color key
plot_euro_education_key <-
  ColorKey(h_ = 126, c_ = 200, l_ = 90, contrast = 0.7, center = c(1/3, 1/3, 1/3)) +
  geom_point(aes(x = p1, y = p2, z = p3), shape = 21, data = euro_education_colors) +
  lline(Lintercept = seq(0.1, 0.9, 0.1), alpha = 0.2) +
  tline(Tintercept = seq(0.1, 0.9, 0.1), alpha = 0.2) +
  rline(Rintercept = seq(0.1, 0.9, 0.1), alpha = 0.2) +
  labs(L = '% Lower secondary\nor less', T = '% Upper secondary', R = '% Tertiary',
       title = 'A',
       caption = 'Ternary diagram displaying the population composition by education level for European regions.') +
  theme_arrowlarge()

ggsave(filename = 'euro_education_key.pdf', path = 'out/', plot = plot_euro_education_key,
       width = 6, height = 6)

# Ternary balance scheme map of regional education levels in Europe -------

euro_education_map <- gghole(left_join(euro_education_colors, euro_geo_nuts2, by = 'id'))

plot_euro_education_map <-
  euro_basemap +
  geom_polygon(aes(x = long, y = lat, group = group, fill = hexsrgb),
               data = euro_education_map$poly) +
  geom_polygon(aes(x = long, y = lat, group = group, fill = hexsrgb),
               data = euro_education_map$hole) +
  scale_fill_identity() +
  labs(title = 'B', caption = 'Ternary balance scheme map displaying the population composition by education level for European regions.')

ggsave(filename = 'euro_education_map.pdf', path = 'out/', plot = plot_euro_education_map,
       width = 6, height = 6)

# TBS map of regional labor force distribution in Europe ------------------

# perform color coding
euro_sectors_colors <-
  ColorMap(as.matrix(euro_sectors[,2:4]),
           h_ = 126, c_ = 200, l_ = 90, contrast = 0.7, center = c(1/3, 1/3, 1/3))
euro_sectors_colors$id <- euro_sectors$id

# generate corresponding color key
plot_euro_sectors_key <-
  ColorKey(h_ = 126, c_ = 200, l_ = 90, contrast = 0.7, center = c(1/3, 1/3, 1/3)) +
  geom_point(aes(x = p1, y = p2, z = p3),
             alpha = 0.2, size = 1,
             data = euro_sectors_colors) +
  labs(L = '% Primary', T = '% Secondary', R = '% Tertiary')

# join geographic data with color information
euro_sectors_map <- gghole(left_join(euro_sectors_colors, euro_geo_nuts2, by = 'id'))

# generate choropleth map
plot_euro_sectors_map <-
  euro_basemap +
  geom_polygon(aes(x = long, y = lat, group = group, fill = hexsrgb),
               data = euro_sectors_map$poly) +
  geom_polygon(aes(x = long, y = lat, group = group, fill = hexsrgb),
               data = euro_sectors_map$hole) +
  scale_fill_identity() +
    annotation_custom(
    ggplotGrob(
      plot_euro_sectors_key +
        theme(plot.background = element_rect(fill = NA, color = NA))
    ),
    xmin = 53e5, xmax = Inf, ymin = 35e5, ymax = Inf) +
  labs(title = 'C', caption = 'Ternary balance scheme map displaying the workforce composition by economic sector for European regions.')

ggsave(filename = 'euro_sectors_map.pdf', path = 'out/', plot = plot_euro_sectors_map,
       width = 6, height = 6)

# Centered TBS map of regional labor force distribution in Europe ---------

# perform color coding
euro_sectors_colors_centered <-
  ColorMap(as.matrix(euro_sectors[,2:4]),
           h_ = 126, c_ = 200, l_ = 90, contrast = 0.7,
           center = c(0.04, 0.24, 0.72))
euro_sectors_colors_centered$id <- euro_sectors$id

# generate corresponding color key
plot_euro_sectors_key_centered <-
  ColorKey(h_ = 126, c_ = 200, l_ = 90, contrast = 0.7, center = c(0.04, 0.24, 0.72)) +
  geom_point(aes(x = p1, y = p2, z = p3),
             alpha = 0.2, size = 1,
             data = euro_sectors_colors) +
  labs(L = '% Primary', T = '% Secondary', R = '% Tertiary')

# join geographic data with color information
euro_sectors_map_centered <- gghole(left_join(euro_sectors_colors_centered, euro_geo_nuts2, by = 'id'))

# generate choropleth map
plot_euro_sectors_map_centered <-
  euro_basemap +
  geom_polygon(aes(x = long, y = lat, group = group, fill = hexsrgb),
               data = euro_sectors_map_centered$poly) +
  geom_polygon(aes(x = long, y = lat, group = group, fill = hexsrgb),
               data = euro_sectors_map_centered$hole) +
  scale_fill_identity() +
    annotation_custom(
    ggplotGrob(
      plot_euro_sectors_key_centered +
        theme(plot.background = element_rect(fill = NA, color = NA))
    ),
    xmin = 53e5, xmax = Inf, ymin = 35e5, ymax = Inf) +
  labs(title = 'D', caption = 'Centered ternary balance scheme map displaying the regional deviation in workforce composition from the European average.')

ggsave(filename = 'euro_sectors_map_centered.pdf', path = 'out/', plot = plot_euro_sectors_map_centered,
       width = 6, height = 6)

# Non centered legend -----------------------------------------------------

legend_style_0 <-
  ColorKey(h_ = 126, c_ = 200, l_ = 90, contrast = 0.7,
           center = c(1/3, 1/3, 1/3)) +
  # grid
  lline(Lintercept = seq(0.1, 0.9, 0.1), alpha = 0.2) +
  tline(Tintercept = seq(0.1, 0.9, 0.1), alpha = 0.2) +
  rline(Rintercept = seq(0.1, 0.9, 0.1), alpha = 0.2) +
  # center
  lline(Lintercept = 0.04) +
  tline(Tintercept = 0.24) +
  rline(Rintercept = 0.72) +
  geom_point(aes(x = p1, y = p2, z = p3),
             alpha = 0.2, size = 1,
             data = euro_sectors_colors) +
  scale_alpha_manual(values = c(`FALSE` = 0.2, `TRUE` = 1)) +
  labs(x = '% primary', y = '% secondary', z = '% ternary')

ggsave(filename = 'legend_style0.pdf', path = 'out/', plot = legend_style_0,
       width = 6, height = 6)

# Transformed data on standard grid ---------------------------------------

legend_style_1 <-
  ColorKey(h_ = 126, c_ = 200, l_ = 90, contrast = 0.7,
           center = c(1/3, 1/3, 1/3)) +
  lline(Lintercept = seq(0.1, 0.9, 0.1), alpha = 0.2) +
  tline(Tintercept = seq(0.1, 0.9, 0.1), alpha = 0.2) +
  rline(Rintercept = seq(0.1, 0.9, 0.1), alpha = 0.2) +
  geom_point(aes(x = cp1, y = cp2, z = cp3),
             alpha = 0.2, size = 1,
             data = euro_sectors_colors_centered) +
  labs(x = '% green', y = '% blue', z = '% red')

ggsave(filename = 'legend_style1.pdf', path = 'out/', plot = legend_style_1,
       width = 6, height = 6)

# Transformed data on transformed gridlines  -----------------

legend_style_2 <-
  ColorKey(h_ = 126, c_ = 200, l_ = 90, contrast = 0.7,
           center = c(1/3, 1/3, 1/3)) +
  geom_segment(aes(x = L_from, xend = L_to,
                   y = T_from, yend = T_to,
                   z = R_from, zend = R_to, alpha = centroid),
               show.legend = FALSE, data = grids$cgrid) +
  geom_point(aes(x = cp1, y = cp2, z = cp3),
             alpha = 0.2, size = 1,
             data = euro_sectors_colors_centered) +
  scale_L_continuous(breaks = grids$cbreaks$L,
                     labels = paste0(round(grids$labels$L, 2), ' (', round(grids$breaks$L, 2), ')')) +
  scale_T_continuous(breaks = grids$cbreaks$T,
                     labels = paste0(round(grids$labels$T, 2), ' (', round(grids$breaks$T, 2), ')')) +
  scale_R_continuous(breaks = grids$cbreaks$R,
                     labels = paste0(round(grids$labels$R, 2), ' (', round(grids$breaks$R, 2), ')')) +
  scale_alpha_manual(values = c(`FALSE` = 0.2, `TRUE` = 1)) +
  labs(x = '%pt. diff. primary (%)', y = '%pt. diff. secondary (%)', z = '%pt. diff. ternary (%)')


ggsave(filename = 'legend_style2.pdf', path = 'out/', plot = legend_style_2,
       width = 6, height = 6)


# Transformed color key ---------------------------------------------------

legend_style_3 <-
  ColorKey(h_ = 126, c_ = 200, l_ = 90, contrast = 0.7,
           center = c(0.04, 0.24, 0.72)) +
  # grid
  lline(Lintercept = seq(0.1, 0.9, 0.1), alpha = 0.2) +
  tline(Tintercept = seq(0.1, 0.9, 0.1), alpha = 0.2) +
  rline(Rintercept = seq(0.1, 0.9, 0.1), alpha = 0.2) +
  # center
  lline(Lintercept = 0.04) +
  tline(Tintercept = 0.24) +
  rline(Rintercept = 0.72) +
  geom_point(aes(x = cp1, y = cp2, z = cp3),
             alpha = 0.2, size = 1,
             data = euro_sectors_colors) +
  labs(x = '% primary', y = '% secondary', z = '% ternary')


ggsave(filename = 'legend_style3.pdf', path = 'out/', plot = legend_style_3,
       width = 6, height = 6)


