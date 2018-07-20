#' ---
#' title: European regional labor force composition 2016
#' author: Jonas Schöley
#' output:
#'   github_document
#' ---

#+echo=FALSE
knitr::opts_chunk$set(warning=FALSE, message=FALSE,
                      fig.width = 12, fig.height = 12)

devtools::install_github('jschoeley/tricolore')
library(tricolore); DemoTricolore()

euro_sectors <- subset(euro_sectors, year == 2016)

centre = tricolore:::Centre(as.matrix(euro_sectors[,3:5]))

# Raw proportions ---------------------------------------------------------

# generate colors based on compositions in `euro_sectors`, default options
tricol <- Tricolore(euro_sectors, 'primary', 'secondary', 'tertiary',
                    show_center = FALSE, hue = 0.35, k = 100)

# merge vector of colors with with map data
euro_sectors$srgb <- tricol$hexsrgb
map_data <- dplyr::left_join(euro_geo_nuts2, euro_sectors, by = c('id' = 'nuts2'))

p1 <-
  euro_basemap +
  geom_polygon(aes(long, lat, group = group, fill = srgb, color = srgb),
               data = map_data) +
  scale_fill_identity() +
  scale_color_identity() +
  annotation_custom(
    ggplotGrob(
      tricol$legend +
        theme(plot.background = element_rect(fill = NA, color = NA))
    ),
    xmin = 53e5, xmax = Inf, ymin = 35e5, ymax = Inf)

# Centered proportions ----------------------------------------------------

# generate colors based on compositions in `euro_sectors`, default options
tricol <- Tricolore(euro_sectors, 'primary', 'secondary', 'tertiary',
                    center = NA, show_center = TRUE, hue = 0.35)

# merge vector of colors with with map data
euro_sectors$srgb <- tricol$hexsrgb
map_data <- dplyr::left_join(euro_geo_nuts2, euro_sectors, by = c('id' = 'nuts2'))

p2 <-
  euro_basemap +
  geom_polygon(aes(long, lat, group = group, fill = srgb, color = srgb),
               data = map_data) +
  scale_fill_identity() +
  scale_color_identity() +
  annotation_custom(
    ggplotGrob(
      tricol$legend +
        theme(plot.background = element_rect(fill = NA, color = NA))
    ),
    xmin = 53e5, xmax = Inf, ymin = 35e5, ymax = Inf)

# Merge -------------------------------------------------------------------

grid.arrange(p1, p2, ncol = 2)
