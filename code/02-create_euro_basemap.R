# Create a flat map of Europe
# Jonas Schöley
# 2020-03-05

library(tidyverse)
library(sf)
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)

eura_sf <-
  # download geospatial data for European, Asian and African countries
  ne_countries(continent = c('europe', 'asia', 'africa'),
               returnclass = 'sf', scale = 50) %>%
  # avoid self-intersection errors
  st_buffer(0) %>%
  # project to crs 3035
  st_transform(crs = 3035) %>%
  # merge into single polygon
  st_union(by_feature = FALSE) %>%
  st_crop(xmin = 25e5, xmax = 75e5, ymin = 13.5e5, ymax = 54.5e5)

# draw a basemap of Europe
euro_basemap <-
  ggplot(eura_sf) +
  geom_sf(
    aes(geometry = geometry),
    color = NA, fill = 'grey90'
  ) +
  coord_sf(expand = FALSE, datum = NA) +
  theme_void()

# save(
#   euro_basemap,
#   file = './data/euro_basemap.RData',
#   compress = 'xz'
# )
