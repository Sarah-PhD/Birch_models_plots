=========================================================
  # Sweden map: ONE point per area (equal size) + non-overlapping labels
  # Requires: area_tbl_out with columns: area, lon, lat (WGS84)
  #=========================================================

library(dplyr)
library(sf)
library(ggplot2)
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthdata)

# Sweden basemap
sweden <- ne_countries(scale = "medium", country = "Sweden", returnclass = "sf")

# Make labels readable (does NOT change grouping)
pretty_area_label <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[._]+", " ") %>%  # dots/underscores -> spaces
    stringr::str_squish()
}

# One point per area (equal size)
areas_sf <- area_tbl_out %>%
  mutate(area_label = pretty_area_label(area)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# We’ll use coordinates explicitly for ggrepel
areas_xy <- areas_sf %>%
  st_drop_geometry() %>%
  bind_cols(as.data.frame(st_coordinates(areas_sf))) %>%
  rename(lon = X, lat = Y)

p_map_labels <- ggplot() +
  geom_sf(data = sweden, fill = "grey95", colour = "grey60") +
  
  # equal-size points
  geom_point(
    data = areas_xy,
    aes(x = lon, y = lat),
    size = 2.8,
    colour = "#1f77b4",
    alpha = 0.95
  ) +
  
  # repelled labels (no overlap as far as possible)
  geom_text_repel(
    data = areas_xy,
    aes(x = lon, y = lat, label = area_label),
    size = 3.2,
    box.padding = 0.6,
    point.padding = 0.35,
    min.segment.length = 0,
    segment.alpha = 0.35,
    max.overlaps = Inf,
    seed = 1
  ) +
  
  theme_classic(base_size = 13) +
  labs(x = NULL, y = NULL, title = "Study areas across Sweden")

p_map_labels

# Optional (often helps A LOT): zoom to Sweden bounds so labels have room
# p_map_labels + coord_sf(xlim = c(10, 25), ylim = c(55, 70))

## Table
library(dplyr)

area_table <- Birch_2425 %>%
  mutate(
    north = suppressWarnings(as.numeric(north)),
    east  = suppressWarnings(as.numeric(east))
  ) %>%
  filter(!is.na(area), !is.na(north), !is.na(east)) %>%
  group_by(area) %>%
  summarise(
    mean_east  = round(mean(east,  na.rm = TRUE), 1),
    mean_north = round(mean(north, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(area)

area_table
print(area_table, n = Inf)
