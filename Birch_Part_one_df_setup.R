#-------------------------------------------------
# Title: Birch analysis 2.0 from the beginning (ÄBIN_compact version)
# Author: Sarah Gore
# Date: 2025-02-03
# Description: Same workflow as your original script, but using ÄBIN_compact
#              (new column titles) + using 2024 & 2025 + fixes:
#              (1) swap north/east ONLY for 2025
#              (2) convert heights > 8 from cm -> meters (divide by 100)
#-------------------------------------------------


#----------------------------#
# 0) Install packages (run once, then comment out)
#----------------------------#
# install.packages(c(
#   "readxl","janitor","visreg","MASS","ggpmisc","performance","ggeffects",
#   "sjPlot","rmarkdown","knitr","ggpattern","tidyverse","lme4","lmerTest"
# ))

ÄBIN_compact <- read.csv(
  "C:/Users/shge0002/Documents/R/R/ÄBIN2025raw/2025-BINdata_raw/ÄBIN_compact",
  stringsAsFactors = FALSE
)

#----------------------------#
# 1) Load libraries (same spirit as your script, but cleaner)
#----------------------------#
library(readxl)
library(ggplot2)
library(janitor)
library(tidyverse)
library(dplyr)
library(visreg)
library(MASS)
library(ggpmisc)
library(lme4)
library(lmerTest)
library(performance)
library(ggeffects)
library(ggpattern)
library(sjPlot)
library(rmarkdown)
library(knitr)

my_colors <- c("#FDE725FF", "#56C667FF", "#238A8DFF", "#3F4788FF")

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise

#----------------------------#
# 2) LOAD YOUR NEW DATAFRAME
#----------------------------#
# IMPORTANT: this must already exist in your environment, or load it above.
# Example: ÄBIN_compact <- read_excel("path.xlsx", sheet="Raw")
Birch <- ÄBIN_compact


#----------------------------#
# 3) Clean names (so Downy.Total -> downy_total, North -> north, etc.)
#----------------------------#
Birch <- Birch %>% clean_names()


#----------------------------#
# 4) Rename / align NEW names to match your original script expectations
#    (so the rest of your script can stay the same)
#----------------------------#
Birch <- Birch %>%
  dplyr::rename(
    # originals you relied on
    area         = area,
    year         = year,
    date         = date,
    stand_number = stand,
    plot         = plot,
    surveyer     = surveyor,
    nearest_tract = nearest_tract,
    productivity = productivity,
    north        = north,
    east         = east,
    pct          = pct,
    half_height  = half_height,
    # pine names
    pine_undamaged           = pine_unbrowsed,
    pine_stems               = pine_stem,
    pine_winter_damage_stems = pine_winter_browsed_stems,
    # spruce names (best match to your later usage)
    spruce_undamaged         = spruce_ub,
    spruce_hh                = spruce_stem,
    # pellet names
    moose_pellets            = moose_pellets,
    red_deer_pellets         = red_deer_pellets,
    small_deer_pellets       = small_deer_pellets,
    reindeer_pellets         = reindeer_pellets,
    wild_boar                = wild_boar
  ) %>%
  # spruce "damaged" column in your original script:
  # your new file has winter browsed stems (closest analogue)
  dplyr::mutate(
    spruce_damaged = spruce_winter_browsed_stems
  )


#----------------------------#
# 5) Keep only the columns you used (same list as your script)
#----------------------------#
Birch_clean <- Birch[, c(
  "area", "year", "date", "stand_number", "productivity", "north", "east",
  "pct", "half_height",
  "pine_undamaged", "pine_stems", "pine_winter_damage_stems",
  "spruce_damaged", "spruce_undamaged", "spruce_hh",
  "downy_total", "downy_height", "downy_damage",
  "silver_total", "silver_height", "silver_damage",
  "rowan_total", "rowan_height", "rowan_damaged",
  "aspen_total", "aspen_height", "aspen_damaged",
  "salix_total", "salix_height", "salix_damaged",
  "oak_total", "oak_height", "oak_damaged",
  "moose_pellets", "red_deer_pellets", "small_deer_pellets",
  "reindeer_pellets", "wild_boar"
)]


#----------------------------#
# 6) Subset to 2024 + 2025
#----------------------------#
Birch_2425 <- subset(Birch_clean, year %in% c(2024, 2025))

#----------------------------#
# 6b) Fix area names (typos)
#----------------------------#
Birch_2425 <- Birch_2425 %>%
  mutate(
    area = recode(area,
                  "Lyksele"    = "Lycksele",
                  "OsterMalma" = "ÖsterMalma",
                  "Växjo"      = "Växjö"
    )
  )

cat("\n--- UNIQUE AREAS (AFTER FIX) ---\n")
print(sort(unique(Birch_2425$area)))

#----------------------------#
# 7) FIX COORDINATES
#    - 0 placeholders -> NA
#    - 2025 raw file has north=lon and east=lat (confirmed)
#    - swap ONLY rows that look reversed in 2025 (robust + safe)
#----------------------------#
Birch_2425 <- Birch_2425 %>%
  mutate(
    north = na_if(as.numeric(north), 0),
    east  = na_if(as.numeric(east),  0),
    
    # detect if columns look reversed (per row)
    north_is_lon = !is.na(north) & north >= 10 & north <= 30,  # looks like Sweden longitude
    east_is_lat  = !is.na(east)  & east  >= 54 & east  <= 70,  # looks like Sweden latitude
    
    swap_row = (year == 2025) & north_is_lon & east_is_lat
  ) %>%
  mutate(
    north_tmp = north,
    east_tmp  = east,
    north = if_else(swap_row, east_tmp,  north_tmp),
    east  = if_else(swap_row, north_tmp, east_tmp)
  ) %>%
  select(-north_tmp, -east_tmp, -north_is_lon, -east_is_lat, -swap_row)

cat("\n--- COORDINATE SUMMARY BY YEAR (AFTER CONDITIONAL SWAP) ---\n")
print(Birch_2425 %>% group_by(year) %>% summarise(
  n = n(),
  north_min = min(north, na.rm = TRUE),
  north_med = median(north, na.rm = TRUE),
  north_max = max(north, na.rm = TRUE),
  east_min  = min(east,  na.rm = TRUE),
  east_med  = median(east, na.rm = TRUE),
  east_max  = max(east,  na.rm = TRUE),
  .groups = "drop"
))

cat("\n--- SANITY CHECK (PROPORTION IN SWEDEN LAT/LON RANGES) ---\n")
print(Birch_2425 %>% group_by(year) %>% summarise(
  prop_north_lat = mean(north >= 54 & north <= 70, na.rm = TRUE),
  prop_east_lon  = mean(east  >= 10 & east  <= 30, na.rm = TRUE),
  .groups = "drop"
))

#----------------------------#
# 8) FIX HEIGHTS: if > 8 assume cm -> convert to meters (/100)
#----------------------------#
height_cols <- c(
  "half_height",
  "downy_height","silver_height",
  "rowan_height","aspen_height","salix_height","oak_height"
)

Birch_2425 <- Birch_2425 %>%
  mutate(
    across(
      all_of(height_cols),
      ~ case_when(
        is.na(.) ~ NA_real_,
        as.numeric(.) > 8 ~ as.numeric(.) / 100,
        TRUE ~ as.numeric(.)
      )
    )
  )

cat("\n--- HEIGHT SUMMARY (AFTER FIX) ---\n")
print(summary(Birch_2425[, height_cols]))

#----------------------------#
# 9) How many plots per stand? (same as your script)
#----------------------------#
table(Birch_2425$stand_number)


#----------------------------#
# 10) Filter stand numbers with more than 10 rows (same as your script)
#----------------------------#
stand_counts <- table(Birch_2425$stand_number)
stands_more_than_10 <- names(stand_counts[stand_counts > 10])
print(stands_more_than_10)


#----------------------------#
# 11) Convert browsing classes to numeric midpoints (same function as your script)
#----------------------------#
convert_damage <- function(x) {
  case_when(
    x == "0" ~ 0,
    x == "≤10" ~ 5,
    x == "11_25" ~ 18,
    x == "26_50" ~ 38,
    x == "51_75" ~ 63,
    x == "76_100" ~ 88,
    TRUE ~ NA_real_
  )
}

damage_cols <- c("downy_damage","silver_damage","rowan_damaged",
                 "aspen_damaged","salix_damaged","oak_damaged")

for (col in damage_cols) {
  new_col <- paste0(col, "_N")
  Birch_2425[[new_col]] <- convert_damage(Birch_2425[[col]])
}

head(Birch_2425)


############              PLOT LEVEL DF                    #####################
# ---------------------------- Create plot_id within each stand ----------------
# You don’t have a plot column, so we create plot_no = 1..n within each stand,
# then plot_id = "standNumber_plotNo" (e.g. "ATT2_1")
# Sorting first ensures numbering is reproducible.
Birch_2425 <- Birch_2425 %>%
  arrange(stand_number, north, east) %>%
  group_by(stand_number) %>%
  mutate(
    plot_no = row_number(),
    plot_id = paste0(stand_number, "_", plot_no)
  ) %>%
  ungroup()

# Check stands with >10 plots (shouldn’t happen, but your data shows some do)
stand_plot_counts <- Birch_2425 %>%
  count(stand_number, name = "n_plots") %>%
  filter(n_plots > 10)
print(stand_plot_counts)




###############                  STAND LEVEL DF             ###################
#----------------------------#
#     Create stand-level means dataframe: stands (same logic as your script)
#     (exclude old categorical damage columns)
#----------------------------#
exclude_cols <- c("downy_damage","silver_damage","rowan_damaged",
                  "aspen_damaged","salix_damaged","oak_damaged")

# Majority-rule helper that never returns NULL (fixes your MA12 issue)
mode_chr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  tab <- table(x)
  names(tab)[which.max(tab)]
}

stands <- Birch_2425 %>%
  dplyr::select(-all_of(exclude_cols)) %>%
  dplyr::group_by(stand_number) %>%
  dplyr::summarise(
    dplyr::across(where(is.numeric), mean, na.rm = TRUE),
    
    # Majority rule for PCT (same as your script)
    pct = ifelse(sum(pct == "N", na.rm = TRUE) > sum(pct == "Y", na.rm = TRUE), "N",
                 ifelse(sum(pct == "N", na.rm = TRUE) < sum(pct == "Y", na.rm = TRUE), "Y", "0")),
    
    # Majority rule for productivity (safe)
    productivity = mode_chr(productivity),
    
    # keep area corresponding to stand_number
    area = dplyr::first(area),
    .groups = "drop"
  )

head(stands)



