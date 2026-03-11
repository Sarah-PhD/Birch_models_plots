###############################################################################
# Project: Birch Browsing Analysis – Mixture, Latitude, Density, and Herbivores
#
# Script: birch_browsing_analysis.R
#
# Author: Sarah Gore
# Date: 11 March 2026
# Location: SLU / Umeå
#
# ---------------------------------------------------------------------------
# Description
# ---------------------------------------------------------------------------
# This script performs the full workflow for analysing browsing damage on
# downy birch (Betula pubescens) and silver birch (Betula pendula) using
# ÄBIN plot data from Sweden.
#
# The analysis investigates whether browsing severity differs between species
# depending on:
#
#   1. Birch mixture (proportion of silver vs downy birch)
#   2. Latitude (north coordinate)
#   3. Stem density
#   4. Herbivore pressure (pellet counts – optional models)
#
# Data structure:
#   - Raw data: abin_compact.csv
#   - Plot-level observations (~10 plots per stand)
#   - Stands nested within geographic areas
#   - Years analysed: 2024 and 2025
#
# Main workflow steps:
#
#   PART A – Data preparation
#       • Import raw ÄBIN dataset
#       • Clean column names
#       • Convert numeric variables
#       • Convert browsing damage classes to midpoint percentages
#       • Filter to 2024–2025
#
#   PART B – Data restructuring
#       • Standardise area names
#       • Create plot-level long dataset
#       • Create stand-year aggregated dataset
#
#   PART C – Statistical modelling
#       • Mixed-effects models for browsing severity
#       • Plot-level models: random effects for area and stand
#       • Stand-level models: random effects for area
#       • Optional pellet × mixture interaction models
#
#   PART D – Visualisation
#       • Model prediction plots
#       • Raw vs model comparisons
#       • Stand-level composition and latitude figures
#
# Key response variable:
#       Browsing severity (damage_N), derived from damage classes:
#
#           0      -> 0%
#           ≤10    -> 5%
#           11–25  -> 18%
#           26–50  -> 38%
#           51–75  -> 63%
#           76–100 -> 88%
#
# Main models test:
#
#       species × prop_s interaction
#
# which evaluates whether the proportion of silver birch alters browsing
# severity differently for downy and silver birch.
#
# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
# The script produces:
#
#   • Mixed model summaries and likelihood tests
#   • Marginal effect plots
#   • Raw vs model comparison figures
#   • Publication-ready two-panel figures
#
##############################################################################################################################################################
# PART A: IMPORT + CLEAN + PREPARE Birch_2425
###############################################################################
install.packages("broom.mixed")
install.packages("sjPlot")
library(sjPlot)
#----------------------------#
# 0) Packages
#----------------------------#
library(dplyr)
library(tidyr)
library(stringr)
library(janitor)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ggeffects)
library(ggplot2)
library(patchwork)

my_colors <- c("#FDE725FF", "#56C667FF", "#238A8DFF", "#3F4788FF")

#----------------------------#
# 1) Import raw data
#----------------------------#
ÄBIN_compact <- read.csv(
  "C:/Users/shge0002/Documents/R/R/ÄBIN2025raw/2025-BINdata_raw/abin_compact.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)



#----------------------------#
# 2) Clean names + rename key columns
#----------------------------#
Birch <- ÄBIN_compact %>%
  clean_names() %>%
  rename(
    stand_number     = stand,
    surveyer         = surveyor,
    pine_undamaged   = pine_unbrowsed,
    spruce_undamaged = spruce_ub
  ) %>%
  mutate(
    spruce_damaged = spruce_winter_damage_stems
  )

#----------------------------#
# 3) Keep only needed columns
#----------------------------#
keep_cols <- c(
  "area",
  "year",
  "date",
  "stand_number",
  "plot",
  "surveyer",
  "nearest_tract",
  "productivity",
  "north",
  "east",
  "pct",
  "half_height",
  
  "pine_undamaged",
  "pine_stems",
  "pine_winter_damage_stems",
  "proportion_pine_damage_winter",
  
  "spruce_undamaged",
  "spruce_stems",
  "spruce_damaged",
  
  "downy_total",
  "downy_height",
  "downy_damage",
  
  "silver_total",
  "silver_height",
  "silver_damage",
  
  "rowan_total",
  "rowan_height",
  "rowan_damage",
  
  "aspen_total",
  "aspen_height",
  "aspen_damage",
  
  "salix_total",
  "salix_height",
  "salix_damage",
  
  "oak_total",
  "oak_height",
  "oak_damage",
  
  "moose_pellets",
  "red_deer_pellets",
  "small_deer_pellets",
  "reindeer_pellets",
  "wild_boar"
)

Birch_clean <- Birch[, keep_cols]

#----------------------------#
# 4) Helper for numeric conversion
#----------------------------#
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

#----------------------------#
# 5) Convert ordinary numeric columns
#----------------------------#
numeric_cols <- c(
  "year",
  "north",
  "east",
  "half_height",
  "pine_undamaged",
  "pine_stems",
  "pine_winter_damage_stems",
  "proportion_pine_damage_winter",
  "spruce_undamaged",
  "spruce_stems",
  "spruce_damaged",
  "downy_total",
  "downy_height",
  "silver_total",
  "silver_height",
  "rowan_total",
  "rowan_height",
  "aspen_total",
  "aspen_height",
  "salix_total",
  "salix_height",
  "oak_total",
  "oak_height",
  "moose_pellets",
  "red_deer_pellets",
  "small_deer_pellets",
  "reindeer_pellets",
  "wild_boar"
)

Birch_clean <- Birch_clean %>%
  mutate(across(all_of(numeric_cols), num))

#----------------------------#
# 6) Convert damage classes to numeric midpoints
#----------------------------#
convert_damage <- function(x) {
  x <- as.character(x)
  
  case_when(
    x == "0"      ~ 0,
    x == "≤10"    ~ 5,
    x == "11_25"  ~ 18,
    x == "26_50"  ~ 38,
    x == "51_75"  ~ 63,
    x == "76_100" ~ 88,
    TRUE ~ NA_real_
  )
}

damage_cols <- c(
  "downy_damage",
  "silver_damage",
  "rowan_damage",
  "aspen_damage",
  "salix_damage",
  "oak_damage"
)

Birch_clean <- Birch_clean %>%
  mutate(
    across(
      all_of(damage_cols),
      convert_damage,
      .names = "{.col}_N"
    )
  )

# checks
table(Birch_clean$downy_damage, Birch_clean$downy_damage_N, useNA = "ifany")
table(Birch_clean$silver_damage, Birch_clean$silver_damage_N, useNA = "ifany")

#----------------------------#
# 7) Restrict to 2024 and 2025
#----------------------------#
Birch_2425 <- Birch_clean %>%
  filter(year %in% c(2024, 2025))

table(Birch_2425$year, useNA = "ifany")
table(Birch_2425$downy_damage_N, useNA = "ifany")
table(Birch_2425$silver_damage_N, useNA = "ifany")


###############################################################################
# CHECK FOR SWAPPED EAST/NORTH COORDINATES AND FIX THEM
###############################################################################

# Expected rough SWEREF 99 TM ranges for Sweden
east_min  <- 200000
east_max  <- 1000000
north_min <- 6000000
north_max <- 8000000

#----------------------------#
# 1) Flag rows that look swapped
#----------------------------#
Birch_2425 <- Birch_2425 %>%
  mutate(
    coord_swap_flag = ifelse(
      !is.na(north) & !is.na(east) &
        north >= east_min  & north <= east_max &
        east  >= north_min & east  <= north_max,
      TRUE, FALSE
    )
  )

# How many suspected swaps?
table(Birch_2425$coord_swap_flag, useNA = "ifany")

# Inspect them
Birch_2425 %>%
  filter(coord_swap_flag) %>%
  select(area, stand_number, plot, east, north) %>%
  arrange(area, stand_number, plot)

#----------------------------#
# 2) Fix swapped rows
#----------------------------#
Birch_2425 <- Birch_2425 %>%
  mutate(
    east_fixed = ifelse(coord_swap_flag, north, east),
    north_fixed = ifelse(coord_swap_flag, east, north)
  ) %>%
  select(-east, -north) %>%
  rename(
    east = east_fixed,
    north = north_fixed
  )

#----------------------------#
# 3) Recheck ranges after fixing
#----------------------------#
summary(Birch_2425$east)
summary(Birch_2425$north)

quantile(Birch_2425$east, probs = c(0, .01, .05, .25, .5, .75, .95, .99, 1), na.rm = TRUE)
quantile(Birch_2425$north, probs = c(0, .01, .05, .25, .5, .75, .95, .99, 1), na.rm = TRUE)



#----------------------------#
# 4) Look for any remaining suspicious rows
#----------------------------#
Birch_2425 %>%
  filter(
    !(east >= east_min & east <= east_max & north >= north_min & north <= north_max)
  ) %>%
  select(area, stand_number, plot, east, north) %>%
  arrange(area, stand_number, plot)


###############################################################################
# FIX WGS84 EAST/NORTH COORDINATES AND CREATE CLEAN NORTH VARIABLE
###############################################################################

# Turn 0,0 into missing
Birch_2425 <- Birch_2425 %>%
  mutate(
    east  = ifelse(east == 0, NA, east),
    north = ifelse(north == 0, NA, north)
  )

# Classify coordinate pattern
Birch_2425 <- Birch_2425 %>%
  mutate(
    coord_type = case_when(
      is.na(east) | is.na(north) ~ "missing",
      
      # correct WGS84: lon, lat
      east >= 10 & east <= 25 & north >= 55 & north <= 70 ~ "wgs84_ok",
      
      # flipped WGS84: lat, lon
      east >= 55 & east <= 70 & north >= 10 & north <= 25 ~ "wgs84_swapped",
      
      TRUE ~ "other"
    )
  )

table(Birch_2425$coord_type, useNA = "ifany")

# Keep originals and fix flipped rows
Birch_2425 <- Birch_2425 %>%
  mutate(
    east_orig  = east,
    north_orig = north,
    east = ifelse(coord_type == "wgs84_swapped", north_orig, east_orig),
    north = ifelse(coord_type == "wgs84_swapped", east_orig, north_orig)
  )

# Recheck classification after fix
Birch_2425 <- Birch_2425 %>%
  mutate(
    coord_type_fixed = case_when(
      is.na(east) | is.na(north) ~ "missing",
      east >= 10 & east <= 25 & north >= 55 & north <= 70 ~ "wgs84",
      TRUE ~ "other"
    )
  )

table(Birch_2425$coord_type_fixed, useNA = "ifany")

# Create modelling north variable from cleaned latitude
Birch_2425 <- Birch_2425 %>%
  mutate(
    north_model = ifelse(coord_type_fixed == "wgs84", north, NA_real_)
  )

summary(Birch_2425$north_model)

# Quick inspection of unresolved rows
Birch_2425 %>%
  filter(coord_type_fixed != "wgs84") %>%
  select(area, stand_number, plot, east_orig, north_orig, east, north, coord_type, coord_type_fixed) %>%
  arrange(area, stand_number, plot)
##############################################################################
#                 CONTINUE ...
##############################################################################


###############################################################################
# PART B: STREAMLINED WORKFLOW: Birch mixture × latitude × density × pellets
###############################################################################

#----------------------------#
# 1) AREA name cleaning
#----------------------------#
standardise_area <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_squish() %>%
    str_replace_all("[\u2010-\u2015]", "-") %>%
    str_replace_all("[_]", "") %>%
    str_replace_all("\\s+", "")
}

area_map <- c(
  "Lyksele"     = "Lycksele",
  "LYKSELE"     = "Lycksele",
  "lycksele"    = "Lycksele",
  "ÖsterMalma"  = "ÖsterMalma",
  "OsterMalma"  = "ÖsterMalma",
  "Ostermalma"  = "ÖsterMalma",
  "Östermalma"  = "ÖsterMalma",
  "Öster Malma" = "ÖsterMalma",
  "Oster Malma" = "ÖsterMalma"
)

fix_area_names <- function(df, area_col = area) {
  df %>%
    mutate(
      area_raw = {{ area_col }},
      area_std = standardise_area({{ area_col }}),
      area = dplyr::recode(area_raw, !!!area_map, .default = area_raw),
      area = dplyr::recode(area, !!!area_map, .default = area),
      area = standardise_area(area)
    ) %>%
    select(-area_raw, -area_std)
}

Birch_2425 <- Birch_2425 %>%
  fix_area_names(area)

sort(unique(Birch_2425$area))
###############################################################################
# PREPARE PELLET VARIABLES FOR MODELLING
# Add this after Birch_2425 is created, before build_plot_long()
###############################################################################

#----------------------------#
# 1) Make sure raw pellet columns are numeric
#----------------------------#
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

Birch_2425 <- Birch_2425 %>%
  mutate(
    moose_pellets      = num(moose_pellets),
    red_deer_pellets   = num(red_deer_pellets),
    small_deer_pellets = num(small_deer_pellets),
    reindeer_pellets   = num(reindeer_pellets)
  )

#----------------------------#
# 2) Replace missing pellet values with 0
#----------------------------#
Birch_2425 <- Birch_2425 %>%
  mutate(
    moose_pellets      = ifelse(is.na(moose_pellets), 0, moose_pellets),
    red_deer_pellets   = ifelse(is.na(red_deer_pellets), 0, red_deer_pellets),
    small_deer_pellets = ifelse(is.na(small_deer_pellets), 0, small_deer_pellets),
    reindeer_pellets   = ifelse(is.na(reindeer_pellets), 0, reindeer_pellets)
  )

#----------------------------#
# 3) Create pellet indices
#----------------------------#
Birch_2425 <- Birch_2425 %>%
  mutate(
    # Moose only
    moose_pellets_only = moose_pellets,
    moose_pellets_log  = log1p(moose_pellets_only),
    moose_present      = ifelse(moose_pellets_only > 0, 1, 0),
    
    # Red deer + small deer combined
    deer_pellets       = red_deer_pellets + small_deer_pellets,
    deer_pellets_log   = log1p(deer_pellets),
    deer_present       = ifelse(deer_pellets > 0, 1, 0),
    
    # Optional all cervids together
    total_cervid_pellets     = moose_pellets + red_deer_pellets + small_deer_pellets,
    total_cervid_pellets_log = log1p(total_cervid_pellets),
    cervid_present           = ifelse(total_cervid_pellets > 0, 1, 0)
  )

#----------------------------#
# 4) Quick checks
#----------------------------#
summary(Birch_2425$moose_pellets_only)
summary(Birch_2425$deer_pellets)
summary(Birch_2425$total_cervid_pellets)

table(Birch_2425$moose_present, useNA = "ifany")
table(Birch_2425$deer_present, useNA = "ifany")
table(Birch_2425$cervid_present, useNA = "ifany")

quantile(Birch_2425$moose_pellets_only, probs = c(0, .25, .5, .75, .9, .95, .99, 1), na.rm = TRUE)
quantile(Birch_2425$deer_pellets, probs = c(0, .25, .5, .75, .9, .95, .99, 1), na.rm = TRUE)



#----------------------------#
# 2) Build PLOT-level long dataset
#----------------------------#
build_plot_long <- function(df, include_pellets = FALSE) {
  
  if (include_pellets) {
    df <- df %>%
      mutate(
        moose_pellets      = num(moose_pellets),
        red_deer_pellets   = num(red_deer_pellets),
        small_deer_pellets = num(small_deer_pellets),
        reindeer_pellets   = if ("reindeer_pellets" %in% names(df)) num(reindeer_pellets) else 0
      ) %>%
      mutate(
        moose_pellets      = ifelse(is.na(moose_pellets), 0, moose_pellets),
        red_deer_pellets   = ifelse(is.na(red_deer_pellets), 0, red_deer_pellets),
        small_deer_pellets = ifelse(is.na(small_deer_pellets), 0, small_deer_pellets),
        reindeer_pellets   = ifelse(is.na(reindeer_pellets), 0, reindeer_pellets),
        total_pellets      = moose_pellets + red_deer_pellets + small_deer_pellets + reindeer_pellets,
        total_pellets_log  = log1p(total_pellets),
        pellet_present     = ifelse(total_pellets > 0, 1, 0)
      )
  }
  
  build_plot_long <- function(df, include_pellets = FALSE) {
    
    wide <- df %>%
      transmute(
        year         = as.factor(year),
        area         = area,
        stand_number = as.factor(stand_number),
        north        = num(north_model),
        north_c      = as.numeric(scale(north, center = TRUE, scale = FALSE)),
        
        downy_total     = num(downy_total),
        silver_total    = num(silver_total),
        downy_damage_N  = num(downy_damage_N),
        silver_damage_N = num(silver_damage_N),
        
        moose_pellets_log        = if (include_pellets) num(moose_pellets_log) else NULL,
        deer_pellets_log         = if (include_pellets) num(deer_pellets_log) else NULL,
        total_cervid_pellets_log = if (include_pellets) num(total_cervid_pellets_log) else NULL,
        moose_present            = if (include_pellets) num(moose_present) else NULL,
        deer_present             = if (include_pellets) num(deer_present) else NULL,
        cervid_present           = if (include_pellets) num(cervid_present) else NULL
      ) %>%
      mutate(
        total_birch = downy_total + silver_total,
        prop_s_plot = ifelse(total_birch > 0, silver_total / total_birch, NA_real_)
      )
    
    long <- wide %>%
      pivot_longer(
        cols = c(downy_total, silver_total, downy_damage_N, silver_damage_N),
        names_to = c("sp", ".value"),
        names_pattern = "(downy|silver)_(total|damage_N)"
      ) %>%
      mutate(
        species    = ifelse(sp == "downy", "Downy", "Silver"),
        stem_count = total
      ) %>%
      select(
        year, area, stand_number,
        north_c,
        species, damage_N, stem_count,
        prop_s_plot,
        any_of(c(
          "moose_pellets_log",
          "deer_pellets_log",
          "total_cervid_pellets_log",
          "moose_present",
          "deer_present",
          "cervid_present"
        ))
      ) %>%
      filter(
        !is.na(year),
        !is.na(area),
        !is.na(stand_number),
        !is.na(damage_N),
        !is.na(stem_count),
        !is.na(prop_s_plot)
      )
    
    long
  }

build_plot_long <- function(df, include_pellets = FALSE) {

  wide <- df %>%
    transmute(
      year         = as.factor(year),
      area         = area,
      stand_number = as.factor(stand_number),
      north        = num(north_model),
      north_c      = as.numeric(scale(north, center = TRUE, scale = FALSE)),

      downy_total     = num(downy_total),
      silver_total    = num(silver_total),
      downy_damage_N  = num(downy_damage_N),
      silver_damage_N = num(silver_damage_N),

      moose_pellets_log        = if (include_pellets) num(moose_pellets_log) else NULL,
      deer_pellets_log         = if (include_pellets) num(deer_pellets_log) else NULL,
      total_cervid_pellets_log = if (include_pellets) num(total_cervid_pellets_log) else NULL,
      moose_present            = if (include_pellets) num(moose_present) else NULL,
      deer_present             = if (include_pellets) num(deer_present) else NULL,
      cervid_present           = if (include_pellets) num(cervid_present) else NULL
    ) %>%
    mutate(
      total_birch = downy_total + silver_total,
      prop_s_plot = ifelse(total_birch > 0, silver_total / total_birch, NA_real_)
    )

  long <- wide %>%
    pivot_longer(
      cols = c(downy_total, silver_total, downy_damage_N, silver_damage_N),
      names_to = c("sp", ".value"),
      names_pattern = "(downy|silver)_(total|damage_N)"
    ) %>%
    mutate(
      species    = ifelse(sp == "downy", "Downy", "Silver"),
      stem_count = total
    ) %>%
    select(
      year, area, stand_number,
      north_c,
      species, damage_N, stem_count,
      prop_s_plot,
      any_of(c(
        "moose_pellets_log",
        "deer_pellets_log",
        "total_cervid_pellets_log",
        "moose_present",
        "deer_present",
        "cervid_present"
      ))
    ) %>%
    filter(
      !is.na(year),
      !is.na(area),
      !is.na(stand_number),
      !is.na(damage_N),
      !is.na(stem_count),
      !is.na(prop_s_plot)
    )

  long
}
birch_plot_long <- build_plot_long(Birch_2425, include_pellets = FALSE)
birch_plot_long_pel <- build_plot_long(Birch_2425, include_pellets = TRUE)
#----------------------------#
# 3) Add unique stand ID + centered stem count
#    THIS is the extra block to include
#----------------------------#
birch_plot_long <- birch_plot_long %>%
  mutate(
    stand_id = interaction(area, stand_number, drop = TRUE)
  ) %>%
  group_by(species) %>%
  mutate(
    stem_count_c = as.numeric(scale(stem_count, center = TRUE, scale = FALSE))
  ) %>%
  ungroup()

birch_plot_long_pel <- birch_plot_long_pel %>%
  mutate(
    stand_id = interaction(area, stand_number, drop = TRUE)
  ) %>%
  group_by(species) %>%
  mutate(
    stem_count_c = as.numeric(scale(stem_count, center = TRUE, scale = FALSE))
  ) %>%
  ungroup()

summary(birch_plot_long$damage_N)
table(birch_plot_long$damage_N, useNA = "ifany")

#----------------------------#
# 4) Build STAND-year long dataset
#----------------------------#
build_stand_long <- function(plot_long, include_pellets = FALSE) {
  
  stand_wide <- plot_long %>%
    group_by(year, area, stand_number) %>%
    summarise(
      mean_north_c = mean(north_c, na.rm = TRUE),
      
      downy_damage_mean  = mean(damage_N[species == "Downy"],  na.rm = TRUE),
      silver_damage_mean = mean(damage_N[species == "Silver"], na.rm = TRUE),
      
      downy_stems_mean   = mean(stem_count[species == "Downy"],  na.rm = TRUE),
      silver_stems_mean  = mean(stem_count[species == "Silver"], na.rm = TRUE),
      
      moose_pellets_log_mean        = if (include_pellets) mean(moose_pellets_log, na.rm = TRUE) else NA_real_,
      deer_pellets_log_mean         = if (include_pellets) mean(deer_pellets_log, na.rm = TRUE) else NA_real_,
      total_cervid_pellets_log_mean = if (include_pellets) mean(total_cervid_pellets_log, na.rm = TRUE) else NA_real_,
      
      moose_present_any  = if (include_pellets) max(moose_present, na.rm = TRUE) else NA_real_,
      deer_present_any   = if (include_pellets) max(deer_present, na.rm = TRUE) else NA_real_,
      cervid_present_any = if (include_pellets) max(cervid_present, na.rm = TRUE) else NA_real_,
      
      .groups = "drop"
    ) %>%
    mutate(
      total_stems_mean = downy_stems_mean + silver_stems_mean,
      prop_s_stand = ifelse(total_stems_mean > 0, silver_stems_mean / total_stems_mean, NA_real_)
    )
  
  stand_long <- stand_wide %>%
    pivot_longer(
      cols = c(downy_stems_mean, silver_stems_mean,
               downy_damage_mean, silver_damage_mean),
      names_to = c("sp", ".value"),
      names_pattern = "(downy|silver)_(stems_mean|damage_mean)"
    ) %>%
    mutate(
      species    = ifelse(sp == "downy", "Downy", "Silver"),
      stem_count = stems_mean,
      damage_N   = damage_mean
    ) %>%
    select(
      year, area, stand_number,
      mean_north_c,
      species, damage_N, stem_count,
      prop_s_stand,
      any_of(c(
        "moose_pellets_log_mean",
        "deer_pellets_log_mean",
        "total_cervid_pellets_log_mean",
        "moose_present_any",
        "deer_present_any",
        "cervid_present_any"
      ))
    ) %>%
    filter(
      !is.na(damage_N),
      !is.na(stem_count),
      !is.na(prop_s_stand)
    )
  
  stand_long
}

birch_stand_long <- build_stand_long(birch_plot_long, include_pellets = FALSE)
birch_stand_long_pel <- build_stand_long(birch_plot_long_pel, include_pellets = TRUE)

birch_stand_long <- birch_stand_long %>%
  group_by(species) %>%
  mutate(
    stem_count_c = as.numeric(scale(stem_count, center = TRUE, scale = FALSE))
  ) %>%
  ungroup()

birch_stand_long_pel <- birch_stand_long_pel %>%
  group_by(species) %>%
  mutate(
    stem_count_c = as.numeric(scale(stem_count, center = TRUE, scale = FALSE))
  ) %>%
  ungroup()

#----------------------------#
# 5) CORE MODELS
#----------------------------#

# Plot-level model
m_plot <- lmer(
  damage_N ~ species * prop_s_plot +
    species * north_c +
    species * stem_count_c +
    year +
    (1 | area) +
    (1 | stand_id),
  data = birch_plot_long,
  REML = TRUE
)

# Stand-level model
m_stand <- lmer(
  damage_N ~ species * prop_s_stand +
    species * mean_north_c +
    species * stem_count_c +
    year +
    (1 | area),
  data = birch_stand_long,
  REML = TRUE
)

summary(m_plot)
r.squaredGLMM(m_plot)
drop1(update(m_plot, REML = FALSE), test = "Chisq")

summary(m_stand)
r.squaredGLMM(m_stand)
drop1(update(m_stand, REML = FALSE), test = "Chisq")

#----------------------------#
# 6) OPTIONAL pellets models
#----------------------------#
birch_downy_plot_pel <- birch_plot_long_pel %>% filter(species == "Downy")
birch_downy_stand_pel <- birch_stand_long_pel %>% filter(species == "Downy")

m_downy_plot_pel <- lmer(
  damage_N ~ prop_s_plot * total_pellets_log +
    north_c + stem_count_c + year +
    (1 | area) + (1 | stand_id),
  data = birch_downy_plot_pel,
  REML = TRUE
)

m_downy_stand_pel <- lmer(
  damage_N ~ prop_s_stand * pellets_log_mean +
    mean_north_c + stem_count_c + year +
    (1 | area),
  data = birch_downy_stand_pel,
  REML = TRUE
)

summary(m_downy_plot_pel)
drop1(update(m_downy_plot_pel, REML = FALSE), test = "Chisq")

summary(m_downy_stand_pel)
drop1(update(m_downy_stand_pel, REML = FALSE), test = "Chisq")

# species model with pellets
Sp <- lmer(
  damage_N ~ species * prop_s_plot * total_pellets_log +
    north_c + stem_count_c + year +
    (1 | area) + (1 | stand_id),
  data = birch_plot_long_pel,
  REML = TRUE
)



###############################################################################
# PART C: FIGURES + PUBLICATION TABLES
# Add this at the END of the script, after all models are fitted
###############################################################################

#----------------------------#
# 0) Packages
#----------------------------#
library(dplyr)
library(ggplot2)
library(ggeffects)
library(patchwork)
library(broom.mixed)
library(lme4)

#----------------------------#
# 1) Figure theme and palette
#----------------------------#
pal <- c("Downy" = "#3BA64C", "Silver" = "#4C5C8C")

theme_pub <- function() {
  theme_classic(base_size = 13) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 12),
      plot.title.position = "plot",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

#----------------------------#
# 2) RAW FIGURES
#----------------------------#

# 2A) Raw damage by species (plot-level)
p_raw_box <- ggplot(
  birch_plot_long,
  aes(x = species, y = damage_N, fill = species)
) +
  geom_boxplot(alpha = 0.75, outlier.alpha = 0.25) +
  scale_fill_manual(values = pal) +
  labs(
    title = "A  Browsing damage by species",
    x = "Species",
    y = "Browsing damage"
  ) +
  theme_pub()

p_raw_box

# 2B) Raw damage vs proportion silver (plot-level)
p_raw_prop_plot <- ggplot(
  birch_plot_long,
  aes(x = prop_s_plot, y = damage_N, colour = species)
) +
  geom_point(alpha = 0.30, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.0) +
  scale_colour_manual(values = pal) +
  labs(
    title = "B  Raw relationship: damage vs proportion silver",
    x = "Proportion silver (plot)",
    y = "Browsing damage"
  ) +
  theme_pub()

p_raw_prop_plot

# 2C) Raw damage vs north (plot-level)
p_raw_north_plot <- ggplot(
  birch_plot_long,
  aes(x = north_c, y = damage_N, colour = species)
) +
  geom_point(alpha = 0.30, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.0) +
  scale_colour_manual(values = pal) +
  labs(
    title = "C  Raw relationship: damage vs latitude",
    x = "North coordinate (centered)",
    y = "Browsing damage"
  ) +
  theme_pub()

p_raw_north_plot

#----------------------------#
# 3) MODEL-BASED FIGURES
#----------------------------#

# 3A) Plot-level composition effect
pred_plot_prop <- ggpredict(
  m_plot,
  terms = c("prop_s_plot [0:1 by=0.05]", "species")
) |> as.data.frame()

p_model_prop_plot <- ggplot() +
  geom_point(
    data = birch_plot_long,
    aes(x = prop_s_plot, y = damage_N, colour = species),
    alpha = 0.20, size = 1.4
  ) +
  geom_ribbon(
    data = pred_plot_prop,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.15, colour = NA
  ) +
  geom_line(
    data = pred_plot_prop,
    aes(x = x, y = predicted, colour = group),
    linewidth = 1.2
  ) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title = "A  Modelled effect of birch mixture (plot level)",
    x = "Proportion silver (plot)",
    y = "Predicted browsing damage"
  ) +
  theme_pub()

p_model_prop_plot

# 3B) Plot-level north effect
pred_plot_north <- ggpredict(
  m_plot,
  terms = c("north_c [all]", "species")
) |> as.data.frame()

p_model_north_plot <- ggplot() +
  geom_point(
    data = birch_plot_long,
    aes(x = north_c, y = damage_N, colour = species),
    alpha = 0.20, size = 1.4
  ) +
  geom_ribbon(
    data = pred_plot_north,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.15, colour = NA
  ) +
  geom_line(
    data = pred_plot_north,
    aes(x = x, y = predicted, colour = group),
    linewidth = 1.2
  ) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title = "B  Modelled latitude effect (plot level)",
    x = "North coordinate (centered)",
    y = "Predicted browsing damage"
  ) +
  theme_pub()

p_model_north_plot

# 3C) Stand-level composition effect
pred_stand_prop <- ggpredict(
  m_stand,
  terms = c("prop_s_stand [0:1 by=0.05]", "species")
) |> as.data.frame()

p_model_prop_stand <- ggplot() +
  geom_point(
    data = birch_stand_long,
    aes(x = prop_s_stand, y = damage_N, colour = species),
    alpha = 0.25, size = 1.6
  ) +
  geom_ribbon(
    data = pred_stand_prop,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.15, colour = NA
  ) +
  geom_line(
    data = pred_stand_prop,
    aes(x = x, y = predicted, colour = group),
    linewidth = 1.2
  ) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title = "C  Modelled effect of birch mixture (stand level)",
    x = "Proportion silver (stand)",
    y = "Predicted browsing damage"
  ) +
  theme_pub()

p_model_prop_stand

# 3D) Stand-level north effect
pred_stand_north <- ggpredict(
  m_stand,
  terms = c("mean_north_c [all]", "species")
) |> as.data.frame()

p_model_north_stand <- ggplot() +
  geom_point(
    data = birch_stand_long,
    aes(x = mean_north_c, y = damage_N, colour = species),
    alpha = 0.25, size = 1.6
  ) +
  geom_ribbon(
    data = pred_stand_north,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.15, colour = NA
  ) +
  geom_line(
    data = pred_stand_north,
    aes(x = x, y = predicted, colour = group),
    linewidth = 1.2
  ) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title = "D  Modelled latitude effect (stand level)",
    x = "North coordinate (centered; stand mean)",
    y = "Predicted browsing damage"
  ) +
  theme_pub()

p_model_north_stand

# 3E) Combine the two main stand-level figures
fig_stand_main <- p_model_prop_stand + p_model_north_stand +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

fig_stand_main

# 3F) Downy-only pellet interaction (plot-level)
pred_downy_plot_pel <- ggpredict(
  m_downy_plot_pel,
  terms = c("prop_s_plot [0:1 by=0.05]", "total_pellets_log [0,0.7,1.4]")
) |> as.data.frame()

p_downy_plot_pel <- ggplot(pred_downy_plot_pel,
                           aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  labs(
    title = "E  Downy birch: mixture × herbivore pressure (plot level)",
    x = "Proportion silver (plot)",
    y = "Predicted browsing damage",
    colour = "Pellets\n(log1p)",
    fill   = "Pellets\n(log1p)"
  ) +
  theme_pub()

p_downy_plot_pel

# 3G) Downy-only pellet interaction (stand-level)
pred_downy_stand_pel <- ggpredict(
  m_downy_stand_pel,
  terms = c("prop_s_stand [0:1 by=0.05]", "pellets_log_mean [0,0.7,1.4]")
) |> as.data.frame()

p_downy_stand_pel <- ggplot(pred_downy_stand_pel,
                            aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  labs(
    title = "F  Downy birch: mixture × herbivore pressure (stand level)",
    x = "Proportion silver (stand)",
    y = "Predicted browsing damage",
    colour = "Pellets\n(log1p)",
    fill   = "Pellets\n(log1p)"
  ) +
  theme_pub()

p_downy_stand_pel

# Optional save commands
# ggsave("Fig_stand_main.png", fig_stand_main, width = 11, height = 4.8, dpi = 600)
# ggsave("Fig_stand_main.pdf", fig_stand_main, width = 11, height = 4.8)
# ggsave("Fig_downy_plot_pel.png", p_downy_plot_pel, width = 6.5, height = 4.8, dpi = 600)
# ggsave("Fig_downy_stand_pel.png", p_downy_stand_pel, width = 6.5, height = 4.8, dpi = 600)

#----------------------------#
# 4) PUBLICATION-READY MODEL TABLES
#----------------------------#

# Helper function for nice p-value formatting
fmt_p <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

# Helper to convert broom output to publication-style fixed effect table
make_fixed_table <- function(model, model_name) {
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      Model = model_name,
      Estimate = round(estimate, 3),
      `Std. Error` = round(std.error, 3),
      `95% CI low` = round(conf.low, 3),
      `95% CI high` = round(conf.high, 3),
      Statistic = round(statistic, 3),
      `P value` = fmt_p(p.value)
    ) %>%
    select(
      Model, term, Estimate, `Std. Error`, `95% CI low`, `95% CI high`,
      Statistic, `P value`
    )
}

# Helper for model summary table
make_model_summary <- function(model, model_name, data_obj) {
  data.frame(
    Model = model_name,
    N = nobs(model),
    Groups = paste(names(lme4::getME(model, "flist")), collapse = ", "),
    AIC = round(AIC(model), 1),
    BIC = round(BIC(model), 1),
    LogLik = round(as.numeric(logLik(model)), 2),
    stringsAsFactors = FALSE
  )
}

#----------------------------#
# 5) Fixed-effects tables
#----------------------------#
tab_m_plot            <- make_fixed_table(m_plot, "Plot-level combined species model")
tab_m_stand           <- make_fixed_table(m_stand, "Stand-level combined species model")
tab_m_downy_plot_pel  <- make_fixed_table(m_downy_plot_pel, "Downy plot-level pellets model")
tab_m_downy_stand_pel <- make_fixed_table(m_downy_stand_pel, "Downy stand-level pellets model")
tab_Sp                <- make_fixed_table(Sp, "Species × mixture × pellets model")

# Print to console for copy/paste
tab_m_plot
tab_m_stand
tab_m_downy_plot_pel
tab_m_downy_stand_pel
tab_Sp

#----------------------------#
# 6) Model overview table
#----------------------------#
tab_model_overview <- bind_rows(
  make_model_summary(m_plot, "Plot-level combined species model", birch_plot_long),
  make_model_summary(m_stand, "Stand-level combined species model", birch_stand_long),
  make_model_summary(m_downy_plot_pel, "Downy plot-level pellets model", birch_plot_long_pel),
  make_model_summary(m_downy_stand_pel, "Downy stand-level pellets model", birch_stand_long_pel),
  make_model_summary(Sp, "Species × mixture × pellets model", birch_plot_long_pel)
)

tab_model_overview

#----------------------------#
# 7) Compact publication tables for Word
#----------------------------#

# Main paper table: only the two core models
tab_main_results <- bind_rows(
  tab_m_plot,
  tab_m_stand
)

tab_main_results

# Supplementary table: pellet models
tab_pellet_results <- bind_rows(
  tab_m_downy_plot_pel,
  tab_m_downy_stand_pel,
  tab_Sp
)

tab_pellet_results

#----------------------------#
# 8) Optional CSV export for Word / Excel
#----------------------------#
# These are very easy to open and paste into a manuscript table
write.csv(tab_main_results, "tab_main_results.csv", row.names = FALSE)
write.csv(tab_pellet_results, "tab_pellet_results.csv", row.names = FALSE)
write.csv(tab_model_overview, "tab_model_overview.csv", row.names = FALSE)

#----------------------------#
# 9) Optional cleaner term labels for publication
#----------------------------#
clean_terms <- function(df) {
  df %>%
    mutate(
      term = recode(
        term,
        "(Intercept)" = "Intercept",
        "speciesSilver" = "Species: Silver",
        "prop_s_plot" = "Proportion silver (plot)",
        "prop_s_stand" = "Proportion silver (stand)",
        "north_c" = "North (plot-centered)",
        "mean_north_c" = "North (stand mean, centered)",
        "stem_count_c" = "Stem count (centered)",
        "total_pellets_log" = "Total pellets (log1p)",
        "pellets_log_mean" = "Mean pellets (log1p)",
        "year2025" = "Year: 2025",
        "speciesSilver:prop_s_plot" = "Silver × proportion silver (plot)",
        "speciesSilver:prop_s_stand" = "Silver × proportion silver (stand)",
        "speciesSilver:north_c" = "Silver × north (plot)",
        "speciesSilver:mean_north_c" = "Silver × north (stand)",
        "speciesSilver:stem_count_c" = "Silver × stem count",
        "prop_s_plot:total_pellets_log" = "Proportion silver × pellets (plot)",
        "prop_s_stand:pellets_log_mean" = "Proportion silver × pellets (stand)",
        "speciesSilver:prop_s_plot:total_pellets_log" = "Silver × proportion silver × pellets"
      )
    )
}

tab_main_results_clean  <- clean_terms(tab_main_results)
tab_pellet_results_clean <- clean_terms(tab_pellet_results)

tab_main_results_clean
tab_pellet_results_clean

write.csv(tab_main_results_clean, "tab_main_results_clean.csv", row.names = FALSE)
write.csv(tab_pellet_results_clean, "tab_pellet_results_clean.csv", row.names = FALSE)


tab_model(
  m_plot,
  m_stand,
  show.ci = TRUE,
  show.se = TRUE,
  show.stat = TRUE,
  show.re.var = TRUE,
  dv.labels = c("Plot-level model", "Stand-level model"),
  file = "model_table_main.html"
)

###############################################################################
# PELLET MODELS: MOOSE ONLY + DEER COMBINED
###############################################################################

#----------------------------#
# 1) Make species-specific subsets
#----------------------------#
birch_downy_plot_pel <- birch_plot_long_pel %>%
  filter(species == "Downy")

birch_downy_stand_pel <- birch_stand_long_pel %>%
  filter(species == "Downy")

#----------------------------#
# 2) MOOSE-ONLY MODELS
#----------------------------#

# Plot level
m_downy_plot_moose <- lmer(
  damage_N ~ prop_s_plot * moose_pellets_log +
    north_c + stem_count_c + year +
    (1 | area) + (1 | stand_id),
  data = birch_downy_plot_pel,
  REML = TRUE
)

# Stand level
m_downy_stand_moose <- lmer(
  damage_N ~ prop_s_stand * moose_pellets_log_mean +
    mean_north_c + stem_count_c + year +
    (1 | area),
  data = birch_downy_stand_pel,
  REML = TRUE
)

summary(m_downy_plot_moose)
r.squaredGLMM(m_downy_plot_moose)
drop1(update(m_downy_plot_moose, REML = FALSE), test = "Chisq")

summary(m_downy_stand_moose)
r.squaredGLMM(m_downy_stand_moose)
drop1(update(m_downy_stand_moose, REML = FALSE), test = "Chisq")

#----------------------------#
# 3) RED DEER + SMALL DEER COMBINED MODELS
#----------------------------#

# Plot level
m_downy_plot_deer <- lmer(
  damage_N ~ prop_s_plot * deer_pellets_log +
    north_c + stem_count_c + year +
    (1 | area) + (1 | stand_id),
  data = birch_downy_plot_pel,
  REML = TRUE
)

# Stand level
m_downy_stand_deer <- lmer(
  damage_N ~ prop_s_stand * deer_pellets_log_mean +
    mean_north_c + stem_count_c + year +
    (1 | area),
  data = birch_downy_stand_pel,
  REML = TRUE
)

summary(m_downy_plot_deer)
r.squaredGLMM(m_downy_plot_deer)
drop1(update(m_downy_plot_deer, REML = FALSE), test = "Chisq")

summary(m_downy_stand_deer)
r.squaredGLMM(m_downy_stand_deer)
drop1(update(m_downy_stand_deer, REML = FALSE), test = "Chisq")

#----------------------------#
# 4) Prediction plots
#----------------------------#
pal_pel <- c("Low" = "#4C5C8C", "Medium" = "#56C667", "High" = "#3BA64C")

theme_pub <- function() {
  theme_classic(base_size = 13) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 12),
      plot.title.position = "plot",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

#----------------------------#
# 4A) Moose-only plot-level figure
#----------------------------#
pred_moose_plot <- ggpredict(
  m_downy_plot_moose,
  terms = c("prop_s_plot [0:1 by=0.05]", "moose_pellets_log [0,0.7,1.4]")
) %>% as.data.frame()

pred_moose_plot$group <- factor(
  pred_moose_plot$group,
  labels = c("Low", "Medium", "High")
)

p_moose_plot <- ggplot(
  pred_moose_plot,
  aes(x = x, y = predicted, colour = group, fill = group)
) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  scale_colour_manual(values = pal_pel) +
  scale_fill_manual(values = pal_pel) +
  labs(
    title = "Downy birch: mixture × moose pellet density (plot level)",
    x = "Proportion silver (plot)",
    y = "Predicted browsing damage",
    colour = "Moose pellets",
    fill   = "Moose pellets"
  ) +
  theme_pub()

p_moose_plot

#----------------------------#
# 4B) Moose-only stand-level figure
#----------------------------#
pred_moose_stand <- ggpredict(
  m_downy_stand_moose,
  terms = c("prop_s_stand [0:1 by=0.05]", "moose_pellets_log_mean [0,0.7,1.4]")
) %>% as.data.frame()

pred_moose_stand$group <- factor(
  pred_moose_stand$group,
  labels = c("Low", "Medium", "High")
)

p_moose_stand <- ggplot(
  pred_moose_stand,
  aes(x = x, y = predicted, colour = group, fill = group)
) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  scale_colour_manual(values = pal_pel) +
  scale_fill_manual(values = pal_pel) +
  labs(
    title = "Downy birch: mixture × moose pellet density (stand level)",
    x = "Proportion silver (stand)",
    y = "Predicted browsing damage",
    colour = "Moose pellets",
    fill   = "Moose pellets"
  ) +
  theme_pub()

p_moose_stand

#----------------------------#
# 4C) Deer-combined plot-level figure
#----------------------------#
pred_deer_plot <- ggpredict(
  m_downy_plot_deer,
  terms = c("prop_s_plot [0:1 by=0.05]", "deer_pellets_log [0,0.7,1.4]")
) %>% as.data.frame()

pred_deer_plot$group <- factor(
  pred_deer_plot$group,
  labels = c("Low", "Medium", "High")
)

p_deer_plot <- ggplot(
  pred_deer_plot,
  aes(x = x, y = predicted, colour = group, fill = group)
) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  scale_colour_manual(values = pal_pel) +
  scale_fill_manual(values = pal_pel) +
  labs(
    title = "Downy birch: mixture × red deer + small deer pellets (plot level)",
    x = "Proportion silver (plot)",
    y = "Predicted browsing damage",
    colour = "Deer pellets",
    fill   = "Deer pellets"
  ) +
  theme_pub()

p_deer_plot

#----------------------------#
# 4D) Deer-combined stand-level figure
#----------------------------#
pred_deer_stand <- ggpredict(
  m_downy_stand_deer,
  terms = c("prop_s_stand [0:1 by=0.05]", "deer_pellets_log_mean [0,0.7,1.4]")
) %>% as.data.frame()

pred_deer_stand$group <- factor(
  pred_deer_stand$group,
  labels = c("Low", "Medium", "High")
)

p_deer_stand <- ggplot(
  pred_deer_stand,
  aes(x = x, y = predicted, colour = group, fill = group)
) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  scale_colour_manual(values = pal_pel) +
  scale_fill_manual(values = pal_pel) +
  labs(
    title = "Downy birch: mixture × red deer + small deer pellets (stand level)",
    x = "Proportion silver (stand)",
    y = "Predicted browsing damage",
    colour = "Deer pellets",
    fill   = "Deer pellets"
  ) +
  theme_pub()

p_deer_stand

#----------------------------#
# 5) Combine figures
#----------------------------#
fig_moose <- p_moose_plot + p_moose_stand +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

fig_deer <- p_deer_plot + p_deer_stand +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

fig_moose
fig_deer

# Optional export
# ggsave("Fig_moose_pellet_models.png", fig_moose, width = 11, height = 4.8, dpi = 600)
# ggsave("Fig_deer_pellet_models.png", fig_deer, width = 11, height = 4.8, dpi = 600)

#----------------------------#
# 6) Publication-style Viewer tables
#----------------------------#
library(sjPlot)

tab_model(
  m_downy_plot_moose,
  m_downy_stand_moose,
  dv.labels = c(
    "Downy plot model: damage_N ~ prop_s_plot * moose_pellets_log + north_c + stem_count_c + year + (1 | area) + (1 | stand_id)",
    "Downy stand model: damage_N ~ prop_s_stand * moose_pellets_log_mean + mean_north_c + stem_count_c + year + (1 | area)"
  ),
  show.ci = TRUE,
  show.se = TRUE,
  show.stat = TRUE,
  show.p = TRUE,
  show.aic = TRUE,
  show.icc = TRUE,
  show.re.var = TRUE,
  show.r2 = TRUE,
  show.obs = TRUE,
  title = "Table. Moose pellet mixed-effects models of downy birch browsing damage"
)

tab_model(
  m_downy_plot_deer,
  m_downy_stand_deer,
  dv.labels = c(
    "Downy plot model: damage_N ~ prop_s_plot * deer_pellets_log + north_c + stem_count_c + year + (1 | area) + (1 | stand_id)",
    "Downy stand model: damage_N ~ prop_s_stand * deer_pellets_log_mean + mean_north_c + stem_count_c + year + (1 | area)"
  ),
  show.ci = TRUE,
  show.se = TRUE,
  show.stat = TRUE,
  show.p = TRUE,
  show.aic = TRUE,
  show.icc = TRUE,
  show.re.var = TRUE,
  show.r2 = TRUE,
  show.obs = TRUE,
  title = "Table. Red deer + small deer pellet mixed-effects models of downy birch browsing damage"
)

# Optional HTML export
tab_model(
  m_downy_plot_moose,
  m_downy_stand_moose,
  dv.labels = c("Downy plot model", "Downy stand model"),
  show.ci = TRUE,
  show.se = TRUE,
  show.stat = TRUE,
  show.p = TRUE,
  show.aic = TRUE,
  show.icc = TRUE,
  show.re.var = TRUE,
  show.r2 = TRUE,
  show.obs = TRUE,
  title = "Table. Moose pellet models",
  file = "moose_pellet_models.html"
)

tab_model(
  m_downy_plot_deer,
  m_downy_stand_deer,
  dv.labels = c("Downy plot model", "Downy stand model"),
  show.ci = TRUE,
  show.se = TRUE,
  show.stat = TRUE,
  show.p = TRUE,
  show.aic = TRUE,
  show.icc = TRUE,
  show.re.var = TRUE,
  show.r2 = TRUE,
  show.obs = TRUE,
  title = "Table. Deer pellet models",
  file = "deer_pellet_models.html"
))

###############################################################################
# PART E: DOES BIRCH AFFECT PINE DAMAGE?
# Tests:
#   1) Does total birch abundance affect pine damage?
#   2) Do downy and silver birch differ in their effect on pine damage?
#   3) Does this relationship vary with latitude?
#   4) Does the relationship differ in northern vs southern stands?
###############################################################################

#----------------------------#
# 0) Packages
#----------------------------#
library(dplyr)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ggplot2)
library(ggeffects)
library(patchwork)
library(sjPlot)

#----------------------------#
# 1) Prepare pine dataset
#    Assumes Birch_2425 already exists and has north_model
#----------------------------#
pine_df <- Birch_2425 %>%
  mutate(
    year = as.factor(year),
    stand_number = as.factor(stand_number),
    
    pine_stems = suppressWarnings(as.numeric(as.character(pine_stems))),
    pine_damage_prop = suppressWarnings(as.numeric(as.character(proportion_pine_damage_winter))),
    
    downy_total = suppressWarnings(as.numeric(as.character(downy_total))),
    silver_total = suppressWarnings(as.numeric(as.character(silver_total))),
    north_model = suppressWarnings(as.numeric(as.character(north_model)))
  ) %>%
  mutate(
    total_birch = downy_total + silver_total
  ) %>%
  filter(
    !is.na(area),
    !is.na(stand_number),
    !is.na(year),
    !is.na(north_model),
    !is.na(pine_stems),
    pine_stems > 0,
    !is.na(pine_damage_prop)
  )

# Quick checks
summary(pine_df[, c("pine_damage_prop", "pine_stems", "downy_total", "silver_total", "total_birch", "north_model")])
table(pine_df$year, useNA = "ifany")

#----------------------------#
# 2) Scale predictors
#----------------------------#
pine_df <- pine_df %>%
  mutate(
    north_c = as.numeric(scale(north_model, center = TRUE, scale = FALSE)),
    downy_c = as.numeric(scale(downy_total, center = TRUE, scale = FALSE)),
    silver_c = as.numeric(scale(silver_total, center = TRUE, scale = FALSE)),
    birch_total_c = as.numeric(scale(total_birch, center = TRUE, scale = FALSE))
  )

#----------------------------#
# 3) Define north vs south
#    Split at the median latitude
#----------------------------#
north_split <- median(pine_df$north_model, na.rm = TRUE)

pine_df <- pine_df %>%
  mutate(
    region = ifelse(north_model > north_split, "North", "South"),
    region = factor(region, levels = c("South", "North"))
  )

table(pine_df$region, useNA = "ifany")

#----------------------------#
# 4) Core models
#----------------------------#

# 4A) Total birch effect
m_pine_total <- lmer(
  pine_damage_prop ~ birch_total_c + north_c + year +
    (1 | area) + (1 | stand_number),
  data = pine_df,
  REML = TRUE
)

# 4B) Species-specific birch effects
m_pine_species <- lmer(
  pine_damage_prop ~ downy_c + silver_c + north_c + year +
    (1 | area) + (1 | stand_number),
  data = pine_df,
  REML = TRUE
)

# 4C) Total birch × latitude
m_pine_total_north <- lmer(
  pine_damage_prop ~ birch_total_c * north_c + year +
    (1 | area) + (1 | stand_number),
  data = pine_df,
  REML = TRUE
)

# 4D) Species-specific birch × latitude
m_pine_species_north <- lmer(
  pine_damage_prop ~ downy_c * north_c + silver_c * north_c + year +
    (1 | area) + (1 | stand_number),
  data = pine_df,
  REML = TRUE
)

#----------------------------#
# 5) Separate north and south models
#----------------------------#
pine_north <- pine_df %>% filter(region == "North")
pine_south <- pine_df %>% filter(region == "South")

m_pine_total_north_only <- lmer(
  pine_damage_prop ~ birch_total_c + year +
    (1 | area) + (1 | stand_number),
  data = pine_north,
  REML = TRUE
)

m_pine_total_south_only <- lmer(
  pine_damage_prop ~ birch_total_c + year +
    (1 | area) + (1 | stand_number),
  data = pine_south,
  REML = TRUE
)

m_pine_species_north_only <- lmer(
  pine_damage_prop ~ downy_c + silver_c + year +
    (1 | area) + (1 | stand_number),
  data = pine_north,
  REML = TRUE
)

m_pine_species_south_only <- lmer(
  pine_damage_prop ~ downy_c + silver_c + year +
    (1 | area) + (1 | stand_number),
  data = pine_south,
  REML = TRUE
)

#----------------------------#
# 6) Summaries and drop1 tests
#----------------------------#
summary(m_pine_total)
r.squaredGLMM(m_pine_total)
drop1(update(m_pine_total, REML = FALSE), test = "Chisq")

summary(m_pine_species)
r.squaredGLMM(m_pine_species)
drop1(update(m_pine_species, REML = FALSE), test = "Chisq")

summary(m_pine_total_north)
r.squaredGLMM(m_pine_total_north)
drop1(update(m_pine_total_north, REML = FALSE), test = "Chisq")

summary(m_pine_species_north)
r.squaredGLMM(m_pine_species_north)
drop1(update(m_pine_species_north, REML = FALSE), test = "Chisq")

summary(m_pine_total_north_only)
summary(m_pine_total_south_only)
summary(m_pine_species_north_only)
summary(m_pine_species_south_only)

#----------------------------#
# 7) Figures
#----------------------------#
theme_pub <- function() {
  theme_classic(base_size = 13) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 12),
      plot.title.position = "plot",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

# 7A) Raw total birch vs pine damage, north vs south
p_raw_total <- ggplot(
  pine_df,
  aes(x = total_birch, y = pine_damage_prop, colour = region)
) +
  geom_point(alpha = 0.30, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.0) +
  labs(
    title = "A  Raw relationship: total birch abundance vs pine damage",
    x = "Total birch stems",
    y = "Pine winter damage (%)"
  ) +
  theme_pub()

p_raw_total

# 7B) Raw downy vs pine damage
p_raw_downy <- ggplot(
  pine_df,
  aes(x = downy_total, y = pine_damage_prop, colour = region)
) +
  geom_point(alpha = 0.30, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.0) +
  labs(
    title = "B  Raw relationship: downy birch abundance vs pine damage",
    x = "Downy birch stems",
    y = "Pine winter damage (%)"
  ) +
  theme_pub()

p_raw_downy

# 7C) Raw silver vs pine damage
p_raw_silver <- ggplot(
  pine_df,
  aes(x = silver_total, y = pine_damage_prop, colour = region)
) +
  geom_point(alpha = 0.30, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.0) +
  labs(
    title = "C  Raw relationship: silver birch abundance vs pine damage",
    x = "Silver birch stems",
    y = "Pine winter damage (%)"
  ) +
  theme_pub()

p_raw_silver

# 7D) Modelled total birch effect
pred_total <- ggpredict(
  m_pine_total,
  terms = c("birch_total_c [all]")
) %>% as.data.frame()

p_model_total <- ggplot() +
  geom_point(
    data = pine_df,
    aes(x = birch_total_c, y = pine_damage_prop),
    alpha = 0.20, size = 1.4
  ) +
  geom_ribbon(
    data = pred_total,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.15
  ) +
  geom_line(
    data = pred_total,
    aes(x = x, y = predicted),
    linewidth = 1.2
  ) +
  labs(
    title = "D  Modelled total birch effect on pine damage",
    x = "Total birch abundance (centered)",
    y = "Predicted pine winter damage (%)"
  ) +
  theme_pub()

p_model_total

# 7E) Modelled species-specific birch effects
pred_species_downy <- ggpredict(
  m_pine_species,
  terms = c("downy_c [all]")
) %>% as.data.frame()

pred_species_silver <- ggpredict(
  m_pine_species,
  terms = c("silver_c [all]")
) %>% as.data.frame()

p_model_downy <- ggplot() +
  geom_point(
    data = pine_df,
    aes(x = downy_c, y = pine_damage_prop),
    alpha = 0.20, size = 1.4
  ) +
  geom_ribbon(
    data = pred_species_downy,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.15
  ) +
  geom_line(
    data = pred_species_downy,
    aes(x = x, y = predicted),
    linewidth = 1.2
  ) +
  labs(
    title = "E  Modelled downy birch effect on pine damage",
    x = "Downy birch abundance (centered)",
    y = "Predicted pine winter damage (%)"
  ) +
  theme_pub()

p_model_silver <- ggplot() +
  geom_point(
    data = pine_df,
    aes(x = silver_c, y = pine_damage_prop),
    alpha = 0.20, size = 1.4
  ) +
  geom_ribbon(
    data = pred_species_silver,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.15
  ) +
  geom_line(
    data = pred_species_silver,
    aes(x = x, y = predicted),
    linewidth = 1.2
  ) +
  labs(
    title = "F  Modelled silver birch effect on pine damage",
    x = "Silver birch abundance (centered)",
    y = "Predicted pine winter damage (%)"
  ) +
  theme_pub()

p_model_downy
p_model_silver

# 7F) Total birch × latitude model
pred_total_north <- ggpredict(
  m_pine_total_north,
  terms = c("birch_total_c [all]", "north_c [-0.5,0,0.5]")
) %>% as.data.frame()

p_total_north <- ggplot(
  pred_total_north,
  aes(x = x, y = predicted, colour = group, fill = group)
) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  labs(
    title = "G  Total birch effect on pine damage at different latitudes",
    x = "Total birch abundance (centered)",
    y = "Predicted pine winter damage (%)",
    colour = "Latitude level",
    fill = "Latitude level"
  ) +
  theme_pub()

p_total_north

# 7G) North vs south separate fits
p_region_total <- ggplot(
  pine_df,
  aes(x = birch_total_c, y = pine_damage_prop, colour = region)
) +
  geom_point(alpha = 0.25, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.0) +
  facet_wrap(~ region, nrow = 1) +
  labs(
    title = "H  Total birch effect on pine damage in northern vs southern stands",
    x = "Total birch abundance (centered)",
    y = "Pine winter damage (%)"
  ) +
  theme_pub()

p_region_total

# 7H) Combined figure
fig_pine_main <- p_model_total + p_total_north + p_region_total +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "top")

fig_pine_main

# Optional save
# ggsave("Fig_pine_birch_effects.png", fig_pine_main, width = 15, height = 5, dpi = 600)

#----------------------------#
# 8) Viewer tables with whole-model information
#----------------------------#
tab_model(
  m_pine_total,
  m_pine_species,
  m_pine_total_north,
  m_pine_species_north,
  dv.labels = c(
    "Total birch model",
    "Species-specific birch model",
    "Total birch × latitude model",
    "Species-specific birch × latitude model"
  ),
  show.ci = TRUE,
  show.se = TRUE,
  show.stat = TRUE,
  show.p = TRUE,
  show.aic = TRUE,
  show.icc = TRUE,
  show.re.var = TRUE,
  show.r2 = TRUE,
  show.obs = TRUE,
  title = "Table. Mixed-effects models of pine browsing damage in relation to birch abundance and latitude"
)

tab_model(
  m_pine_total_north_only,
  m_pine_total_south_only,
  m_pine_species_north_only,
  m_pine_species_south_only,
  dv.labels = c(
    "North stands: total birch model",
    "South stands: total birch model",
    "North stands: species-specific model",
    "South stands: species-specific model"
  ),
  show.ci = TRUE,
  show.se = TRUE,
  show.stat = TRUE,
  show.p = TRUE,
  show.aic = TRUE,
  show.icc = TRUE,
  show.re.var = TRUE,
  show.r2 = TRUE,
  show.obs = TRUE,
  title = "Table. North vs south mixed-effects models of pine browsing damage"
)

# Optional HTML export
tab_model(
  m_pine_total,
  m_pine_species,
  m_pine_total_north,
  m_pine_species_north,
  dv.labels = c(
    "Total birch model",
    "Species-specific birch model",
    "Total birch × latitude model",
    "Species-specific birch × latitude model"
  ),
  show.ci = TRUE,
  show.se = TRUE,
  show.stat = TRUE,
  show.p = TRUE,
  show.aic = TRUE,
  show.icc = TRUE,
  show.re.var = TRUE,
  show.r2 = TRUE,
  show.obs = TRUE,
  title = "Table. Mixed-effects models of pine browsing damage in relation to birch abundance and latitude",
  file = "pine_birch_models.html"
)

tab_model(
  m_pine_total_north_only,
  m_pine_total_south_only,
  m_pine_species_north_only,
  m_pine_species_south_only,
  dv.labels = c(
    "North stands: total birch model",
    "South stands: total birch model",
    "North stands: species-specific model",
    "South stands: species-specific model"
  ),
  show.ci = TRUE,
  show.se = TRUE,
  show.stat = TRUE,
  show.p = TRUE,
  show.aic = TRUE,
  show.icc = TRUE,
  show.re.var = TRUE,
  show.r2 = TRUE,
  show.obs = TRUE,
  title = "Table. North vs south mixed-effects models of pine browsing damage",
  file = "pine_birch_models_north_south.html"
)

