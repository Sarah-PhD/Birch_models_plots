################################################################################
# STREAMLINED WORKFLOW: Birch mixture × latitude × density × (optional) pellets
# Data: Birch_2425  (plot-level rows; 10 plots per stand; stands nested in area)
#
# Main questions:
#   1) Does proportion silver (prop_s) affect browsing severity differently for
#      Downy vs Silver?   -> species × prop_s interaction
#   2) Does this hold at plot level AND stand-year level?
#   3) Optional: Is the mixture effect conditional on herbivore pressure (pellets)?
#
# Key design decisions (from our previous chats):
#   - Use prop_s only (NOT prop_d; collinear)
#   - Random structure:
#       Plot-level:  (1|area) + (1|stand_number)
#       Stand-level: (1|area)
#   - Year as fixed effect (2024/2025)
#   - Keep code DRY: one builder for long plot data, one builder for long stand data
#   - Clean/standardise AREA names to avoid misspellings splitting groups
################################################################################
my_colors <- c("#FDE725FF", "#56C667FF", "#238A8DFF", "#3F4788FF")
#----------------------------#
# 0) Packages
#----------------------------#
library(dplyr)
library(tidyr)
library(stringr)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ggeffects)
library(ggplot2)

#----------------------------#
# 1) AREA name cleaning (IMPORTANT)
#----------------------------#
# 1A) Minimal standardisation (always do this)
standardise_area <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_squish() %>%
    str_replace_all("[\u2010-\u2015]", "-") %>%   # hyphen variants -> "-"
    str_replace_all("[_]", "") %>%               # remove underscores if they occur
    str_replace_all("\\s+", "")                  # remove spaces (keeps names consistent)
}

# 1B) Explicit corrections for known misspellings.
# Edit this list as you discover new variants in table(Birch_2425$area).
area_map <- c(
  # known from earlier: Lycksele was misspelled as Lyksele
  "Lyksele"     = "Lycksele",
  "LYKSELE"     = "Lycksele",
  "lycksele"    = "Lycksele",
  
  # common formatting variants (examples — keep/edit what applies to your data)
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
      # apply explicit recode on BOTH raw and standardised variants:
      area = dplyr::recode(area_raw, !!!area_map, .default = area_raw),
      area = dplyr::recode(area, !!!area_map, .default = area),
      # final standardisation to remove whitespace/formatting differences
      area = standardise_area(area)
    ) %>%
    select(-area_raw, -area_std)
}

# Apply cleaning
Birch_2425 <- Birch_2425 %>%
  fix_area_names(area)

# Quick audit: look for suspicious near-duplicates after cleaning
# (If you see duplicates that should be merged, add them to area_map above.)
sort(unique(Birch_2425$area))

#----------------------------#
# 2) Helper: robust numeric conversion
#----------------------------#
num <- function(x) suppressWarnings(as.numeric(x))

#----------------------------#
# 3) Build PLOT-level long dataset (core variables + optional pellets)
#----------------------------#
build_plot_long <- function(df, include_pellets = FALSE) {
  
  # Optional pellet index (robust: NA -> 0, then log1p)
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
  
  wide <- df %>%
    transmute(
      year         = as.factor(year),
      area         = area,
      stand_number = as.factor(stand_number),
      north        = num(north),
      north_c      = as.numeric(scale(north, center = TRUE, scale = FALSE)),
      
      downy_total      = num(downy_total),
      silver_total     = num(silver_total),
      downy_damage_N   = num(downy_damage_N),
      silver_damage_N  = num(silver_damage_N),
      
      # optional pellets:
      total_pellets_log = if (include_pellets) num(total_pellets_log) else NULL,
      pellet_present    = if (include_pellets) num(pellet_present) else NULL
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
      any_of(c("total_pellets_log", "pellet_present"))
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
# 4) Build STAND-year long dataset (aggregated from plot-wide)
#----------------------------#
build_stand_long <- function(plot_long, include_pellets = FALSE) {
  
  # First rebuild a stand-year wide summary from plot_long
  stand_wide <- plot_long %>%
    group_by(year, area, stand_number) %>%
    summarise(
      mean_north_c = mean(north_c, na.rm = TRUE),
      
      # species-specific means at stand-year
      downy_damage_mean  = mean(damage_N[species == "Downy"],  na.rm = TRUE),
      silver_damage_mean = mean(damage_N[species == "Silver"], na.rm = TRUE),
      
      downy_stems_mean   = mean(stem_count[species == "Downy"],  na.rm = TRUE),
      silver_stems_mean  = mean(stem_count[species == "Silver"], na.rm = TRUE),
      
      pellets_log_mean   = if (include_pellets) mean(total_pellets_log, na.rm = TRUE) else NA_real_,
      pellet_present_any = if (include_pellets) max(pellet_present, na.rm = TRUE) else NA_real_,
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
      any_of(c("pellets_log_mean", "pellet_present_any"))
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

#----------------------------#
# 5) CORE MODELS (NO HEIGHT)
#----------------------------#

# 5A) Plot-level combined species model (main test = species:prop_s_plot)
m_plot <- lmer(
  damage_N ~ species * prop_s_plot +
    species * north_c +
    species * stem_count +
    year +
    (1 | area) +
    (1 | stand_number),
  data = birch_plot_long,
  REML = TRUE
)

# 5B) Stand-level combined species model (main test = species:prop_s_stand)
m_stand <- lmer(
  damage_N ~ species * prop_s_stand +
    species * mean_north_c +
    species * stem_count +
    year +
    (1 | area),
  data = birch_stand_long,
  REML = TRUE
)

# Summaries + R2 + ML drop tests
summary(m_plot)
r.squaredGLMM(m_plot)
drop1(update(m_plot, REML = FALSE), test = "Chisq")

summary(m_stand)
r.squaredGLMM(m_stand)
drop1(update(m_stand, REML = FALSE), test = "Chisq")

#----------------------------#
# 6) OPTIONAL: Pellets × mixture (Downy-only; plot + stand)
#----------------------------#
# Use this if pellets are central to the mechanism.
birch_downy_plot_pel <- birch_plot_long_pel %>% filter(species == "Downy")
birch_downy_stand_pel <- birch_stand_long_pel %>% filter(species == "Downy")

m_downy_plot_pel <- lmer(
  damage_N ~ prop_s_plot * total_pellets_log +
    north_c + stem_count + year +
    (1 | area) + (1 | stand_number),
  data = birch_downy_plot_pel,
  REML = TRUE
)

m_downy_stand_pel <- lmer(
  damage_N ~ prop_s_stand * pellets_log_mean +
    mean_north_c + stem_count + year +
    (1 | area),
  data = birch_downy_stand_pel,
  REML = TRUE
)

summary(m_downy_plot_pel)
drop1(update(m_downy_plot_pel, REML = FALSE), test = "Chisq")

summary(m_downy_stand_pel)
drop1(update(m_downy_stand_pel, REML = FALSE), test = "Chisq")

#________________________________
#6a)
#Model for species not individual

Sp <- lmer(damage_N ~ species * prop_s_plot * total_pellets_log
                            + north_c + stem_count + year
                            + (1 | area) + (1 | stand_number))
#---------------------------------
#6b) Graphs for the proportion plus pellets


library(ggeffects)

pred_plot <- ggpredict(
  m_downy_plot_pel,
  terms = c("total_pellets_log", "prop_s_plot [0, 0.5, 1]")
)

library(ggplot2)

ggplot(pred_plot, aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, colour = NA) +
  scale_colour_manual(values = c("#3BA64C", "#66C26F", "#1F7F35")) +
  scale_fill_manual(values = c("#3BA64C", "#66C26F", "#1F7F35")) +
  labs(
    x = "Log pellet density",
    y = "Predicted Downy browsing severity",
    colour = "Silver proportion",
    fill = "Silver proportion"
  ) +
  theme_classic()

#Stand levek

pred_stand <- ggpredict(
  m_downy_stand_pel,
  terms = c("pellets_log_mean", "prop_s_stand [0, 0.5, 1]")
)

#----------------------------#
# 7) CLEAN PREDICTION PLOTS (publication-ready)
#----------------------------#
theme_pub <- function() {
  theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 11)
    )
}

plot_prop_effect <- function(model, prop_term, title_x) {
  ggpredict(model, terms = c(paste0(prop_term, " [0:1 by=0.05]"), "species")) %>%
    ggplot(aes(x = x, y = predicted, colour = group, fill = group)) +
    geom_line(linewidth = 1.1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
    labs(x = title_x, y = "Predicted browsing damage", colour = "Species", fill = "Species") +
    theme_pub()
}

p_plot_prop  <- plot_prop_effect(m_plot,  "prop_s_plot",  "Proportion silver (plot)")
p_stand_prop <- plot_prop_effect(m_stand, "prop_s_stand", "Proportion silver (stand)")

p_plot_prop
p_stand_prop

# Optional: Downy-only pellets interaction plot at meaningful pellet levels
pel_levels <- c(0, 0.7, 1.4)  # log1p: 0=no pellets; ~0.7≈1; ~1.4≈3 (roughly)
pred_downy_pel <- ggpredict(
  m_downy_plot_pel,
  terms = c("prop_s_plot [0:1 by=0.05]",
            paste0("total_pellets_log [", paste(pel_levels, collapse=","), "]"))
)

p_downy_pel <- ggplot(pred_downy_pel, aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  labs(
    x = "Proportion silver (plot)",
    y = "Predicted DOWNY browsing damage",
    colour = "Total pellets\n(log1p)",
    fill   = "Total pellets\n(log1p)"
  ) +
  theme_pub()

p_downy_pel

################################################################################
# What to report (key terms):
#   Plot-level mixture test:  speciesSilver:prop_s_plot
#   Stand-level mixture test: speciesSilver:prop_s_stand
#   Downy pellets mechanism:  prop_s_plot:total_pellets_log (and stand equivalent)
################################################################################


birch_stand_long %>%
  filter(prop_s_stand >= 0.75) %>%
  summarise(n_stand_years = n())


birch_stand_long %>%
  filter(prop_s_stand >= 0.75) %>%
  distinct(stand_number) %>%
  summarise(n_stands = n())

birch_stand_long %>%
  group_by(stand_number) %>%
  summarise(min_prop = min(prop_s_stand, na.rm = TRUE)) %>%
  filter(min_prop >= 0.75) %>%
  summarise(n_stands = n())

birch_stand_long %>%
  mutate(class = case_when(
    prop_s_stand < 0 ~ "Downy-dominated (<25% silver)",
    prop_s_stand < 1 ~ "Mixed (25–75% silver)",
    TRUE ~ "Silver-dominated (≥75%)"
  )) %>%
  count(class)

##########################################################################
##########    GRAPHS    ###########################################
############################################################################
#Damage vs Proportion Silver (stand level
library(ggeffects)
library(ggplot2)

# Model-based predictions
pred_prop <- ggpredict(
  m_stand,
  terms = c("prop_s_stand [0:1 by=0.05]", "species")
)

p_prop_stand <- ggplot() +
  
  # Raw data
  geom_point(
    data = birch_stand_long,
    aes(x = prop_s_stand, y = damage_N, colour = species),
    alpha = 0.4
  ) +
  
  # Model trend lines
  geom_line(
    data = pred_prop,
    aes(x = x, y = predicted, colour = group),
    linewidth = 1.2
  ) +
  
  geom_ribbon(
    data = pred_prop,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.15,
    colour = NA
  ) +
  
  scale_colour_manual(values = c("Downy" = "#3BA64C", "Silver" = "#4C5C8C")) +
  scale_fill_manual(values = c("Downy" = "#3BA64C", "Silver" = "#4C5C8C")) +
  
  labs(
    x = "Proportion silver (stand)",
    y = "Browsing damage (stand-level mean)",
    colour = "Species",
    fill = "Species"
  ) +
  
  theme_classic(base_size = 13)

p_prop_stand

#Damage vs North Coordinate (stand level)

pred_north <- ggpredict(
  m_stand,
  terms = c("mean_north_c [all]", "species")
)

p_north_stand <- ggplot() +
  
  geom_point(
    data = birch_stand_long,
    aes(x = mean_north_c, y = damage_N, colour = species),
    alpha = 0.4
  ) +
  
  geom_line(
    data = pred_north,
    aes(x = x, y = predicted, colour = group),
    linewidth = 1.2
  ) +
  
  geom_ribbon(
    data = pred_north,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.15,
    colour = NA
  ) +
  
  scale_colour_manual(values = c("Downy" = "#3BA64C", "Silver" = "#4C5C8C")) +
  scale_fill_manual(values = c("Downy" = "#3BA64C", "Silver" = "#4C5C8C")) +
  
  labs(
    x = "North coordinate (centered, stand mean)",
    y = "Browsing damage (stand-level mean)",
    colour = "Species",
    fill = "Species"
  ) +
  
  theme_classic(base_size = 13)

p_north_stand

ggplot(birch_stand_long,
       aes(x = prop_s_stand,
           y = damage_N,
           colour = species)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic()


#=========================================================
# 2-panel, journal-ready figure (stand level)
# Panel A: Damage ~ proportion silver (stand)
# Panel B: Damage ~ north (stand mean, centered)
#
# Uses your stand-level mixed model object: m_stand
# Data: birch_stand_long
#
# Requirements you asked for:
#  - 2 panels in one figure
#  - same colour palette as before
#  - clean journal style
#  - optional "mechanism figure" emphasising Silver × North (below)
#=========================================================

library(ggplot2)
library(ggeffects)
library(patchwork)

#----- Your palette (same as earlier) -----
pal <- c("Downy" = "#3BA64C", "Silver" = "#4C5C8C")

#----- Consistent journal theme -----
theme_journal <- function() {
  theme_classic(base_size = 13) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 12),
      plot.title.position = "plot",
      panel.grid = element_blank()
    )
}

#=========================================================
# Panel A: Prop silver
#=========================================================
pred_prop <- ggpredict(
  m_stand,
  terms = c("prop_s_stand [0:1 by=0.05]", "species")
) |> as.data.frame()

pA <- ggplot() +
  # raw points
  geom_point(
    data = birch_stand_long,
    aes(x = prop_s_stand, y = damage_N, colour = species),
    alpha = 0.35, size = 1.8
  ) +
  # model ribbon + line
  geom_ribbon(
    data = pred_prop,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.14, colour = NA
  ) +
  geom_line(
    data = pred_prop,
    aes(x = x, y = predicted, colour = group),
    linewidth = 1.25
  ) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(
    title = "A  Stand composition effect",
    x = "Proportion silver (stand)",
    y = "Browsing damage (stand mean)"
  ) +
  theme_journal()

#=========================================================
# Panel B: North coordinate (mean_north_c)
#=========================================================
pred_north <- ggpredict(
  m_stand,
  terms = c("mean_north_c [all]", "species")
) |> as.data.frame()

pB <- ggplot() +
  geom_point(
    data = birch_stand_long,
    aes(x = mean_north_c, y = damage_N, colour = species),
    alpha = 0.35, size = 1.8
  ) +
  geom_ribbon(
    data = pred_north,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.14, colour = NA
  ) +
  geom_line(
    data = pred_north,
    aes(x = x, y = predicted, colour = group),
    linewidth = 1.25
  ) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title = "B  Latitude effect (species-specific)",
    x = "North (centered; stand mean)",
    y = "Browsing damage (stand mean)"
  ) +
  theme_journal()

#=========================================================
# Combine into one 2-panel figure
#=========================================================
p_two_panel <- pA + pB +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

p_two_panel

# Optional: export high-res
# ggsave("Fig_stand_prop_and_north.png", p_two_panel, width = 11, height = 4.8, dpi = 600)
# ggsave("Fig_stand_prop_and_north.pdf", p_two_panel, width = 11, height = 4.8)

#=========================================================
# Add rugs (data density) AND optional faceting by year
# Works for BOTH panels (prop_s and north)
#=========================================================

library(ggplot2)
library(ggeffects)
library(patchwork)

pal <- c("Downy" = "#3BA64C", "Silver" = "#4C5C8C")

theme_journal <- function() {
  theme_classic(base_size = 13) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 12),
      plot.title.position = "plot",
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

#=========================================================
# A) Damage ~ prop_s_stand (with rugs + optional facet by year)
#=========================================================
pred_prop <- ggpredict(
  m_stand,
  terms = c("prop_s_stand [0:1 by=0.05]", "species", "year")
) |> as.data.frame()

# If you prefer NOT faceted, you can switch to terms = c("prop_s_stand ...", "species")
# and set facet_year <- FALSE below.

facet_year <- TRUE  # <--- toggle

pA <- ggplot() +
  # raw points
  geom_point(
    data = birch_stand_long,
    aes(x = prop_s_stand, y = damage_N, colour = species),
    alpha = 0.30, size = 1.7
  ) +
  # x-density rugs (raw data)
  geom_rug(
    data = birch_stand_long,
    aes(x = prop_s_stand, colour = species),
    sides = "b", alpha = 0.18, linewidth = 0.6
  ) +
  # model ribbon + line
  geom_ribbon(
    data = pred_prop,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.14, colour = NA
  ) +
  geom_line(
    data = pred_prop,
    aes(x = x, y = predicted, colour = group),
    linewidth = 1.2
  ) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(
    title = "A  Stand composition effect",
    x = "Proportion silver (stand)",
    y = "Browsing damage (stand mean)"
  ) +
  theme_journal()

if (facet_year) {
  pA <- pA + facet_wrap(~ facet, nrow = 1)  # facet column from ggpredict is called 'facet'
}

#=========================================================
# B) Damage ~ mean_north_c (with rugs + optional facet by year)
#=========================================================
pred_north <- ggpredict(
  m_stand,
  terms = c("mean_north_c [all]", "species", "year")
) |> as.data.frame()

pB <- ggplot() +
  geom_point(
    data = birch_stand_long,
    aes(x = mean_north_c, y = damage_N, colour = species),
    alpha = 0.30, size = 1.7
  ) +
  geom_rug(
    data = birch_stand_long,
    aes(x = mean_north_c, colour = species),
    sides = "b", alpha = 0.18, linewidth = 0.6
  ) +
  geom_ribbon(
    data = pred_north,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.14, colour = NA
  ) +
  geom_line(
    data = pred_north,
    aes(x = x, y = predicted, colour = group),
    linewidth = 1.2
  ) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title = "B  Latitude effect (species-specific)",
    x = "North (centered; stand mean)",
    y = "Browsing damage (stand mean)"
  ) +
  theme_journal()

if (facet_year) {
  pB <- pB + facet_wrap(~ facet, nrow = 1)
}

#=========================================================
# Combine into one figure
#=========================================================
p_two_panel <- pA + pB +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

p_two_panel

# Optional export:
# ggsave("Fig_stand_prop_and_north_faceted_rug.png",
#        p_two_panel, width = 12.5, height = 5.2, dpi = 600)
# ggsave("Fig_stand_prop_and_north_faceted_rug.pdf",


#RAW GRAPH

ggplot(birch_stand_long,
       aes(x = mean_north_c,
           y = damage_N,
           colour = species)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic()
#        p_two_panel, width = 12.5, height = 5.2)

############
#SIDE BY SIDE 
#=========================================================
# SIDE-BY-SIDE: Raw vs Mixed-model predictions (stand level)
# For BOTH x-axes:
#   (1) Proportion silver (prop_s_stand)
#   (2) North coordinate (mean_north_c)
#
# Left = RAW (points + lm smooth)
# Right = MODEL (points + mixed-model marginal effect line + CI)
# Facet by year (2024 vs 2025) in BOTH, clean journal style
#=========================================================

library(ggplot2)
library(ggeffects)
library(patchwork)
library(dplyr)

pal <- c("Downy" = "#3BA64C", "Silver" = "#4C5C8C")

theme_journal <- function() {
  theme_classic(base_size = 13) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 11)
    )
}

# Ensure year is a factor for faceting
birch_stand_long <- birch_stand_long %>% mutate(year = as.factor(year))

#=========================================================
# 1) PROPORTION: Raw vs Model
#=========================================================

# RAW plot (lm smooth per species, per year)
p_prop_raw <- ggplot(birch_stand_long,
                     aes(x = prop_s_stand, y = damage_N, colour = species)) +
  geom_point(alpha = 0.30, size = 1.7) +
  geom_rug(aes(x = prop_s_stand), sides = "b", alpha = 0.15, linewidth = 0.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.0) +
  facet_wrap(~ year, nrow = 1) +
  scale_colour_manual(values = pal) +
  labs(
    title = "RAW: points + lm smooth",
    x = "Proportion silver (stand)",
    y = "Browsing damage (stand mean)"
  ) +
  theme_journal()

# MODEL plot (mixed-model marginal predictions)
pred_prop <- ggpredict(
  m_stand,
  terms = c("prop_s_stand [0:1 by=0.05]", "species", "year")
) %>% as.data.frame()

p_prop_model <- ggplot() +
  geom_point(
    data = birch_stand_long,
    aes(x = prop_s_stand, y = damage_N, colour = species),
    alpha = 0.25, size = 1.6
  ) +
  geom_rug(
    data = birch_stand_long,
    aes(x = prop_s_stand, colour = species),
    sides = "b", alpha = 0.12, linewidth = 0.6
  ) +
  geom_ribbon(
    data = pred_prop,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.14, colour = NA
  ) +
  geom_line(
    data = pred_prop,
    aes(x = x, y = predicted, colour = group),
    linewidth = 1.2
  ) +
  facet_wrap(~ facet, nrow = 1) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title = "MODEL: mixed-model marginal effect (±95% CI)",
    x = "Proportion silver (stand)",
    y = "Browsing damage (stand mean)"
  ) +
  theme_journal()

fig_prop_raw_vs_model <- p_prop_raw + p_prop_model +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

fig_prop_raw_vs_model

#=========================================================
# 2) NORTH: Raw vs Model
#=========================================================

# RAW plot (lm smooth per species, per year)
p_north_raw <- ggplot(birch_stand_long,
                      aes(x = mean_north_c, y = damage_N, colour = species)) +
  geom_point(alpha = 0.30, size = 1.7) +
  geom_rug(aes(x = mean_north_c), sides = "b", alpha = 0.15, linewidth = 0.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.0) +
  facet_wrap(~ year, nrow = 1) +
  scale_colour_manual(values = pal) +
  labs(
    title = "RAW: points + lm smooth",
    x = "North (centered; stand mean)",
    y = "Browsing damage (stand mean)"
  ) +
  theme_journal()

# MODEL plot (mixed-model marginal predictions)
pred_north <- ggpredict(
  m_stand,
  terms = c("mean_north_c [all]", "species", "year")
) %>% as.data.frame()

p_north_model <- ggplot() +
  geom_point(
    data = birch_stand_long,
    aes(x = mean_north_c, y = damage_N, colour = species),
    alpha = 0.25, size = 1.6
  ) +
  geom_rug(
    data = birch_stand_long,
    aes(x = mean_north_c, colour = species),
    sides = "b", alpha = 0.12, linewidth = 0.6
  ) +
  geom_ribbon(
    data = pred_north,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.14, colour = NA
  ) +
  geom_line(
    data = pred_north,
    aes(x = x, y = predicted, colour = group),
    linewidth = 1.2
  ) +
  facet_wrap(~ facet, nrow = 1) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title = "MODEL: mixed-model marginal effect (±95% CI)",
    x = "North (centered; stand mean)",
    y = "Browsing damage (stand mean)"
  ) +
  theme_journal()

fig_north_raw_vs_model <- p_north_raw + p_north_model +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

fig_north_raw_vs_model

# Optional exports:
# ggsave("Fig_prop_raw_vs_model.png", fig_prop_raw_vs_model, width = 13, height = 4.8, dpi = 600)
# ggsave("Fig_north_raw_vs_model.png", fig_north_raw_vs_model, width = 13, height = 4.8, dpi = 600)
# ggsave("Fig_prop_raw_vs_model.pdf", fig_prop_raw_vs_model, width = 13, height = 4.8)
# ggsave("Fig_north_raw_vs_model.pdf", fig_north_raw_vs_model, width = 13, height = 4.8)

