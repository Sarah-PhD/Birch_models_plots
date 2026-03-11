################################################################################
# ONE FLOW: Sensitivity analysis (exclude Åtvidaberg, Malå, Barksätter)
# Author: Sarah Gore
# Purpose:
#   Re-run your key birch + pine models after removing three areas that look
#   potentially influential in the year-by-area plot:
#     Åtvidaberg, Malå, Barksätter
#
# Assumptions (objects already exist from your workflow):
#   - Birch_2425         : plot-level wide data (cleaned; 2024/2025 only)
#   - birch_plot_long    : long birch plot-level df (species rows; damage_N etc.)
#   - birch_stand_long   : long birch stand-year df (species rows; damage_N etc.)
#   - pine_stand_binom   : stand-year pine df with undamaged + total_pellets_log
#   - pine_plot          : plot-level pine df with undamaged + total_pellets_log
#
# Output:
#   - Models refit on filtered data
#   - Side-by-side fixed-effect comparisons vs original models
################################################################################

#----------------------------#
# 0) Packages
#----------------------------#
library(dplyr)
library(lme4)
library(lmerTest)
library(car)

#----------------------------#
# 1) Define areas to exclude
#----------------------------#
areas_exclude <- c("Åtvidaberg", "Malå", "Barksätter")

# Quick audit: confirm spelling matches your data
cat("\n--- Areas present in Birch_2425 ---\n")
print(sort(unique(Birch_2425$area)))

cat("\n--- Excluding these areas ---\n")
print(areas_exclude)

#----------------------------#
# 2) Filter ALL relevant datasets (keep naming consistent)
#----------------------------#

# Birch (plot-level long)
birch_plot_long_sens <- birch_plot_long %>%
  filter(!area %in% areas_exclude)

# Birch (stand-year long)
birch_stand_long_sens <- birch_stand_long %>%
  filter(!area %in% areas_exclude)

# Pine (stand-year binomial df)
pine_stand_binom_sens <- pine_stand_binom %>%
  filter(!area %in% areas_exclude)

# Pine (plot-level pine df)
pine_plot_sens <- pine_plot %>%
  filter(!area %in% areas_exclude)

# Sanity checks
cat("\n--- N rows after exclusion ---\n")
print(tibble(
  birch_plot_long  = nrow(birch_plot_long_sens),
  birch_stand_long = nrow(birch_stand_long_sens),
  pine_stand_binom = nrow(pine_stand_binom_sens),
  pine_plot        = nrow(pine_plot_sens)
))

cat("\n--- Areas remaining (pine stand) ---\n")
print(sort(unique(pine_stand_binom_sens$area)))

#----------------------------#
# 3) Refit BIRCH models (your core mixed models)
#    NOTE: These match your earlier "core" models.
#----------------------------#

# 3A) Birch PLOT-level combined species model
#     Key tests typically:
#       species:prop_s_plot   (mixture effect by species)
#       species:north_c       (latitude effect by species)
m_birch_plot_sens <- lmer(
  damage_N ~ species * prop_s_plot +
    species * north_c +
    species * stem_count +
    year +
    (1 | area) + (1 | stand_number),
  data = birch_plot_long_sens,
  REML = TRUE
)

cat("\n================ BIRCH: Plot-level (SENS) ================\n")
print(summary(m_birch_plot_sens))
cat("\n--- Type III tests (BIRCH plot SENS) ---\n")
print(car::Anova(m_birch_plot_sens, type = 3))

# 3B) Birch STAND-year combined species model
#     Uses stand-year variables you created:
#       prop_s_stand, mean_north_c, stem_count, year
m_birch_stand_sens <- lmer(
  damage_N ~ species * prop_s_stand +
    species * mean_north_c +
    species * stem_count +
    year +
    (1 | area),
  data = birch_stand_long_sens,
  REML = TRUE
)

cat("\n================ BIRCH: Stand-year (SENS) ================\n")
print(summary(m_birch_stand_sens))
cat("\n--- Type III tests (BIRCH stand SENS) ---\n")
print(car::Anova(m_birch_stand_sens, type = 3))

#----------------------------#
# 4) Refit PINE models (binomial GLMMs)
#    Response = cbind(damaged, undamaged) = proportion with correct weighting
#----------------------------#

# Ensure factors (safe)
pine_stand_binom_sens <- pine_stand_binom_sens %>%
  mutate(
    year = factor(year),
    stand_number = factor(stand_number)
  )

pine_plot_sens <- pine_plot_sens %>%
  mutate(
    year = factor(year),
    stand_number = factor(stand_number)
  )

# 4A) Pine STAND-year model (full; includes pellets + birch predictors)
m_pine_stand_sens <- glmer(
  cbind(pine_winter_damage_stems, undamaged) ~
    scale(birch_density) +
    birch_prop_s +
    total_pellets_log +
    north_c +
    year +
    (1 | stand_number),
  family = binomial,
  data = pine_stand_binom_sens,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5))
)

cat("\n================ PINE: Stand-year (SENS) ================\n")
print(summary(m_pine_stand_sens))
cat("\n--- Type III tests (PINE stand SENS) ---\n")
print(car::Anova(m_pine_stand_sens, type = 3))

# 4B) Pine PLOT-level model (full; includes pellets + birch predictors)
m_pine_plot_sens <- glmer(
  cbind(pine_winter_damage_stems, undamaged) ~
    scale(birch_density_plot) +
    birch_prop_s_plot +
    total_pellets_log +
    north_c +
    year +
    (1 | stand_number),
  family = binomial,
  data = pine_plot_sens,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5))
)

cat("\n================ PINE: Plot-level (SENS) ================\n")
print(summary(m_pine_plot_sens))
cat("\n--- Type III tests (PINE plot SENS) ---\n")
print(car::Anova(m_pine_plot_sens, type = 3))

#----------------------------#
# 5) Compare FIXED effects: original vs sensitivity
#    (Assumes your original models exist in environment with these names)
#      - m_birch_plot
#      - m_birch_stand
#      - m_pine_binom_full    (stand-year pine full)
#      - m_pine_plot_full2    (plot-level pine full)
# If any are named differently, rename below.
#----------------------------#

cat("\n================ FIXED EFFECT COMPARISONS ================\n")

# Helper: compact fixed-effect table
fe_tbl <- function(mod) {
  co <- summary(mod)$coefficients
  out <- data.frame(term = rownames(co), co, row.names = NULL)
  names(out) <- sub("Pr\\(>\\|z\\|\\)", "p", names(out))
  names(out) <- sub("Pr\\(>\\|t\\|\\)", "p", names(out))
  out
}

# ---- BIRCH plot ----
if (exists("m_birch_plot")) {
  cat("\n--- Birch plot: ORIGINAL vs SENS ---\n")
  print(list(
    ORIGINAL = fe_tbl(m_birch_plot),
    SENS     = fe_tbl(m_birch_plot_sens)
  ))
} else {
  cat("\nNOTE: m_birch_plot not found; skipping Birch plot comparison.\n")
}

# ---- BIRCH stand ----
if (exists("m_stand")) {
  cat("\n--- Birch stand: ORIGINAL (m_stand) vs SENS ---\n")
  print(list(
    ORIGINAL = fe_tbl(m_stand),
    SENS     = fe_tbl(m_birch_stand_sens)
  ))
} else if (exists("m_birch_stand")) {
  cat("\n--- Birch stand: ORIGINAL vs SENS ---\n")
  print(list(
    ORIGINAL = fe_tbl(m_birch_stand),
    SENS     = fe_tbl(m_birch_stand_sens)
  ))
} else {
  cat("\nNOTE: No original birch stand model found (m_stand or m_birch_stand); skipping.\n")
}

# ---- PINE stand ----
if (exists("m_pine_binom_full")) {
  cat("\n--- Pine stand: ORIGINAL vs SENS ---\n")
  print(list(
    ORIGINAL = fe_tbl(m_pine_binom_full),
    SENS     = fe_tbl(m_pine_stand_sens)
  ))
} else {
  cat("\nNOTE: m_pine_binom_full not found; skipping Pine stand comparison.\n")
}

# ---- PINE plot ----
if (exists("m_pine_plot_full2")) {
  cat("\n--- Pine plot: ORIGINAL vs SENS ---\n")
  print(list(
    ORIGINAL = fe_tbl(m_pine_plot_full2),
    SENS     = fe_tbl(m_pine_plot_sens)
  ))
} else {
  cat("\nNOTE: m_pine_plot_full2 not found; skipping Pine plot comparison.\n")
}

cat("\nDONE: sensitivity models fitted + comparisons printed.\n")
################################################################################

################################################################################
# ONE FLOW: Visualise RAW damage data between years (pine + birch)
# Author: Sarah Gore
#
# Goal:
#   Show raw distributions (not model predictions) for 2024 vs 2025.
#
# Assumes you already have:
#   - pine_plot        (plot-level pine df with pine_winter_damage_stems, pine_stems, year, area, stand_number)
#   - birch_plot_long  (plot-level long birch df with damage_N, species, year, area, stand_number)
#
# Output:
#   Figure 1: Pine raw plot-level proportions by year (jitter + box/violin)
#   Figure 2: Birch raw plot-level damage by year, faceted by species
#   Figure 3: Area-level RAW points between years (pine + birch) in facets
################################################################################

library(dplyr)
library(ggplot2)
library(tidyr)

#----------------------------#
# 0) Prep: compute raw pine proportion per plot
#----------------------------#
pine_plot_raw <- pine_plot %>%
  mutate(
    year = factor(year),
    area = factor(area),
    stand_number = factor(stand_number),
    pine_prop = pine_winter_damage_stems / pine_stems
  ) %>%
  filter(
    is.finite(pine_prop),
    pine_stems > 0
  )

# For birch, ensure factors
birch_plot_raw <- birch_plot_long %>%
  mutate(
    year = factor(year),
    area = factor(area),
    stand_number = factor(stand_number)
  ) %>%
  filter(!is.na(damage_N))

#----------------------------#
# 1) FIGURE 1: Pine raw damage proportion by year (distribution + raw points)
#   - violin shows full distribution
#   - box shows median/IQR
#   - jitter shows every observation (raw)
#----------------------------#
p_pine_year_raw <- ggplot(pine_plot_raw, aes(x = year, y = pine_prop)) +
  geom_violin(trim = FALSE, alpha = 0.25) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.45) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.25, size = 1.1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Pine winter browsing (RAW, plot level)",
    x = NULL,
    y = "Proportion winter-damaged pine stems"
  ) +
  theme_classic(base_size = 13)

p_pine_year_raw

#----------------------------#
# 2) FIGURE 2: Birch raw damage by year, split by species
#   - This shows the year shift AND the species difference in raw values
#----------------------------#
p_birch_year_raw <- ggplot(birch_plot_raw, aes(x = year, y = damage_N)) +
  geom_violin(trim = FALSE, alpha = 0.25) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.45) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.18, size = 1.0) +
  facet_wrap(~ species, nrow = 1, scales = "free_y") +
  labs(
    title = "Birch browsing severity (RAW, plot level)",
    x = NULL,
    y = "Damage severity (damage_N)"
  ) +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

p_birch_year_raw

#----------------------------#
# 3) FIGURE 3: Area-level RAW points (means per area-year)
#   - This is still “raw summaries”: no modelling, just aggregation.
#   - Cleaner than dual axes: use facets for pine vs birch.
#----------------------------#

# Pine area-year mean prop (raw)
pine_area_year_raw <- pine_plot_raw %>%
  group_by(area, year) %>%
  summarise(pine_prop = mean(pine_prop, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "Pine: winter damage proportion",
         value = pine_prop)

# Birch area-year mean (raw) — average across species within area-year
birch_area_year_raw <- birch_plot_raw %>%
  group_by(area, year) %>%
  summarise(birch_damage = mean(damage_N, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "Birch: browsing severity (damage_N)",
         value = birch_damage)

area_year_raw <- bind_rows(
  pine_area_year_raw %>% select(area, year, metric, value),
  birch_area_year_raw %>% select(area, year, metric, value)
)

p_area_year_raw <- ggplot(area_year_raw,
                          aes(x = reorder(area, value, FUN = mean),
                              y = value,
                              colour = year,
                              group = year)) +
  geom_point(position = position_dodge(width = 0.45), size = 2.2) +
  geom_line(aes(group = interaction(area, metric)),
            colour = "grey60", linewidth = 0.4, alpha = 0.6) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(
    title = "RAW area-year summaries (no model)",
    x = "Area",
    y = NULL,
    colour = "Year"
  ) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

p_area_year_raw
################################################################################
birch_long_plot %>%
  count(area, stand_number, plot_id) %>%
  count(n)

