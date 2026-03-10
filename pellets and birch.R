################################################################################
# Plotting pellet × silver proportion interaction
# (A) Downy-only models (plot + stand) + plots
# (B) Recommended full model including species (3-way interaction) + plot
################################################################################

# Packages
library(dplyr)
library(lme4)
library(lmerTest)
library(ggeffects)
library(ggplot2)
library(patchwork)

# --- Your palette (from your figure) ---
COL_DOWNY  <- "#3BA64C"  # green
COL_SILVER <- "#4C5C8C"  # blue

# Optional: helper for 3 shades of a base color (no extra packages)
shade3 <- function(col) {
  grDevices::colorRampPalette(c(col, "white"))(4)[1:3]
}

################################################################################
# 1) Fit the Downy-only models (as you already do)
################################################################################

birch_downy_plot_pel  <- birch_plot_long_pel  %>% filter(species == "Downy")
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

# Summaries / tests (optional)
summary(m_downy_plot_pel)
drop1(update(m_downy_plot_pel, REML = FALSE), test = "Chisq")

summary(m_downy_stand_pel)
drop1(update(m_downy_stand_pel, REML = FALSE), test = "Chisq")

################################################################################
# 2) Plot the interaction: PREDICTED lines across pellets at prop_s = 0, 0.5, 1
################################################################################

# ---- Plot-level predictions ----
pred_downy_plot <- ggpredict(
  m_downy_plot_pel,
  terms = c("total_pellets_log", "prop_s_plot [0, 0.5, 1]")
)

# Three green shades for the three prop_s lines
greens <- shade3(COL_DOWNY)
names(greens) <- unique(pred_downy_plot$group)

p_downy_plot <- ggplot(pred_downy_plot, aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = greens) +
  scale_fill_manual(values = greens) +
  labs(
    title  = "Plot level (Downy only): pellets × silver proportion",
    x      = "Total pellets (log scale)",
    y      = "Predicted browsing severity (damage_N)",
    colour = "prop_s_plot",
    fill   = "prop_s_plot"
  ) +
  theme_classic()

# ---- Stand-level predictions ----
pred_downy_stand <- ggpredict(
  m_downy_stand_pel,
  terms = c("pellets_log_mean", "prop_s_stand [0, 0.5, 1]")
)

greens2 <- shade3(COL_DOWNY)
names(greens2) <- unique(pred_downy_stand$group)

p_downy_stand <- ggplot(pred_downy_stand, aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = greens2) +
  scale_fill_manual(values = greens2) +
  labs(
    title  = "Stand level (Downy only): pellets × silver proportion",
    x      = "Mean pellets (log scale)",
    y      = "Predicted browsing severity (damage_N)",
    colour = "prop_s_stand",
    fill   = "prop_s_stand"
  ) +
  theme_classic()

# ---- Combine into one figure (top = plot, bottom = stand) ----
(p_downy_plot / p_downy_stand) + plot_layout(heights = c(1, 1))

################################################################################
# 3) Recommended: ONE model including species (tests if the mechanism differs)
################################################################################

# This assumes birch_plot_long_pel has BOTH species with columns:
# damage_N, prop_s_plot, total_pellets_log, north_c, stem_count, year, area, stand_number, species

m_species_plot_pel <- lmer(
  damage_N ~ species * prop_s_plot * total_pellets_log +
    species * north_c +
    species * stem_count +
    year +
    (1 | area) + (1 | stand_number),
  data = birch_plot_long_pel,
  REML = TRUE
)

summary(m_species_plot_pel)

# Likelihood ratio style (term-deletion) test for the 3-way interaction
drop1(update(m_species_plot_pel, REML = FALSE), test = "Chisq")

################################################################################
# 4) Plot the species model: pellets × prop_s, faceted by species
################################################################################

pred_species_plot <- ggpredict(
  m_species_plot_pel,
  terms = c("total_pellets_log", "prop_s_plot [0, 0.5, 1]", "species")
)

# Set colours by species; line type by prop_s group (clean + readable)
p_species <- ggplot(pred_species_plot,
                    aes(x = x, y = predicted,
                        colour = facet, linetype = group, fill = facet)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~facet) +
  scale_colour_manual(values = c("Downy" = COL_DOWNY, "Silver" = COL_SILVER)) +
  scale_fill_manual(values   = c("Downy" = COL_DOWNY, "Silver" = COL_SILVER)) +
  labs(
    title    = "Plot level: species-specific pellets × silver proportion effect",
    x        = "Total pellets (log scale)",
    y        = "Predicted browsing severity (damage_N)",
    colour   = "Species",
    fill     = "Species",
    linetype = "prop_s_plot"
  ) +
  theme_classic()

p_species

################################################################################
# Notes:
# - If your species labels are "Downy birch" / "Silver birch" instead of "Downy"/"Silver",
#   change the names in scale_*_manual accordingly.
################################################################################

birch_silver_plot_pel <- birch_plot_long_pel %>%
  filter(species == "Silver")

m_silver_plot_pel <- lmer(
  damage_N ~ prop_s_plot * total_pellets_log +
    north_c + stem_count + year +
    (1 | area) + (1 | stand_number),
  data = birch_silver_plot_pel,
  REML = TRUE
)
summary(m_silver_plot_pel)

drop1(update(m_silver_plot_pel, REML = FALSE), test = "Chisq")




birch_silver_stand_pel <- birch_stand_long_pel %>%
  filter(species == "Silver")

m_silver_stand_pel <- lmer(
  damage_N ~ prop_s_stand * pellets_log_mean +
    mean_north_c + stem_count + year +
    (1 | area),
  data = birch_silver_stand_pel,
  REML = TRUE
)

summary(m_silver_stand_pel)

drop1(update(m_silver_stand_pel, REML = FALSE), test = "Chisq")



## Downy
birch_plot_long_pel <- birch_plot_long_pel %>%
  mutate(prop_d_plot = 1 - prop_s_plot)

birch_silver_plot_pel <- birch_plot_long_pel %>%
  filter(species == "Silver")

m_silver_downyprop <- lmer(
  damage_N ~ prop_d_plot * total_pellets_log +
    north_c + stem_count + year +
    (1 | area) + (1 | stand_number),
  data = birch_silver_plot_pel,
  REML = TRUE
)
summary(m_silver_downyprop)

drop1(update(m_silver_downyprop, REML = FALSE), test = "Chisq")


# Packages
library(ggeffects)
library(ggplot2)

# Your palette (from earlier)
COL_SILVER <- "#4C5C8C"

# 1) Predictions: pellets gradient, at prop_d = 0, 0.5, 1
pred_silver_facil <- ggpredict(
  m_silver_downyprop,
  terms = c("total_pellets_log", "prop_d_plot [0, 0.5, 1]")
)

# 2) Plot
ggplot(pred_silver_facil, aes(x = x, y = predicted, linetype = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.18) +
  geom_line(linewidth = 1.2, colour = COL_SILVER) +
  labs(
    title    = "Silver birch browsing shows facilitation by Downy proportion under high herbivore load",
    subtitle = "Predicted silver damage across pellet density at low/medium/high Downy proportion",
    x        = "Pellet density (log scale)",
    y        = "Predicted Silver browsing severity (damage_N)",
    linetype = "Downy proportion (prop_d_plot)"
  ) +
  theme_classic()

ggplot() +
  geom_point(
    data = birch_silver_plot_pel,
    aes(x = total_pellets_log, y = damage_N),
    alpha = 0.20,
    size = 1
  ) +
  geom_ribbon(
    data = pred_silver_facil,
    aes(x = x, ymin = conf.low, ymax = conf.high, group = group),
    alpha = 0.18
  ) +
  geom_line(
    data = pred_silver_facil,
    aes(x = x, y = predicted, linetype = group, group = group),
    linewidth = 1.2,
    colour = COL_SILVER
  ) +
  labs(
    title    = "Silver browsing increases with herbivore load most strongly when Downy is abundant",
    x        = "Pellet density (log scale)",
    y        = "Silver browsing severity (damage_N)",
    linetype = "Downy proportion (prop_d_plot)"
  ) +
  theme_classic()

#################################################################################
#                           PINE                                              #
##############################################################################

