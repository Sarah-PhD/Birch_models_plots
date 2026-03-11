library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpmisc)
install.packages("see")
#----------------------------#
# A) Make plot-level long df
#----------------------------#
birch_long_plot <- Birch_2425 %>%
  select(plot_id, stand_number, area, year, north,
         downy_total, silver_total) %>%
  pivot_longer(
    cols = c(downy_total, silver_total),
    names_to = "species",
    values_to = "stems"
  ) %>%
  mutate(
    species = recode(species,
                     downy_total = "Downy birch",
                     silver_total = "Silver birch"),
    stems = as.numeric(stems),
    north = as.numeric(north)
  ) %>%
  filter(!is.na(north), !is.na(stems))

birch_long_plot <- birch_long_plot %>%
  mutate(
    north_c = scale(north, center = TRUE, scale = TRUE)[,1]
  )

# Optional: quick look
summary(birch_long_plot$stems)
table(birch_long_plot$species)

#----------------------------#
# B) Linear model:
#    does slope differ by species?
#----------------------------#
m_plot <- lm(stems ~ species * north, data = birch_long_plot)
summary(m_plot)

# This ANOVA table gives the p-value for species:north (your key test)
anova(m_plot)


anova(m_plot)["species:north", ]

#------------------------------------#
#  glmer
# -----------------------------------#

##############################
#CTUAL MODEL

m_nb <- glmer.nb(
  stems ~ species * north_c + factor(year) +
    (1 | area/stand_number),
  data = birch_long_plot,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
summary(m_nb)
car::Anova(m_nb, type = 3)
##############       PLOT: STEMS VS NORTH             ####################
library(dplyr)
library(ggplot2)

#-----------------------------------------------------------
# 1) Stand-level dataset (one point per stand per species)
#    Averaged across both years + plots
#-----------------------------------------------------------
birch_stand <- birch_long_plot %>%
  group_by(area, stand_number, species) %>%
  summarise(
    latitude = mean(north, na.rm = TRUE),     # real latitude (°N)
    stems_mean = mean(stems, na.rm = TRUE),   # mean stems per plot, summarised to stand level
    .groups = "drop"
  )

#-----------------------------------------------------------
# 2) Journal-style plot (two lines)
#-----------------------------------------------------------
p_journal_lat <- ggplot(birch_stand, aes(x = latitude, y = stems_mean, color = species)) +
  geom_point(size = 1.4) +                               # smaller points (no alpha)
  geom_smooth(method = "lm", se = FALSE, linewidth = 2.2) +  # thicker lines
  coord_cartesian(ylim = c(0, 120), expand = FALSE) +     # <- adjust 120 as you like
  scale_x_continuous(breaks = seq(57, 66, by = 1)) +
  scale_y_continuous(breaks = seq(0, 120, by = 20)) +
  scale_color_manual(values = my_colors[c(2,4)]) +
  labs(
    x = "Latitude (°N)",
    y = "Mean stems per plot (stand summary)",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "top",
    legend.text = element_text(size = 11),
    plot.margin = margin(8, 10, 8, 10)
  )

p_journal_lat

###############################################################################
#                      AREA REPRESENTATION ON GRAPH                           #
###############################################################################
library(dplyr)
library(ggplot2)

# area label positions
area_positions <- birch_stand %>%
  group_by(area) %>%
  summarise(latitude = mean(latitude, na.rm = TRUE), .groups = "drop") %>%
  arrange(latitude)

# pick ~6 areas spread across the latitude gradient
n_labels <- 6
label_areas <- area_positions %>%
  slice(round(seq(1, n(), length.out = n_labels)))

p <- ggplot(birch_stand, aes(x = latitude, y = stems_mean, color = species)) +
  geom_point(size = 1.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 2.2) +
  coord_cartesian(ylim = c(0, 120), clip = "off", expand = FALSE) +
  scale_x_continuous(breaks = seq(57, 66, by = 1)) +
  scale_y_continuous(breaks = seq(0, 120, by = 20)) +
  scale_color_manual(values = my_colors[c(2,4)]) +
  labs(x = "Latitude (°N)", y = "Mean stems per stand", color = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.margin = margin(10, 10, 45, 10) # bottom room for labels
  ) +
  geom_text(
    data = label_areas,
    aes(x = latitude, y = -10, label = area),
    inherit.aes = FALSE,
    angle = 45,
    size = 4,
    vjust =1
  )

p

###############################################################################
####                   Stand level                                        #####
###############################################################################

library(dplyr)

birch_stand <- birch_long_plot %>%
  group_by(area, stand_number, year, species) %>%
  summarise(
    stems_mean = mean(stems, na.rm = TRUE),
    latitude   = mean(north, na.rm = TRUE),
    .groups = "drop"
  )

summary(birch_stand$stems_mean)


####
# FIT MODEL
####
library(lme4)

birch_stand <- birch_stand %>%
  mutate(
    north_c = as.numeric(scale(latitude)),
    year_f  = factor(year)
  )

m_nb_stand <- glmer.nb(
  stems_mean ~ species * north_c + year_f +
    (1 | area),
  data = birch_stand,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(m_nb_stand)
car::Anova(m_nb_stand, type = 3)


###########################################################################
#               STand level with sums 
############################################################################

library(dplyr)

# Stand-level MEAN stems per plot (continuous; good for Gaussian)
birch_stand_mean <- birch_long_plot %>%
  group_by(area, stand_number, year, species) %>%
  summarise(
    stems_mean = mean(stems, na.rm = TRUE),
    latitude   = mean(north, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    north_c = as.numeric(scale(latitude)),
    year_f  = factor(year)
  )

# Stand-level SUM stems across plots (integer counts; good for NB)
birch_stand_sum <- birch_long_plot %>%
  group_by(area, stand_number, year, species) %>%
  summarise(
    stems_sum = sum(stems, na.rm = TRUE),
    latitude  = mean(north, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    north_c = as.numeric(scale(latitude)),
    year_f  = factor(year)
  )


####
# Model comparison
####
# Gausian transformed
library(lme4)

m_gauss_mean <- lmer(
  log1p(stems_mean) ~ species * north_c + year_f + (1 | area),
  data = birch_stand_mean,
  REML = FALSE
)

summary(m_gauss_mean)
car::Anova(m_gauss_mean, type = 3)


#gausian, non transformed (less reccommended)

m_gauss_mean_raw <- lmer(
  stems_mean ~ species * north_c + year_f + (1 | area),
  data = birch_stand_mean,
  REML = FALSE
)
AIC(m_gauss_mean, m_gauss_mean_raw)

# Negative binomial GLMM (on sum counts)
library(lme4)

m_nb_sum <- glmer.nb(
  stems_sum ~ species * north_c + year_f + (1 | area),
  data = birch_stand_sum,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(m_nb_sum)
car::Anova(m_nb_sum, type = 3)


#Comparison: They cannot be compared using AIC as they are on count or mean
performance::check_model(m_gauss_mean)

m_pois_sum <- glmer(
  stems_sum ~ species * north_c + year_f + (1 | area),
  family = poisson,
  data = birch_stand_sum
)
performance::check_overdispersion(m_pois_sum)  # should show why NB is needed

### Didnt work so ttying this way instead to compare
par(mfrow = c(1,2))
plot(fitted(m_gauss_mean), resid(m_gauss_mean), xlab="Fitted", ylab="Residuals")
abline(h = 0, lty = 2)
qqnorm(resid(m_gauss_mean)); qqline(resid(m_gauss_mean))
par(mfrow = c(1,1))
library(lme4)
library(dplyr)

# If not already created:
birch_stand_sum <- birch_long_plot %>%
  group_by(area, stand_number, year, species) %>%
  summarise(
    stems_sum = sum(stems, na.rm = TRUE),
    latitude  = mean(north, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    north_c = as.numeric(scale(latitude)),
    year_f  = factor(year)
  )

m_pois_sum <- glmer(
  stems_sum ~ species * north_c + year_f + (1 | area),
  family = poisson,
  data = birch_stand_sum,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

performance::check_overdispersion(m_pois_sum)

m_nb_sum <- glmer.nb(
  stems_sum ~ species * north_c + year_f + (1 | area),
  data = birch_stand_sum,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(m_nb_sum)
car::Anova(m_nb_sum, type = 3)

birch_long_plot %>%
  distinct(area, stand_number, plot_id) %>%   # use your actual plot ID column
  count(area, stand_number) %>%
  count(n)
#############################################################################
# PLOTS
########################################################################

library(dplyr)
library(ggplot2)
library(ggeffects)

# ----------------------------
# 0) Ensure factors + scaled north are consistent
# ----------------------------
birch_long_plot <- birch_long_plot %>%
  mutate(
    year_f  = factor(year),
    north_c = as.numeric(scale(north))   # make sure plot-level scaling is truly from latitude
  )

# Stand-level SUM dataset (integer counts; best for NB)
birch_stand_sum <- birch_long_plot %>%
  group_by(area, stand_number, year, species) %>%
  summarise(
    stems_sum = sum(stems, na.rm = TRUE),
    latitude  = mean(north, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    year_f  = factor(year),
    north_c = as.numeric(scale(latitude))
  )

# ----------------------------
# 1) Fit models (or reuse if already fitted)
# ----------------------------
m_nb_plot <- glmer.nb(
  stems ~ species * north_c + year_f + (1 | area) + (1 | stand_number),
  data = birch_long_plot,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

m_nb_sum <- glmer.nb(
  stems_sum ~ species * north_c + year_f + (1 | area),
  data = birch_stand_sum,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# ----------------------------
# 2) Predictions for BOTH years
# ----------------------------
pred_plot <- ggpredict(m_nb_plot, terms = c("north_c [all]", "species", "year_f")) %>%
  as.data.frame() %>%
  mutate(scale = "Plot level (NB; stems per plot)")

pred_sum <- ggpredict(m_nb_sum, terms = c("north_c [all]", "species", "year_f")) %>%
  as.data.frame() %>%
  mutate(scale = "Stand level (NB; total stems per stand)")

# ----------------------------
# 3) Convert north_c (scaled) back to latitude (°N)
#    NOTE: scaling differs between plot model and stand model
# ----------------------------
lat_mu_plot <- mean(birch_long_plot$north, na.rm = TRUE)
lat_sd_plot <- sd(birch_long_plot$north, na.rm = TRUE)

lat_mu_stand <- mean(birch_stand_sum$latitude, na.rm = TRUE)
lat_sd_stand <- sd(birch_stand_sum$latitude, na.rm = TRUE)

pred_plot <- pred_plot %>%
  mutate(latitude = lat_mu_plot + (x * lat_sd_plot))

pred_sum <- pred_sum %>%
  mutate(latitude = lat_mu_stand + (x * lat_sd_stand))

# Combine
pred_all <- bind_rows(pred_plot, pred_sum) %>%
  mutate(
    year = as.character(facet),     # ggpredict stores year factor in 'facet'
    species = group                # ggpredict stores species in 'group'
  )

# ----------------------------
# 4) Publication-style plot
# ----------------------------
p_pred <- ggplot(pred_all, aes(x = latitude, y = predicted, color = species, fill = species)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.18, color = NA) +
  geom_line(linewidth = 2.2) +
  facet_grid(scale ~ year, scales = "free_y") +
  scale_x_continuous(breaks = seq(57, 66, by = 1)) +
  scale_color_manual(values = my_colors[c(2,4)]) +
  scale_fill_manual(values = my_colors[c(2,4)]) +
  labs(
    x = "Latitude (°N)",
    y = "Predicted stems",
    color = NULL,
    fill  = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

p_pred

