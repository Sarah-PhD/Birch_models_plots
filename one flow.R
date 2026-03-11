################################################################################
# ONE FLOW: Add pine columns (and pine damage proportion) into your existing
# long birch modelling dataframe (e.g., birch_plot_long_pel)
#
# Assumptions:
# - Pine columns exist in Birch_2425 (plot-level wide-ish source):
#     pine_stems
#     pine_winter_damage_stems
# - Your long birch df (birch_plot_long_pel) already has identifiers:
#     area, stand_number, year
#   and ideally a plot identifier (e.g., plot_id, plot, plot_number, etc.)
#
# Goal:
# - Add pine_stems, pine_winter_damage_stems
# - Add pine_damage_prop = pine_winter_damage_stems / pine_stems
################################################################################

library(dplyr)

# --- 0) Decide your join keys -------------------------------------------------
# Use the most specific unique plot identifier you have.
# Common options: plot_id, plot, plot_number, provyta_id, etc.

# Example (CHANGE if your column is named differently):
plot_key <- "plot_id"

# Check whether both dfs have it
has_plot_key_in_Birch2425 <- plot_key %in% names(Birch_2425)
has_plot_key_in_long      <- plot_key %in% names(birch_plot_long_pel)

# If you DON'T have a plot_id in both, you MUST NOT join plot-level pine to
# plot-level birch using only stand/year/area (it will replicate rows).
# In that case, aggregate pine to stand-year FIRST (I include that option below).

# --- 1) Create a pine table at the correct grain ------------------------------

if (has_plot_key_in_Birch2425 && has_plot_key_in_long) {
  
  # 1A) PLOT-LEVEL pine table (one row per plot)
  pine_plot_tbl <- Birch_2425 %>%
    transmute(
      area,
      stand_number,
      year,
      !!plot_key := .data[[plot_key]],
      pine_stems               = as.numeric(pine_stems),
      pine_winter_damage_stems = as.numeric(pine_winter_damage_stems),
      pine_damage_prop = ifelse(
        pine_stems > 0,
        pine_winter_damage_stems / pine_stems,
        NA_real_
      )
    )
  
  # Optional sanity check: damaged <= total
  stopifnot(all(pine_plot_tbl$pine_winter_damage_stems <= pine_plot_tbl$pine_stems, na.rm = TRUE))
  
  # --- 2) Join into your long df --------------------------------------------
  birch_plot_long_pel2 <- birch_plot_long_pel %>%
    left_join(
      pine_plot_tbl,
      by = c("area", "stand_number", "year", plot_key)
    )
  
} else {
  
  # 1B) STAND-YEAR pine table (safe if you lack plot_id in one/both dfs)
  pine_stand_tbl <- Birch_2425 %>%
    group_by(area, stand_number, year) %>%
    summarise(
      pine_stems               = sum(as.numeric(pine_stems), na.rm = TRUE),
      pine_winter_damage_stems = sum(as.numeric(pine_winter_damage_stems), na.rm = TRUE),
      pine_damage_prop = ifelse(
        pine_stems > 0,
        pine_winter_damage_stems / pine_stems,
        NA_real_
      ),
      .groups = "drop"
    )
  
  stopifnot(all(pine_stand_tbl$pine_winter_damage_stems <= pine_stand_tbl$pine_stems, na.rm = TRUE))
  
  # --- 2) Join into your long df --------------------------------------------
  # This will repeat the stand-year pine value across the 10 plots (and across
  # species rows), which is OK if that is what you intend.
  birch_plot_long_pel2 <- birch_plot_long_pel %>%
    left_join(
      pine_stand_tbl,
      by = c("area", "stand_number", "year")
    )
}

# --- 3) Quick confirmation ----------------------------------------------------
birch_plot_long_pel2 %>%
  summarise(
    n = n(),
    pine_missing = sum(is.na(pine_damage_prop)),
    pine_prop_min = min(pine_damage_prop, na.rm = TRUE),
    pine_prop_max = max(pine_damage_prop, na.rm = TRUE)
  )

with(birch_long_plot, any(duplicated(interaction(area, stand_number))))
# should be FALSE if each stand exists in only one area (it usually is FALSE)

with(birch_long_plot, any(duplicated(stand_number)))
# if TRUE, then stand_number repeats across areas -> your current random effect is wrong
tab <- with(birch_long_plot, tapply(area, stand_number, function(x) length(unique(x))))
summary(tab)
any(tab > 1)
