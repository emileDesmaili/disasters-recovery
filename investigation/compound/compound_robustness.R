# compound_robustness.R
# Additional robustness checks for compound vs isolated shock analysis
# 1. Continuous intensity interaction (wind speed × had_recent_shock)
# 2. Cumulative exposure measure (number of shocks in past 3 years)
# 3. Dropping always-compound countries (CHN, PHL, AUS, MEX, JPN)
# 4. Period heterogeneity within compound/isolated

library(tidyverse)
library(fixest)
library(haven)
library(countrycode)

source("../../emileRegs.R")
setFixest_notes(FALSE)

# ------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------
if (!file.exists("treated_panel.rds")) {
  source("compound_diagnostics.R")
}
d_treated <- readRDS("treated_panel.rds")

d_treated <- d_treated %>%
  group_by(countrycode) %>%
  mutate(
    shock_lag1  = lag(maxwind_95, 1, default = 0L),
    shock_lag2  = lag(maxwind_95, 2, default = 0L),
    shock_lag3  = lag(maxwind_95, 3, default = 0L),
    compound_strict = as.integer(maxwind_95 == 1 & shock_lag1 == 1),
    isolated_strict = as.integer(maxwind_95 == 1 & shock_lag1 == 0),
    had_recent_shock = as.factor(shock_lag1),
    # Cumulative: number of shocks in t-1 to t-3
    recent_shock_count = shock_lag1 + shock_lag2 + shock_lag3,
    recent_shock_cat = case_when(
      recent_shock_count == 0 ~ "0 recent",
      recent_shock_count == 1 ~ "1 recent",
      recent_shock_count >= 2 ~ "2+ recent"
    ),
    recent_shock_cat = factor(recent_shock_cat,
                              levels = c("0 recent", "1 recent", "2+ recent")),
    region = countrycode(countrycode, "iso3c", "region"),
    region = ifelse(is.na(region), "Other", region)
  ) %>%
  ungroup()

outcome  <- "loggdp"
horizon  <- 10
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year
controls <- "l(gdp_diff,1:2) + l(maxwind_95,1:2)"

# ==================================================================
# ROBUSTNESS 1: Continuous wind speed × recent shock interaction
# ==================================================================
cat("=== Robustness 1: Continuous intensity × recent shock ===\n")

# Standardize wind speed among treated
d_treated <- d_treated %>%
  mutate(
    maxwind_sd = maxwind / sd(maxwind[maxwind > 0], na.rm = TRUE)
  )

irf_r1 <- lp_panel_inter(
  data = d_treated, outcome = outcome, main_var = "maxwind_sd",
  interact_var = "had_recent_shock",
  controls = "l(gdp_diff,1:2) + l(maxwind_sd,1:2)",
  horizon = horizon,
  fe = "countrycode[year] + countrycode[year2] + year",
  panel_id = panel_id, vcov_formula = vcov_fm
) %>%
  mutate(
    category = case_when(
      str_detect(category, "^0") ~ "Isolated",
      str_detect(category, "^1") ~ "Compound",
      TRUE ~ category
    )
  )
irf_r1$spec <- "Continuous Intensity × Recent Shock"

# ==================================================================
# ROBUSTNESS 2: Cumulative exposure (0, 1, 2+ shocks in past 3 yrs)
# ==================================================================
cat("=== Robustness 2: Cumulative exposure categories ===\n")

irf_r2 <- lp_panel_inter(
  data = d_treated, outcome = outcome, main_var = "maxwind_95",
  interact_var = "recent_shock_cat",
  controls = controls,
  horizon = horizon,
  fe = "countrycode[year] + countrycode[year2] + year",
  panel_id = panel_id, vcov_formula = vcov_fm
) %>%
  mutate(
    category = str_replace_all(category, "_", " ")
  )
irf_r2$spec <- "Cumulative Exposure (past 3 years)"

# ==================================================================
# ROBUSTNESS 3: Drop always-compound countries
# ==================================================================
cat("=== Robustness 3: Drop always-compound countries ===\n")

# Countries where >90% of shocks are compound
country_shock <- readRDS("country_shock_summary.rds")
always_compound <- country_shock %>%
  filter(pct_compound >= 90, n_shocks >= 10) %>%
  pull(countrycode)

cat(sprintf("Dropping %d always-compound countries: %s\n",
            length(always_compound), paste(always_compound, collapse = ", ")))

d_drop <- d_treated %>%
  filter(!(countrycode %in% always_compound))

cat(sprintf("Remaining: %d countries, %d obs\n",
            n_distinct(d_drop$countrycode), nrow(d_drop)))
cat(sprintf("Shocks remaining: %d (compound: %d, isolated: %d)\n",
            sum(d_drop$maxwind_95),
            sum(d_drop$compound_strict),
            sum(d_drop$isolated_strict)))

irf_r3 <- lp_panel_inter(
  data = d_drop, outcome = outcome, main_var = "maxwind_95",
  interact_var = "had_recent_shock",
  controls = controls,
  horizon = horizon,
  fe = "countrycode[year] + countrycode[year2] + year",
  panel_id = panel_id, vcov_formula = vcov_fm
) %>%
  mutate(
    category = case_when(
      str_detect(category, "^0") ~ "Isolated",
      str_detect(category, "^1") ~ "Compound",
      TRUE ~ category
    )
  )
irf_r3$spec <- "Drop Always-Compound Countries"

# ==================================================================
# ROBUSTNESS 4: Triple interaction — compound × period
# ==================================================================
cat("=== Robustness 4: Compound × Period interaction ===\n")

d_treated <- d_treated %>%
  mutate(
    compound_period = paste0(
      ifelse(shock_lag1 == 1, "Compound", "Isolated"),
      " × ",
      period
    ),
    compound_period = factor(compound_period,
                             levels = c("Isolated × Pre-1990",
                                        "Isolated × Post-1990",
                                        "Compound × Pre-1990",
                                        "Compound × Post-1990"))
  )

irf_r4 <- lp_panel_inter(
  data = d_treated, outcome = outcome, main_var = "maxwind_95",
  interact_var = "compound_period",
  controls = controls,
  horizon = horizon,
  fe = "countrycode[year] + countrycode[year2] + year",
  panel_id = panel_id, vcov_formula = vcov_fm
) %>%
  mutate(
    category = str_replace_all(category, "_", "-") %>%
      str_replace("Pre-1990", "Pre-1990") %>%
      str_replace("Post-1990", "Post-1990")
  )
irf_r4$spec <- "Compound × Period"

# ==================================================================
# ROBUSTNESS 5: Country FE only (no trends, for maximum power)
# ==================================================================
cat("=== Robustness 5: Country FE only (no trends) ===\n")

irf_r5 <- lp_panel_inter(
  data = d_treated, outcome = outcome, main_var = "maxwind_95",
  interact_var = "had_recent_shock",
  controls = controls,
  horizon = horizon,
  fe = "countrycode + year",
  panel_id = panel_id, vcov_formula = vcov_fm
) %>%
  mutate(
    category = case_when(
      str_detect(category, "^0") ~ "Isolated",
      str_detect(category, "^1") ~ "Compound",
      TRUE ~ category
    )
  )
irf_r5$spec <- "Country + Year FE (no trends)"

# ------------------------------------------------------------------
# Print summary at h=5
# ------------------------------------------------------------------
cat("\n========== ROBUSTNESS RESULTS AT h=5 ==========\n")

all_rob <- bind_rows(irf_r1, irf_r2, irf_r3, irf_r5)
for (s in unique(all_rob$spec)) {
  cat(sprintf("\n--- %s ---\n", s))
  h5 <- all_rob %>% filter(spec == s, horizon == 5)
  for (cat_name in unique(h5$category)) {
    row <- h5 %>% filter(category == cat_name)
    cat(sprintf("  %-35s: %.3f (SE=%.3f) [%.3f, %.3f]\n",
                cat_name, row$irf_mean, row$se, row$irf_down, row$irf_up))
  }
}

cat("\n--- Compound × Period ---\n")
h5_r4 <- irf_r4 %>% filter(horizon == 5)
for (cat_name in unique(h5_r4$category)) {
  row <- h5_r4 %>% filter(category == cat_name)
  if (nrow(row) > 0)
    cat(sprintf("  %-40s: %.3f (SE=%.3f)\n",
                cat_name, row$irf_mean, row$se))
}

# Save results
saveRDS(list(
  r1_continuous  = irf_r1,
  r2_cumulative  = irf_r2,
  r3_drop_always = irf_r3,
  r4_triple      = irf_r4,
  r5_no_trends   = irf_r5
), "robustness_results.rds")
cat("\nSaved: robustness_results.rds\n")
cat("Done.\n")
