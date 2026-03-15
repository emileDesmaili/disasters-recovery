# compound_estimation.R
# Estimate differential effects of compound vs isolated shocks
# Multiple specifications and robustness checks

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

# Also load full panel for comparison
pwt  <- read_stata("../../raw_data/pwt_clean.dta")
tcs  <- read_stata("../../raw_data/ibtracs_clean.dta")
data_full <- pwt %>% left_join(tcs, by = c("year", "countrycode"))

knot_to_ms <- 0.514444
data_full <- data_full %>%
  mutate(
    maxwind  = replace_na(max_ann_wind_i_nob, 0) * knot_to_ms,
    loggdp   = ln_real_gdp_usd_pc,
    gdp_diff = growth_real_gdp_usd_pc,
    year2    = year^2,
    period   = ifelse(year <= 1990, "Pre-1990", "Post-1990")
  ) %>%
  filter(year >= 1969) %>%
  mutate(
    maxwind_p95 = quantile(maxwind[year >= 1970 & year <= 2014], 0.95, na.rm = TRUE),
    maxwind_95 = as.integer(maxwind_p95 > 0 & maxwind >= maxwind_p95),
    region = countrycode(countrycode, "iso3c", "region"),
    region = ifelse(is.na(region), "Other", region)
  )

# Ensure compound variables exist
d_treated <- d_treated %>%
  group_by(countrycode) %>%
  mutate(
    shock_lag1  = lag(maxwind_95, 1, default = 0L),
    shock_lag2  = lag(maxwind_95, 2, default = 0L),
    compound_strict = as.integer(maxwind_95 == 1 & shock_lag1 == 1),
    compound_broad  = as.integer(maxwind_95 == 1 & (shock_lag1 == 1 | shock_lag2 == 1)),
    isolated_strict = as.integer(maxwind_95 == 1 & shock_lag1 == 0),
    isolated_broad  = as.integer(maxwind_95 == 1 & shock_lag1 == 0 & shock_lag2 == 0),
    region = countrycode(countrycode, "iso3c", "region"),
    region = ifelse(is.na(region), "Other", region)
  ) %>%
  ungroup()

cat(sprintf("Treated panel: %d obs, %d countries\n",
            nrow(d_treated), n_distinct(d_treated$countrycode)))
cat(sprintf("Compound (strict): %d, Isolated: %d\n",
            sum(d_treated$compound_strict), sum(d_treated$isolated_strict)))

# ------------------------------------------------------------------
# Estimation helpers
# ------------------------------------------------------------------
outcome  <- "loggdp"
horizon  <- 10
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

# Custom LP function for two separate shock variables
lp_two_shocks <- function(data, outcome, shock1, shock2,
                           controls = NULL, horizon = 10,
                           fe = "countrycode[year] + countrycode[year2] + year",
                           panel_id = c("countrycode", "year"),
                           vcov_formula = DK ~ year) {
  irf_list <- vector("list", horizon + 1)
  rhs_controls <- if (!is.null(controls)) paste0(" + ", controls) else ""

  for (h in 0:horizon) {
    fml <- as.formula(paste0(
      "f(", outcome, ", ", h, ") - l(", outcome, ", 1) ~ ",
      shock1, " + ", shock2,
      rhs_controls, " | ", fe
    ))

    mod <- feols(fml, data = data, panel.id = panel_id, vcov = vcov_formula)

    for (var in c(shock1, shock2)) {
      beta <- coef(mod)[var]
      se   <- sqrt(vcov(mod)[var, var])
      irf_list[[length(irf_list) + 1]] <- data.frame(
        horizon  = h,
        term     = var,
        irf_mean = beta,
        se       = se,
        irf_down = beta - 1.96 * se,
        irf_up   = beta + 1.96 * se,
        stringsAsFactors = FALSE
      )
    }

    # Test equality
    diff <- coef(mod)[shock1] - coef(mod)[shock2]
    v1 <- vcov(mod)[shock1, shock1]
    v2 <- vcov(mod)[shock2, shock2]
    cov12 <- vcov(mod)[shock1, shock2]
    diff_se <- sqrt(v1 + v2 - 2 * cov12)
    irf_list[[length(irf_list) + 1]] <- data.frame(
      horizon  = h,
      term     = "Difference",
      irf_mean = diff,
      se       = diff_se,
      irf_down = diff - 1.96 * diff_se,
      irf_up   = diff + 1.96 * diff_se,
      stringsAsFactors = FALSE
    )
  }
  bind_rows(irf_list)
}

# ==================================================================
# SPECIFICATION 1: Separate compound & isolated shock indicators
# (ever-treated sample only, year FE)
# ==================================================================
cat("\n=== Spec 1: Separate shock indicators, treated sample, year FE ===\n")

controls_1 <- "l(gdp_diff,1:2) + l(maxwind_95,1:2)"

irf_spec1 <- lp_two_shocks(
  data = d_treated, outcome = outcome,
  shock1 = "isolated_strict", shock2 = "compound_strict",
  controls = controls_1, horizon = horizon,
  fe = "countrycode[year] + countrycode[year2] + year",
  panel_id = panel_id, vcov_formula = vcov_fm
)
irf_spec1$spec <- "Treated Only + Year FE"

# ==================================================================
# SPECIFICATION 2: Same but with region x year FE
# ==================================================================
cat("=== Spec 2: Separate shock indicators, treated sample, region x year FE ===\n")

irf_spec2 <- lp_two_shocks(
  data = d_treated, outcome = outcome,
  shock1 = "isolated_strict", shock2 = "compound_strict",
  controls = controls_1, horizon = horizon,
  fe = "countrycode[year] + countrycode[year2] + region^year",
  panel_id = panel_id, vcov_formula = vcov_fm
)
irf_spec2$spec <- "Treated Only + Region×Year FE"

# ==================================================================
# SPECIFICATION 3: Full sample with compound/isolated indicators
# ==================================================================
cat("=== Spec 3: Full sample (add never-treated as controls) ===\n")

data_full <- data_full %>%
  group_by(countrycode) %>%
  mutate(
    shock_lag1 = lag(maxwind_95, 1, default = 0L),
    compound_strict = as.integer(maxwind_95 == 1 & shock_lag1 == 1),
    isolated_strict = as.integer(maxwind_95 == 1 & shock_lag1 == 0)
  ) %>%
  ungroup()

irf_spec3 <- lp_two_shocks(
  data = data_full, outcome = outcome,
  shock1 = "isolated_strict", shock2 = "compound_strict",
  controls = controls_1, horizon = horizon,
  fe = "countrycode[year] + countrycode[year2] + region^year",
  panel_id = panel_id, vcov_formula = vcov_fm
)
irf_spec3$spec <- "Full Sample + Region×Year FE"

# ==================================================================
# SPECIFICATION 4: Interaction approach — shock × had_recent_shock
# ==================================================================
cat("=== Spec 4: Interaction approach (shock × recent_shock indicator) ===\n")

d_treated <- d_treated %>%
  mutate(
    had_recent_shock = as.factor(shock_lag1)
  )

irf_spec4 <- lp_panel_inter(
  data = d_treated, outcome = outcome, main_var = "maxwind_95",
  interact_var = "had_recent_shock",
  controls = "l(gdp_diff,1:2) + l(maxwind_95,1:2)",
  horizon = horizon,
  fe = "countrycode[year] + countrycode[year2] + year",
  panel_id = panel_id, vcov_formula = vcov_fm
) %>%
  mutate(
    category = case_when(
      str_detect(category, "^0") ~ "Isolated (no shock in t-1)",
      str_detect(category, "^1") ~ "Compound (shock in t-1)",
      TRUE ~ category
    )
  )
irf_spec4$spec <- "Interaction: Treated + Year FE"

# Same with region x year FE
irf_spec4r <- lp_panel_inter(
  data = d_treated, outcome = outcome, main_var = "maxwind_95",
  interact_var = "had_recent_shock",
  controls = "l(gdp_diff,1:2) + l(maxwind_95,1:2)",
  horizon = horizon,
  fe = "countrycode[year] + countrycode[year2] + region^year",
  panel_id = panel_id, vcov_formula = vcov_fm
) %>%
  mutate(
    category = case_when(
      str_detect(category, "^0") ~ "Isolated (no shock in t-1)",
      str_detect(category, "^1") ~ "Compound (shock in t-1)",
      TRUE ~ category
    )
  )
irf_spec4r$spec <- "Interaction: Treated + Region×Year FE"

# ==================================================================
# SPECIFICATION 5: Broader compound definition (t-1 or t-2)
# ==================================================================
cat("=== Spec 5: Broader compound definition (shock in t-1 or t-2) ===\n")

d_treated <- d_treated %>%
  mutate(
    had_recent_shock_broad = as.factor(as.integer(shock_lag1 == 1 | shock_lag2 == 1))
  )

irf_spec5 <- lp_panel_inter(
  data = d_treated, outcome = outcome, main_var = "maxwind_95",
  interact_var = "had_recent_shock_broad",
  controls = "l(gdp_diff,1:2) + l(maxwind_95,1:2)",
  horizon = horizon,
  fe = "countrycode[year] + countrycode[year2] + year",
  panel_id = panel_id, vcov_formula = vcov_fm
) %>%
  mutate(
    category = case_when(
      str_detect(category, "^0") ~ "Isolated (no shock in t-1 or t-2)",
      str_detect(category, "^1") ~ "Compound (shock in t-1 or t-2)",
      TRUE ~ category
    )
  )
irf_spec5$spec <- "Broad Compound: Treated + Year FE"

# ==================================================================
# SPECIFICATION 6: Baseline IRF (pooled, no compound split)
# ==================================================================
cat("=== Spec 6: Baseline (pooled) for comparison ===\n")

irf_baseline <- lp_panel(
  data = d_treated, outcome = outcome, main_var = "maxwind_95",
  controls = "l(gdp_diff,1:2) + l(maxwind_95,1:2)",
  horizon = horizon,
  fe = "countrycode[year] + countrycode[year2] + year",
  panel_id = panel_id, vcov_formula = vcov_fm
)
irf_baseline$spec <- "Pooled (no compound split)"

# ------------------------------------------------------------------
# Print results at h=5
# ------------------------------------------------------------------
cat("\n========== RESULTS AT h=5 ==========\n")

for (spec_name in c("Treated Only + Year FE",
                     "Treated Only + Region×Year FE",
                     "Full Sample + Region×Year FE")) {
  cat(sprintf("\n--- %s ---\n", spec_name))
  df <- bind_rows(irf_spec1, irf_spec2, irf_spec3) %>%
    filter(spec == spec_name, horizon == 5)
  for (t in c("isolated_strict", "compound_strict", "Difference")) {
    row <- df %>% filter(term == t)
    if (nrow(row) > 0)
      cat(sprintf("  %-20s: %.3f (SE=%.3f) [%.3f, %.3f]\n",
                  t, row$irf_mean, row$se, row$irf_down, row$irf_up))
  }
}

cat("\n--- Interaction approach at h=5 ---\n")
for (spec_name in unique(c(irf_spec4$spec, irf_spec4r$spec))) {
  cat(sprintf("\n  %s:\n", spec_name))
  df <- bind_rows(irf_spec4, irf_spec4r) %>%
    filter(spec == spec_name, horizon == 5)
  for (cat_name in unique(df$category)) {
    row <- df %>% filter(category == cat_name)
    if (nrow(row) > 0)
      cat(sprintf("    %-35s: %.3f (SE=%.3f)\n",
                  cat_name, row$irf_mean, row$se))
  }
}

# ------------------------------------------------------------------
# Save all results
# ------------------------------------------------------------------
saveRDS(list(
  spec1_separate_yearfe    = irf_spec1,
  spec2_separate_regionfe  = irf_spec2,
  spec3_fullsample_regionfe = irf_spec3,
  spec4_interact_yearfe    = irf_spec4,
  spec4r_interact_regionfe = irf_spec4r,
  spec5_broad_yearfe       = irf_spec5,
  baseline                 = irf_baseline
), "estimation_results.rds")
cat("\nSaved: estimation_results.rds\n")
cat("Done.\n")
