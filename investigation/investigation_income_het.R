# investigation_income_het.R
# Finer-grained income heterogeneity in cyclone-GDP IRFs
# Bins countries into quartiles of GDP per capita, runs period-interacted LPs
# Also identifies which countries show positive pre-1990 responses

library(tidyverse)
library(fixest)
library(haven)
source("../emileRegs.R")

# ------------------------------------------------------------------
# 1. Load and prepare data (same pipeline as lp_gdp_shocks.R)
# ------------------------------------------------------------------
pwt <- read_stata("../raw_data/pwt_clean.dta")
tcs <- read_stata("../raw_data/ibtracs_clean.dta")
data <- pwt %>% left_join(tcs, by = c("year", "countrycode"))

knot_to_ms <- 0.514444

data <- data %>%
  mutate(
    maxwind = replace_na(max_ann_wind_i_nob, 0) * knot_to_ms,
    energy  = replace_na(sum_ann_energy_i_nob, 0),
    nlands  = replace_na(sum_a_lands_nob, 0),
    loggdp   = ln_real_gdp_usd_pc,
    gdp_diff = growth_real_gdp_usd_pc,
    year2    = year^2,
    period   = ifelse(year <= 1990, "Pre-1990", "Post-1990")
  ) %>%
  filter(year >= 1969)

# Global 95th-percentile thresholds
data <- data %>%
  mutate(
    maxwind_p95 = quantile(maxwind[year >= 1970 & year <= 2014], 0.95, na.rm = TRUE)
  ) %>%
  mutate(
    maxwind_95 = as.integer(maxwind_p95 > 0 & maxwind >= maxwind_p95)
  )

# LP settings
outcome  <- "loggdp"
horizon  <- 10
make_controls <- function(shock_var) {
  paste0("l(gdp_diff,1:2) + l(", shock_var, ",1:2)")
}
fe       <- "countrycode[year] + countrycode[year2] + year"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

irf_theme <- theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# ------------------------------------------------------------------
# 2. Create income quartiles based on country mean GDP per capita
# ------------------------------------------------------------------
country_gdp <- data %>%
  filter(!is.na(loggdp)) %>%
  group_by(countrycode) %>%
  summarise(
    mean_loggdp = mean(loggdp, na.rm = TRUE),
    mean_gdp_pc = exp(mean(loggdp, na.rm = TRUE)),
    n_shocks_95 = sum(maxwind_95, na.rm = TRUE),
    total_wind  = sum(maxwind, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    income_quartile = ntile(mean_loggdp, 4),
    income_quartile_label = factor(
      income_quartile,
      levels = 1:4,
      labels = c("Q1 (Poorest)", "Q2", "Q3", "Q4 (Richest)")
    )
  )

data <- data %>%
  left_join(
    country_gdp %>% select(countrycode, income_quartile, income_quartile_label,
                            mean_loggdp, mean_gdp_pc),
    by = "countrycode"
  )

# Print quartile composition
cat("\n============================================================\n")
cat("Income Quartile Composition\n")
cat("============================================================\n\n")

quartile_summary <- country_gdp %>%
  group_by(income_quartile_label) %>%
  summarise(
    n_countries = n(),
    min_gdp_pc = min(mean_gdp_pc),
    max_gdp_pc = max(mean_gdp_pc),
    n_with_shocks = sum(n_shocks_95 > 0),
    total_shocks = sum(n_shocks_95),
    .groups = "drop"
  )
print(as.data.frame(quartile_summary), row.names = FALSE)

# ------------------------------------------------------------------
# 3. Run period-interacted LP within each quartile
# ------------------------------------------------------------------
cat("\n============================================================\n")
cat("Period-Interacted LP by Income Quartile\n")
cat("============================================================\n\n")

irfs_quartile <- map_dfr(1:4, function(q) {
  label <- levels(country_gdp$income_quartile_label)[q]
  cat(sprintf("  Estimating LP for %s...\n", label))
  sub_data <- data %>% filter(income_quartile == q)

  tryCatch({
    lp_panel_inter(
      data = sub_data,
      outcome = outcome,
      main_var = "maxwind_95",
      interact_var = "period",
      controls = make_controls("maxwind_95"),
      horizon = horizon,
      fe = fe,
      panel_id = panel_id,
      vcov_formula = vcov_fm
    ) %>%
      mutate(income_group = label)
  }, error = function(e) {
    cat(sprintf("    Error for %s: %s\n", label, e$message))
    return(NULL)
  })
})

# Clean category names
irfs_quartile <- irfs_quartile %>%
  mutate(category = ifelse(grepl("Pre", category), "Pre-1990", "Post-1990"))

# Print summary at h=5
cat("\n--- IRF at h=5 by income quartile and period ---\n")
irfs_quartile %>%
  filter(horizon == 5) %>%
  select(income_group, category, irf_mean, se, irf_down, irf_up) %>%
  arrange(income_group, category) %>%
  as.data.frame() %>%
  print(digits = 4, row.names = FALSE)

# ------------------------------------------------------------------
# 4. Plot: 4-panel figure, one per quartile
# ------------------------------------------------------------------
fig_quartile <- irfs_quartile %>%
  mutate(
    category = factor(category, levels = c("Pre-1990", "Post-1990")),
    income_group = factor(income_group,
                          levels = c("Q1 (Poorest)", "Q2", "Q3", "Q4 (Richest)"))
  ) %>%
  ggplot(aes(x = horizon, y = irf_mean, color = category, fill = category)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~income_group, ncol = 4) +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_fill_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "Horizon (years)", y = "GDP response",
    color = "Period", fill = "Period"
  ) +
  irf_theme +
  theme(legend.position = "top")

ggsave("figures/irf_income_quartiles_by_period.png", fig_quartile,
       width = 14, height = 4.5, dpi = 300)

# ------------------------------------------------------------------
# 5. Identify countries with positive pre-1990 responses
# ------------------------------------------------------------------
cat("\n\n============================================================\n")
cat("Countries with Positive Pre-1990 IRFs (by country)\n")
cat("============================================================\n\n")

# Run country-level LPs for cyclone-exposed countries in top 2 quartiles
exposed_rich <- country_gdp %>%
  filter(income_quartile >= 3, n_shocks_95 >= 2) %>%
  arrange(desc(n_shocks_95))

cat(sprintf("Running country-level LPs for %d high-income exposed countries...\n\n",
            nrow(exposed_rich)))

country_irfs <- map_dfr(exposed_rich$countrycode, function(cc) {
  sub_data <- data %>% filter(countrycode == cc)

  # Run separate pre and post LPs (not interacted - too few obs per country)
  map_dfr(c("Pre-1990", "Post-1990"), function(per) {
    per_data <- sub_data %>% filter(period == per)

    # Need enough shock observations to estimate
    if (sum(per_data$maxwind_95, na.rm = TRUE) < 1) return(NULL)

    tryCatch({
      lp_panel(
        data = per_data,
        outcome = outcome,
        main_var = "maxwind_95",
        controls = "l(gdp_diff,1:2)",
        horizon = min(horizon, nrow(per_data) - 5),
        fe = "year",
        panel_id = panel_id,
        vcov_formula = NW ~ year
      ) %>%
        mutate(countrycode = cc, period = per,
               mean_gdp_pc = exposed_rich$mean_gdp_pc[exposed_rich$countrycode == cc],
               n_shocks = exposed_rich$n_shocks_95[exposed_rich$countrycode == cc])
    }, error = function(e) {
      return(NULL)
    })
  })
})

if (nrow(country_irfs) > 0) {
  # Summarize at h=5
  country_summary <- country_irfs %>%
    filter(horizon == 5) %>%
    select(countrycode, period, irf_mean, se, mean_gdp_pc, n_shocks) %>%
    pivot_wider(names_from = period, values_from = c(irf_mean, se),
                names_glue = "{period}_{.value}") %>%
    arrange(desc(`Pre-1990_irf_mean`))

  cat("--- Country-level IRF at h=5 (Q3+Q4 countries with >=2 shocks) ---\n")
  cat("    Sorted by Pre-1990 response (most positive first)\n\n")
  print(as.data.frame(country_summary), digits = 3, row.names = FALSE)

  # Flag countries with positive pre-1990 response
  positive_pre <- country_summary %>%
    filter(!is.na(`Pre-1990_irf_mean`) & `Pre-1990_irf_mean` > 0)

  if (nrow(positive_pre) > 0) {
    cat(sprintf("\n\nCountries with POSITIVE Pre-1990 h=5 response: %s\n",
                paste(positive_pre$countrycode, collapse = ", ")))
    cat(sprintf("  Mean GDP per capita range: $%.0f - $%.0f\n",
                min(positive_pre$mean_gdp_pc), max(positive_pre$mean_gdp_pc)))
  }
}

# ------------------------------------------------------------------
# 6. Also run with quintiles for even finer resolution
# ------------------------------------------------------------------
cat("\n\n============================================================\n")
cat("Quintile Analysis (5 bins)\n")
cat("============================================================\n\n")

country_gdp <- country_gdp %>%
  mutate(
    income_quintile = ntile(mean_loggdp, 5),
    income_quintile_label = factor(
      income_quintile,
      levels = 1:5,
      labels = c("Q1 (Poorest)", "Q2", "Q3", "Q4", "Q5 (Richest)")
    )
  )

data <- data %>%
  left_join(
    country_gdp %>% select(countrycode, income_quintile, income_quintile_label),
    by = "countrycode"
  )

irfs_quintile <- map_dfr(1:5, function(q) {
  label <- levels(country_gdp$income_quintile_label)[q]
  cat(sprintf("  Estimating LP for %s...\n", label))
  sub_data <- data %>% filter(income_quintile == q)

  tryCatch({
    lp_panel_inter(
      data = sub_data,
      outcome = outcome,
      main_var = "maxwind_95",
      interact_var = "period",
      controls = make_controls("maxwind_95"),
      horizon = horizon,
      fe = fe,
      panel_id = panel_id,
      vcov_formula = vcov_fm
    ) %>%
      mutate(income_group = label)
  }, error = function(e) {
    cat(sprintf("    Error for %s: %s\n", label, e$message))
    return(NULL)
  })
})

irfs_quintile <- irfs_quintile %>%
  mutate(category = ifelse(grepl("Pre", category), "Pre-1990", "Post-1990"))

fig_quintile <- irfs_quintile %>%
  mutate(
    category = factor(category, levels = c("Pre-1990", "Post-1990")),
    income_group = factor(income_group,
                          levels = c("Q1 (Poorest)", "Q2", "Q3", "Q4", "Q5 (Richest)"))
  ) %>%
  ggplot(aes(x = horizon, y = irf_mean, color = category, fill = category)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~income_group, ncol = 5) +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_fill_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "Horizon (years)", y = "GDP response",
    color = "Period", fill = "Period"
  ) +
  irf_theme +
  theme(legend.position = "top")

ggsave("figures/irf_income_quintiles_by_period.png", fig_quintile,
       width = 16, height = 4.5, dpi = 300)

# Print summary table
cat("\n--- IRF at h=5 by income quintile and period ---\n")
irfs_quintile %>%
  filter(horizon == 5) %>%
  select(income_group, category, irf_mean, se) %>%
  mutate(irf_pct = paste0(round(100 * irf_mean, 2), "%")) %>%
  arrange(income_group, category) %>%
  as.data.frame() %>%
  print(digits = 4, row.names = FALSE)

# ------------------------------------------------------------------
# 7. Exposure composition: which countries are hit in each period?
# ------------------------------------------------------------------
cat("\n\n============================================================\n")
cat("Exposure Composition by Period\n")
cat("============================================================\n\n")

# Which countries experience p95 shocks in each period?
exposure_by_period <- data %>%
  filter(maxwind_95 == 1) %>%
  group_by(countrycode, period) %>%
  summarise(
    n_shocks = n(),
    mean_wind = mean(maxwind, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    country_gdp %>% select(countrycode, mean_gdp_pc, income_quartile_label),
    by = "countrycode"
  )

# Pivot to show pre vs post side by side
exposure_wide <- exposure_by_period %>%
  pivot_wider(
    names_from = period,
    values_from = c(n_shocks, mean_wind),
    names_glue = "{period}_{.value}"
  ) %>%
  mutate(
    across(starts_with("Pre-1990"), ~replace_na(.x, 0)),
    across(starts_with("Post-1990"), ~replace_na(.x, 0)),
    exposed_both   = `Pre-1990_n_shocks` > 0 & `Post-1990_n_shocks` > 0,
    exposed_pre_only  = `Pre-1990_n_shocks` > 0 & `Post-1990_n_shocks` == 0,
    exposed_post_only = `Pre-1990_n_shocks` == 0 & `Post-1990_n_shocks` > 0
  ) %>%
  arrange(desc(`Pre-1990_n_shocks` + `Post-1990_n_shocks`))

cat(sprintf("Countries with p95 shocks in BOTH periods: %d\n",
            sum(exposure_wide$exposed_both)))
cat(sprintf("Countries with p95 shocks PRE-1990 ONLY:   %d\n",
            sum(exposure_wide$exposed_pre_only)))
cat(sprintf("Countries with p95 shocks POST-1990 ONLY:  %d\n",
            sum(exposure_wide$exposed_post_only)))

cat("\n--- Countries exposed in both periods ---\n")
exposure_wide %>%
  filter(exposed_both) %>%
  select(countrycode, income_quartile_label, mean_gdp_pc,
         `Pre-1990_n_shocks`, `Post-1990_n_shocks`) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n--- Countries exposed ONLY pre-1990 ---\n")
exposure_wide %>%
  filter(exposed_pre_only) %>%
  select(countrycode, income_quartile_label, mean_gdp_pc,
         `Pre-1990_n_shocks`) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n--- Countries exposed ONLY post-1990 ---\n")
exposure_wide %>%
  filter(exposed_post_only) %>%
  select(countrycode, income_quartile_label, mean_gdp_pc,
         `Post-1990_n_shocks`) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# Composition by income quartile
cat("\n--- Exposure composition by income quartile and period ---\n")
exposure_composition <- data %>%
  filter(maxwind_95 == 1) %>%
  group_by(period, income_quartile_label) %>%
  summarise(
    n_shock_events = n(),
    n_countries_hit = n_distinct(countrycode),
    mean_wind = mean(maxwind, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(period, income_quartile_label)
print(as.data.frame(exposure_composition), digits = 3, row.names = FALSE)

# Figure: shock counts by quartile and period
fig_composition <- exposure_composition %>%
  mutate(
    period = factor(period, levels = c("Pre-1990", "Post-1990")),
    income_quartile_label = factor(income_quartile_label,
                                   levels = c("Q1 (Poorest)", "Q2", "Q3", "Q4 (Richest)"))
  ) %>%
  ggplot(aes(x = income_quartile_label, y = n_shock_events, fill = period)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(aes(label = n_countries_hit, y = n_shock_events + 1),
            position = position_dodge(width = 0.9), size = 3.5, vjust = 0) +
  scale_fill_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  labs(
    x = "Income Quartile", y = "Number of p95 shock events",
    fill = "Period",
    caption = "Numbers above bars = distinct countries hit"
  ) +
  irf_theme +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 11))

ggsave("figures/exposure_composition_by_income.png", fig_composition,
       width = 9, height = 5, dpi = 300)

# ------------------------------------------------------------------
# 8. Robustness: restrict to countries exposed in BOTH periods
# ------------------------------------------------------------------
cat("\n\n============================================================\n")
cat("Robustness: Restrict to countries exposed in both periods\n")
cat("============================================================\n\n")

both_exposed <- exposure_wide %>%
  filter(exposed_both) %>%
  pull(countrycode)

cat(sprintf("Running LP on %d countries with shocks in both periods...\n\n",
            length(both_exposed)))

data_both <- data %>% filter(countrycode %in% both_exposed)

irfs_both_only <- lp_panel_inter(
  data = data_both,
  outcome = outcome,
  main_var = "maxwind_95",
  interact_var = "period",
  controls = make_controls("maxwind_95"),
  horizon = horizon,
  fe = fe,
  panel_id = panel_id,
  vcov_formula = vcov_fm
) %>%
  mutate(
    category = ifelse(grepl("Pre", category), "Pre-1990", "Post-1990"),
    sample = "Both-period exposed"
  )

# Compare to full sample
irfs_full <- lp_panel_inter(
  data = data,
  outcome = outcome,
  main_var = "maxwind_95",
  interact_var = "period",
  controls = make_controls("maxwind_95"),
  horizon = horizon,
  fe = fe,
  panel_id = panel_id,
  vcov_formula = vcov_fm
) %>%
  mutate(
    category = ifelse(grepl("Pre", category), "Pre-1990", "Post-1990"),
    sample = "Full sample"
  )

irfs_comparison <- bind_rows(irfs_full, irfs_both_only)

cat("--- Comparison at h=5 ---\n")
irfs_comparison %>%
  filter(horizon == 5) %>%
  select(sample, category, irf_mean, se, irf_down, irf_up) %>%
  arrange(sample, category) %>%
  as.data.frame() %>%
  print(digits = 4, row.names = FALSE)

fig_comparison <- irfs_comparison %>%
  mutate(
    category = factor(category, levels = c("Pre-1990", "Post-1990")),
    sample = factor(sample, levels = c("Full sample", "Both-period exposed"))
  ) %>%
  ggplot(aes(x = horizon, y = irf_mean, color = category, fill = category,
             linetype = sample)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.10, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_fill_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_linetype_manual(values = c("Full sample" = "solid",
                                    "Both-period exposed" = "dashed")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "Horizon (years)", y = "GDP response",
    color = "Period", fill = "Period", linetype = "Sample"
  ) +
  irf_theme +
  theme(legend.position = "top")

ggsave("figures/irf_both_period_exposed_comparison.png", fig_comparison,
       width = 9, height = 5, dpi = 300)

cat("\n\nFigures saved:\n")
cat("  - figures/irf_income_quartiles_by_period.png\n")
cat("  - figures/irf_income_quintiles_by_period.png\n")
cat("  - figures/exposure_composition_by_income.png\n")
cat("  - figures/irf_both_period_exposed_comparison.png\n")
cat("\nDone.\n")
