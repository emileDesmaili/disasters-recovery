# investigation_composition.R
# Investigate whether sample composition changes pre- vs post-1990
# explain the divergent local projection IRFs.
#
# Key question: Do different (poorer, more vulnerable) countries enter
# the sample post-1990, and does restricting to a balanced panel
# eliminate the pre/post divergence?

library(tidyverse)
library(fixest)
library(haven)
source("../emileRegs.R")

cat("============================================================\n")
cat("  SAMPLE COMPOSITION INVESTIGATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------------
# 1. Load and construct data identically to lp_gdp_shocks.R
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

data <- data %>%
  mutate(
    maxwind_p95 = quantile(maxwind[year >= 1970 & year <= 2014], 0.95, na.rm = TRUE),
    energy_p95  = quantile(energy[year >= 1970 & year <= 2014],  0.95, na.rm = TRUE),
    nlands_p95  = quantile(nlands[year >= 1970 & year <= 2014],  0.95, na.rm = TRUE)
  ) %>%
  mutate(
    maxwind_95 = as.integer(maxwind_p95 > 0 & maxwind >= maxwind_p95),
    energy_95  = as.integer(energy_p95  > 0 & energy  >= energy_p95),
    nlands_95  = as.integer(nlands_p95  > 0 & nlands  >= nlands_p95)
  )

data <- data %>%
  mutate(
    maxwind_sd = maxwind / sd(maxwind[maxwind > 0], na.rm = TRUE),
    energy_sd  = energy  / sd(energy[energy > 0],   na.rm = TRUE),
    nlands_sd  = nlands  / sd(nlands[nlands > 0],   na.rm = TRUE)
  )

# Compute real GDP per capita in levels for comparison
# (loggdp is ln_real_gdp_usd_pc, so exponentiate)
data <- data %>%
  mutate(gdp_pc = exp(loggdp))

# ------------------------------------------------------------------
# 2. Count unique countries pre- vs post-1990
# ------------------------------------------------------------------
cat("============================================================\n")
cat("  PART 1: Country Counts by Period\n")
cat("============================================================\n\n")

# Use non-missing loggdp observations (the estimation sample)
est_sample <- data %>% filter(!is.na(loggdp))

countries_pre  <- est_sample %>% filter(period == "Pre-1990") %>%
  pull(countrycode) %>% unique()
countries_post <- est_sample %>% filter(period == "Post-1990") %>%
  pull(countrycode) %>% unique()

cat(sprintf("Countries with GDP data pre-1990 (1969-1990):  %d\n", length(countries_pre)))
cat(sprintf("Countries with GDP data post-1990 (1991+):     %d\n", length(countries_post)))
cat(sprintf("Countries in BOTH periods:                     %d\n",
            length(intersect(countries_pre, countries_post))))
cat(sprintf("Countries ONLY pre-1990:                       %d\n",
            length(setdiff(countries_pre, countries_post))))
cat(sprintf("Countries ONLY post-1990:                      %d\n",
            length(setdiff(countries_post, countries_pre))))

# ------------------------------------------------------------------
# 3. Identify countries unique to each period
# ------------------------------------------------------------------
cat("\n============================================================\n")
cat("  PART 2: Countries Appearing ONLY in One Period\n")
cat("============================================================\n\n")

only_pre  <- sort(setdiff(countries_pre, countries_post))
only_post <- sort(setdiff(countries_post, countries_pre))
both      <- sort(intersect(countries_pre, countries_post))

cat("Countries ONLY in pre-1990:\n")
if (length(only_pre) > 0) {
  cat(paste("  ", only_pre, collapse = "\n"), "\n")
} else {
  cat("  (none)\n")
}

cat("\nCountries ONLY in post-1990:\n")
if (length(only_post) > 0) {
  cat(paste("  ", only_post, collapse = "\n"), "\n")
} else {
  cat("  (none)\n")
}

# ------------------------------------------------------------------
# 4. Compare average GDP per capita across periods
# ------------------------------------------------------------------
cat("\n============================================================\n")
cat("  PART 3: Average GDP per Capita by Period\n")
cat("============================================================\n\n")

gdp_by_period <- est_sample %>%
  group_by(period) %>%
  summarise(
    n_obs     = n(),
    n_countries = n_distinct(countrycode),
    mean_gdp_pc = mean(gdp_pc, na.rm = TRUE),
    median_gdp_pc = median(gdp_pc, na.rm = TRUE),
    sd_gdp_pc = sd(gdp_pc, na.rm = TRUE),
    .groups = "drop"
  )

print(gdp_by_period)

# Also compare using only the first observation per country in each period
# to avoid weighting by panel length
gdp_first_obs <- est_sample %>%
  group_by(countrycode, period) %>%
  summarise(mean_gdp_pc = mean(gdp_pc, na.rm = TRUE), .groups = "drop")

cat("\nCross-country average GDP p.c. (averaging within country first):\n")
gdp_first_obs %>%
  group_by(period) %>%
  summarise(
    n_countries = n(),
    mean_gdp_pc = mean(mean_gdp_pc, na.rm = TRUE),
    median_gdp_pc = median(mean_gdp_pc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# ------------------------------------------------------------------
# 5. Are "entering" countries (post-1990 only) systematically
#    poorer or more cyclone-exposed?
# ------------------------------------------------------------------
cat("\n============================================================\n")
cat("  PART 4: Characteristics of Entering vs Continuing Countries\n")
cat("============================================================\n\n")

# Tag each country
est_sample <- est_sample %>%
  mutate(
    country_group = case_when(
      countrycode %in% only_pre  ~ "Only Pre-1990",
      countrycode %in% only_post ~ "Only Post-1990",
      TRUE                       ~ "Both Periods"
    )
  )

# Compare characteristics in the POST-1990 period specifically
post_comparison <- est_sample %>%
  filter(period == "Post-1990") %>%
  group_by(country_group) %>%
  summarise(
    n_countries = n_distinct(countrycode),
    mean_gdp_pc = mean(gdp_pc, na.rm = TRUE),
    median_gdp_pc = median(gdp_pc, na.rm = TRUE),
    mean_maxwind = mean(maxwind, na.rm = TRUE),
    frac_shock_95 = mean(maxwind_95, na.rm = TRUE),
    mean_nlands = mean(nlands, na.rm = TRUE),
    .groups = "drop"
  )

cat("Comparison of country groups in the POST-1990 period:\n")
print(post_comparison)

# Country-level averages (to avoid panel-length weighting)
country_level_post <- est_sample %>%
  filter(period == "Post-1990") %>%
  group_by(countrycode, country_group) %>%
  summarise(
    mean_gdp_pc = mean(gdp_pc, na.rm = TRUE),
    mean_maxwind = mean(maxwind, na.rm = TRUE),
    frac_shock_95 = mean(maxwind_95, na.rm = TRUE),
    total_shocks_95 = sum(maxwind_95, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nCountry-level averages (one obs per country), post-1990:\n")
country_level_post %>%
  group_by(country_group) %>%
  summarise(
    n_countries = n(),
    mean_gdp_pc = mean(mean_gdp_pc, na.rm = TRUE),
    median_gdp_pc = median(mean_gdp_pc, na.rm = TRUE),
    mean_maxwind = mean(mean_maxwind, na.rm = TRUE),
    mean_frac_shock = mean(frac_shock_95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# t-test: are entering countries poorer?
if (length(only_post) >= 2) {
  entering  <- country_level_post %>% filter(country_group == "Only Post-1990") %>% pull(mean_gdp_pc)
  continuing <- country_level_post %>% filter(country_group == "Both Periods") %>% pull(mean_gdp_pc)

  cat("\nt-test: GDP p.c. of entering vs continuing countries (post-1990):\n")
  tt <- t.test(entering, continuing)
  cat(sprintf("  Entering mean:    $%.0f\n", mean(entering, na.rm = TRUE)))
  cat(sprintf("  Continuing mean:  $%.0f\n", mean(continuing, na.rm = TRUE)))
  cat(sprintf("  Difference:       $%.0f\n", mean(entering, na.rm = TRUE) - mean(continuing, na.rm = TRUE)))
  cat(sprintf("  t-stat:           %.3f\n", tt$statistic))
  cat(sprintf("  p-value:          %.4f\n", tt$p.value))
}

# Also check cyclone exposure
if (length(only_post) >= 2) {
  entering_wind  <- country_level_post %>% filter(country_group == "Only Post-1990") %>% pull(mean_maxwind)
  continuing_wind <- country_level_post %>% filter(country_group == "Both Periods") %>% pull(mean_maxwind)

  cat("\nt-test: Mean max wind of entering vs continuing countries (post-1990):\n")
  tt2 <- t.test(entering_wind, continuing_wind)
  cat(sprintf("  Entering mean wind:    %.3f m/s\n", mean(entering_wind, na.rm = TRUE)))
  cat(sprintf("  Continuing mean wind:  %.3f m/s\n", mean(continuing_wind, na.rm = TRUE)))
  cat(sprintf("  t-stat:                %.3f\n", tt2$statistic))
  cat(sprintf("  p-value:               %.4f\n", tt2$p.value))
}

# ------------------------------------------------------------------
# 6. Which post-1990-only countries actually experience shocks?
# ------------------------------------------------------------------
cat("\n============================================================\n")
cat("  PART 5: Shock Events Among Entering Countries\n")
cat("============================================================\n\n")

entering_shocks <- est_sample %>%
  filter(country_group == "Only Post-1990", maxwind_95 == 1) %>%
  select(countrycode, year, maxwind, loggdp, gdp_pc) %>%
  arrange(countrycode, year)

cat(sprintf("Number of shock events (maxwind_95=1) among entering countries: %d\n",
            nrow(entering_shocks)))
cat(sprintf("Number of entering countries experiencing at least one shock: %d\n",
            n_distinct(entering_shocks$countrycode)))

if (nrow(entering_shocks) > 0) {
  cat("\nShock events among entering countries:\n")
  print(as.data.frame(entering_shocks), row.names = FALSE)
}

# ------------------------------------------------------------------
# 7. Balanced panel LP: restrict to countries in BOTH periods
# ------------------------------------------------------------------
cat("\n============================================================\n")
cat("  PART 6: Balanced Panel LP (Countries in Both Periods Only)\n")
cat("============================================================\n\n")

data_balanced <- data %>% filter(countrycode %in% both)

cat(sprintf("Full sample countries:     %d\n", n_distinct(data$countrycode)))
cat(sprintf("Balanced sample countries: %d\n", n_distinct(data_balanced$countrycode)))
cat(sprintf("Full sample obs:           %d\n", nrow(data)))
cat(sprintf("Balanced sample obs:       %d\n", nrow(data_balanced)))

# LP settings (same as original)
outcome  <- "loggdp"
horizon  <- 10
make_controls <- function(shock_var) {
  paste0("l(gdp_diff,1:2) + l(", shock_var, ",1:2)")
}
fe       <- "countrycode[year] + countrycode[year2] + year"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

# --- A. Full sample period-split (reproduce original) ---
cat("\nEstimating period-split IRF on FULL sample...\n")
irfs_full <- lp_panel_inter(
  data = data, outcome = outcome, main_var = "maxwind_95",
  interact_var = "period", controls = make_controls("maxwind_95"),
  horizon = horizon, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
)

cat("\nFull sample IRFs by period (maxwind_95):\n")
irfs_full %>%
  select(horizon, category, irf_mean, se, irf_down, irf_up) %>%
  pivot_wider(
    names_from = category,
    values_from = c(irf_mean, se, irf_down, irf_up),
    names_glue = "{category}_{.value}"
  ) %>%
  print(n = 20)

# --- B. Balanced panel period-split ---
cat("\nEstimating period-split IRF on BALANCED sample...\n")
irfs_balanced <- lp_panel_inter(
  data = data_balanced, outcome = outcome, main_var = "maxwind_95",
  interact_var = "period", controls = make_controls("maxwind_95"),
  horizon = horizon, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
)

cat("\nBalanced sample IRFs by period (maxwind_95):\n")
irfs_balanced %>%
  select(horizon, category, irf_mean, se, irf_down, irf_up) %>%
  pivot_wider(
    names_from = category,
    values_from = c(irf_mean, se, irf_down, irf_up),
    names_glue = "{category}_{.value}"
  ) %>%
  print(n = 20)

# --- C. Compare peak (horizon 5-10) effects ---
cat("\n============================================================\n")
cat("  PART 7: Summary Comparison of Peak Effects (h=5-10)\n")
cat("============================================================\n\n")

summarise_peak <- function(df, label) {
  df %>%
    filter(horizon >= 5 & horizon <= 10) %>%
    group_by(category) %>%
    summarise(
      mean_irf = mean(irf_mean),
      min_irf  = min(irf_mean),
      .groups = "drop"
    ) %>%
    mutate(sample = label)
}

peak_comparison <- bind_rows(
  summarise_peak(irfs_full, "Full Sample"),
  summarise_peak(irfs_balanced, "Balanced Panel")
)

print(peak_comparison)

# --- D. Plot comparison ---
cat("\nGenerating comparison plot...\n")

irfs_full$sample <- "Full Sample"
irfs_balanced$sample <- "Balanced Panel"
irfs_compare <- bind_rows(irfs_full, irfs_balanced)

# Clean up category names
irfs_compare <- irfs_compare %>%
  mutate(
    category = str_replace_all(category, "_", "-"),
    category = case_when(
      str_detect(category, "Pre")  ~ "Pre-1990",
      str_detect(category, "Post") ~ "Post-1990",
      TRUE ~ category
    ),
    category = factor(category, levels = c("Pre-1990", "Post-1990"))
  )

p <- ggplot(irfs_compare, aes(x = horizon, y = irf_mean,
                               color = category, fill = category)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~sample) +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_fill_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "Horizon (years)",
    y = "GDP response to maxwind_95 shock",
    color = "Period", fill = "Period",
    title = "Period-Split IRFs: Full Sample vs Balanced Panel",
    subtitle = "Does restricting to countries present in both periods eliminate the divergence?"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave("figures/investigation_composition.png", p, width = 12, height = 5, dpi = 300)
cat("Plot saved to figures/investigation_composition.png\n")

# ------------------------------------------------------------------
# 8. Final verdict
# ------------------------------------------------------------------
cat("\n============================================================\n")
cat("  VERDICT\n")
cat("============================================================\n\n")

# Compare the pre-post gap in full vs balanced
gap_full <- irfs_full %>%
  filter(horizon >= 5 & horizon <= 10) %>%
  group_by(category) %>%
  summarise(mean_irf = mean(irf_mean), .groups = "drop") %>%
  pivot_wider(names_from = category, values_from = mean_irf)

gap_balanced <- irfs_balanced %>%
  filter(horizon >= 5 & horizon <= 10) %>%
  group_by(category) %>%
  summarise(mean_irf = mean(irf_mean), .groups = "drop") %>%
  pivot_wider(names_from = category, values_from = mean_irf)

# The column names come from the factor levels after cleaning
pre_col_full  <- names(gap_full)[str_detect(names(gap_full), "Pre")]
post_col_full <- names(gap_full)[str_detect(names(gap_full), "Post")]
pre_col_bal   <- names(gap_balanced)[str_detect(names(gap_balanced), "Pre")]
post_col_bal  <- names(gap_balanced)[str_detect(names(gap_balanced), "Post")]

if (length(pre_col_full) > 0 & length(post_col_full) > 0) {
  gap_f <- gap_full[[post_col_full]] - gap_full[[pre_col_full]]
  cat(sprintf("Full sample:     Post-Pre gap at h=5-10: %.4f (%.2f pp)\n", gap_f, gap_f * 100))
}

if (length(pre_col_bal) > 0 & length(post_col_bal) > 0) {
  gap_b <- gap_balanced[[post_col_bal]] - gap_balanced[[pre_col_bal]]
  cat(sprintf("Balanced panel:  Post-Pre gap at h=5-10: %.4f (%.2f pp)\n", gap_b, gap_b * 100))
}

if (exists("gap_f") & exists("gap_b")) {
  pct_explained <- (1 - gap_b / gap_f) * 100
  cat(sprintf("\nFraction of gap explained by composition: %.1f%%\n", pct_explained))

  if (abs(pct_explained) > 50) {
    cat("\n=> CONCLUSION: Sample composition appears to be a MAJOR driver of the\n")
    cat("   divergent IRFs. The entering countries substantially change the result.\n")
  } else if (abs(pct_explained) > 20) {
    cat("\n=> CONCLUSION: Sample composition explains a MODERATE share of the divergence.\n")
    cat("   It contributes but is not the sole explanation.\n")
  } else {
    cat("\n=> CONCLUSION: Sample composition does NOT appear to drive the divergence.\n")
    cat("   The pre/post gap persists even on the balanced panel.\n")
    cat("   The puzzle is about within-country changes over time, not composition.\n")
  }
}

cat("\nDone.\n")
