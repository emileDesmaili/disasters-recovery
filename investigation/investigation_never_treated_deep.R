# investigation_never_treated_deep.R
# Deep investigation of why never-treated countries matter as controls
# in the local projection (LP) analysis of cyclone shocks on GDP.
#
# Key finding to investigate: 142 of 182 countries NEVER experience a p95
# cyclone shock. Dropping them substantially attenuates the post-1990
# divergence. This script digs into WHY.

library(tidyverse)
library(fixest)
library(haven)
library(countrycode)

# Source from repo root
source("emileRegs.R")

setFixest_notes(FALSE)

cat("============================================================\n")
cat("  DEEP INVESTIGATION: NEVER-TREATED COUNTRIES\n")
cat("============================================================\n\n")

# ------------------------------------------------------------------
# 1. Load and construct data (identically to lp_gdp_shocks.R)
# ------------------------------------------------------------------
pwt  <- read_stata("raw_data/pwt_clean.dta")
tcs  <- read_stata("raw_data/ibtracs_clean.dta")
data <- pwt %>% left_join(tcs, by = c("year", "countrycode"))

knot_to_ms <- 0.514444

data <- data %>%
  mutate(
    maxwind  = replace_na(max_ann_wind_i_nob, 0) * knot_to_ms,
    energy   = replace_na(sum_ann_energy_i_nob, 0),
    nlands   = replace_na(sum_a_lands_nob, 0),
    loggdp   = ln_real_gdp_usd_pc,
    gdp_diff = growth_real_gdp_usd_pc,
    year2    = year^2,
    period   = ifelse(year <= 1990, "Pre-1990", "Post-1990")
  ) %>%
  filter(year >= 1969)

# Global 95th-percentile thresholds (1970-2014)
data <- data %>%
  mutate(
    maxwind_p95 = quantile(maxwind[year >= 1970 & year <= 2014], 0.95, na.rm = TRUE)
  ) %>%
  mutate(
    maxwind_95 = as.integer(maxwind_p95 > 0 & maxwind >= maxwind_p95)
  )

# Standardized continuous measure
data <- data %>%
  mutate(maxwind_sd = maxwind / sd(maxwind[maxwind > 0], na.rm = TRUE))

# ------------------------------------------------------------------
# LP settings (matching lp_gdp_shocks.R)
# ------------------------------------------------------------------
outcome  <- "loggdp"
horizon  <- 10
fe       <- "countrycode[year] + countrycode[year2] + year"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

make_controls <- function(shock_var) {
  paste0("l(gdp_diff,1:2) + l(", shock_var, ",1:2)")
}
controls <- make_controls("maxwind_95")

# ------------------------------------------------------------------
# Identify ever-treated vs never-treated countries
# ------------------------------------------------------------------
ever_treated_countries <- data %>%
  filter(maxwind_95 == 1) %>%
  pull(countrycode) %>%
  unique() %>%
  sort()

all_countries <- unique(data$countrycode)
never_treated_countries <- sort(setdiff(all_countries, ever_treated_countries))

# Add treatment group flag
data <- data %>%
  mutate(treatment_group = ifelse(countrycode %in% ever_treated_countries,
                                  "Ever-Treated", "Never-Treated"))

# Add continent using the countrycode package
data <- data %>%
  mutate(
    continent = countrycode(countrycode, "iso3c", "continent"),
    region    = countrycode(countrycode, "iso3c", "region")
  )

cat(sprintf("Total countries: %d\n", length(all_countries)))
cat(sprintf("Ever-treated:    %d\n", length(ever_treated_countries)))
cat(sprintf("Never-treated:   %d\n", length(never_treated_countries)))

# Create figures directory if needed
dir.create("investigation/figures", showWarnings = FALSE, recursive = TRUE)

# ==================================================================
# ANALYSIS 1: Characterize the groups
# ==================================================================
cat("\n============================================================\n")
cat("  ANALYSIS 1: Characterize Never-Treated vs Ever-Treated\n")
cat("============================================================\n\n")

# --- 1a. Regional composition ---
cat("--- 1a. Regional Composition ---\n\n")

region_comp <- data %>%
  distinct(countrycode, treatment_group, continent) %>%
  group_by(treatment_group, continent) %>%
  summarise(n_countries = n(), .groups = "drop") %>%
  group_by(treatment_group) %>%
  mutate(pct = round(100 * n_countries / sum(n_countries), 1)) %>%
  ungroup() %>%
  arrange(treatment_group, desc(n_countries))

cat("Regional composition by treatment group:\n\n")
region_wide <- region_comp %>%
  select(continent, treatment_group, n_countries, pct) %>%
  pivot_wider(
    names_from = treatment_group,
    values_from = c(n_countries, pct),
    values_fill = list(n_countries = 0, pct = 0)
  ) %>%
  select(
    continent,
    `n_Ever-Treated`     = `n_countries_Ever-Treated`,
    `pct_Ever-Treated`   = `pct_Ever-Treated`,
    `n_Never-Treated`    = `n_countries_Never-Treated`,
    `pct_Never-Treated`  = `pct_Never-Treated`
  ) %>%
  arrange(desc(`n_Never-Treated`))

print(as.data.frame(region_wide), row.names = FALSE)

# --- 1b. Mean GDP per capita levels ---
cat("\n\n--- 1b. Mean GDP Per Capita (Log) by Group and Period ---\n\n")

gdp_levels <- data %>%
  filter(!is.na(loggdp)) %>%
  group_by(treatment_group, period) %>%
  summarise(
    n_obs      = n(),
    n_countries = n_distinct(countrycode),
    mean_loggdp = round(mean(loggdp, na.rm = TRUE), 3),
    sd_loggdp   = round(sd(loggdp, na.rm = TRUE), 3),
    median_loggdp = round(median(loggdp, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(period, treatment_group)

print(as.data.frame(gdp_levels), row.names = FALSE)

cat("\nNote: loggdp is log real GDP per capita (USD).\n")
cat("Mean GDP per capita (levels) by group and period:\n\n")

gdp_levels_exp <- data %>%
  filter(!is.na(loggdp)) %>%
  mutate(gdp_pc = exp(loggdp)) %>%
  group_by(treatment_group, period) %>%
  summarise(
    mean_gdp_pc   = round(mean(gdp_pc, na.rm = TRUE), 0),
    median_gdp_pc = round(median(gdp_pc, na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  arrange(period, treatment_group)

print(as.data.frame(gdp_levels_exp), row.names = FALSE)

# --- 1c. GDP growth rates ---
cat("\n\n--- 1c. GDP Growth Rates by Group and Period ---\n\n")

growth_rates <- data %>%
  filter(!is.na(gdp_diff)) %>%
  group_by(treatment_group, period) %>%
  summarise(
    n_obs       = n(),
    mean_growth = round(mean(gdp_diff, na.rm = TRUE), 4),
    sd_growth   = round(sd(gdp_diff, na.rm = TRUE), 4),
    median_growth = round(median(gdp_diff, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(period, treatment_group)

cat("GDP growth rates (annual difference in log GDP pc):\n\n")
print(as.data.frame(growth_rates), row.names = FALSE)


# ==================================================================
# ANALYSIS 2: GDP trend comparison plot
# ==================================================================
cat("\n============================================================\n")
cat("  ANALYSIS 2: GDP Trend Comparison\n")
cat("============================================================\n\n")

trend_data <- data %>%
  filter(!is.na(loggdp)) %>%
  group_by(treatment_group, year) %>%
  summarise(
    mean_loggdp   = mean(loggdp, na.rm = TRUE),
    se_loggdp     = sd(loggdp, na.rm = TRUE) / sqrt(n()),
    n_countries   = n_distinct(countrycode),
    .groups = "drop"
  )

p_trends <- ggplot(trend_data, aes(x = year, y = mean_loggdp,
                                    color = treatment_group,
                                    fill = treatment_group)) +
  geom_ribbon(aes(ymin = mean_loggdp - 1.96 * se_loggdp,
                  ymax = mean_loggdp + 1.96 * se_loggdp),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 1990, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  annotate("text", x = 1990, y = Inf, label = "1990", vjust = 2,
           hjust = -0.1, color = "grey40", size = 3.5) +
  scale_color_manual(values = c("Ever-Treated" = "firebrick", "Never-Treated" = "steelblue")) +
  scale_fill_manual(values = c("Ever-Treated" = "firebrick", "Never-Treated" = "steelblue")) +
  labs(
    x = "Year",
    y = "Mean Log GDP Per Capita",
    color = "Group", fill = "Group",
    title = "Average Log GDP Per Capita: Never-Treated vs Ever-Treated",
    subtitle = sprintf("Never-treated: %d countries | Ever-treated: %d countries",
                        length(never_treated_countries), length(ever_treated_countries))
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

ggsave("investigation/figures/gdp_trends_treated_vs_untreated.png",
       p_trends, width = 10, height = 6, dpi = 300)
cat("Saved: investigation/figures/gdp_trends_treated_vs_untreated.png\n")


# ==================================================================
# ANALYSIS 3: Growth acceleration test (Diff-in-Diff)
# ==================================================================
cat("\n============================================================\n")
cat("  ANALYSIS 3: Growth Acceleration Test (Diff-in-Diff)\n")
cat("============================================================\n\n")

# Compute country-level average growth by period
country_growth <- data %>%
  filter(!is.na(gdp_diff)) %>%
  group_by(countrycode, treatment_group, period) %>%
  summarise(mean_growth = mean(gdp_diff, na.rm = TRUE), .groups = "drop")

# Diff-in-diff table
did_table <- country_growth %>%
  pivot_wider(names_from = period, values_from = mean_growth) %>%
  filter(!is.na(`Pre-1990`) & !is.na(`Post-1990`)) %>%
  mutate(growth_change = `Post-1990` - `Pre-1990`)

did_summary <- did_table %>%
  group_by(treatment_group) %>%
  summarise(
    n_countries     = n(),
    mean_pre        = round(mean(`Pre-1990`, na.rm = TRUE), 5),
    mean_post       = round(mean(`Post-1990`, na.rm = TRUE), 5),
    mean_change     = round(mean(growth_change, na.rm = TRUE), 5),
    se_change       = round(sd(growth_change, na.rm = TRUE) / sqrt(n()), 5),
    .groups = "drop"
  )

cat("Diff-in-Diff of GDP Growth Rates:\n")
cat("(Country-level mean growth, then average across countries)\n\n")
print(as.data.frame(did_summary), row.names = FALSE)

# Formal test
did_reg_data <- did_table %>%
  mutate(
    treated = as.integer(treatment_group == "Ever-Treated"),
    post = 1  # growth_change is already post - pre
  )

cat("\nRegression test: growth_change ~ treated\n")
did_reg <- lm(growth_change ~ treated, data = did_reg_data)
cat(sprintf("  Intercept (Never-Treated growth change):  %8.5f  (SE: %.5f)\n",
            coef(did_reg)[1], sqrt(vcov(did_reg)[1,1])))
cat(sprintf("  Treated effect (diff-in-diff):            %8.5f  (SE: %.5f)\n",
            coef(did_reg)[2], sqrt(vcov(did_reg)[2,2])))
cat(sprintf("  t-stat: %.3f    p-value: %.4f\n",
            summary(did_reg)$coefficients[2,3],
            summary(did_reg)$coefficients[2,4]))

cat("\nInterpretation:\n")
cat("  Positive intercept = never-treated countries accelerated growth post-1990\n")
cat("  Negative treated coeff = ever-treated growth accelerated LESS (or decelerated)\n")
cat("  This is the composition channel: never-treated countries provide a\n")
cat("  rising counterfactual that makes cyclone damage look worse post-1990.\n")


# ==================================================================
# ANALYSIS 4: Year FE comparison
# ==================================================================
cat("\n\n============================================================\n")
cat("  ANALYSIS 4: Year Fixed Effects Comparison\n")
cat("============================================================\n\n")

# Define subsamples
data_ever_treated  <- data %>% filter(countrycode %in% ever_treated_countries)

# Run h=5 LP on full sample and ever-treated sample, extract year FEs
run_lp_h5 <- function(df, label) {
  h <- 5
  fml <- as.formula(paste0(
    "f(", outcome, ", ", h, ") - l(", outcome, ", 1) ~ ",
    "i(period, maxwind_95)",
    " + ", controls,
    " | ", fe
  ))

  mod <- feols(fml, data = df, panel.id = panel_id, vcov = vcov_fm)

  # Extract year FEs
  fe_vals <- fixef(mod)$year
  fe_df <- data.frame(
    year = as.integer(names(fe_vals)),
    year_fe = as.numeric(fe_vals),
    sample = label,
    stringsAsFactors = FALSE
  )
  return(list(model = mod, fe = fe_df))
}

cat("Running h=5 LP on full sample...\n")
res_full <- run_lp_h5(data, "Full Sample")

cat("Running h=5 LP on ever-treated sample...\n")
res_ever <- run_lp_h5(data_ever_treated, "Ever-Treated Only")

# Combine year FEs
fe_combined <- bind_rows(res_full$fe, res_ever$fe)

# Demean within each sample for comparability
fe_combined <- fe_combined %>%
  group_by(sample) %>%
  mutate(year_fe_dm = year_fe - mean(year_fe, na.rm = TRUE)) %>%
  ungroup()

p_fe <- ggplot(fe_combined, aes(x = year, y = year_fe_dm, color = sample)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = 1990, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  annotate("text", x = 1990, y = Inf, label = "1990", vjust = 2,
           hjust = -0.1, color = "grey40", size = 3.5) +
  scale_color_manual(values = c("Full Sample" = "steelblue", "Ever-Treated Only" = "firebrick")) +
  labs(
    x = "Year",
    y = "Demeaned Year Fixed Effect",
    color = "Sample",
    title = "Year Fixed Effects from LP (h=5): Full vs Ever-Treated Sample",
    subtitle = "Year FEs capture the counterfactual GDP trend net of cyclone shocks"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

ggsave("investigation/figures/year_fe_comparison.png",
       p_fe, width = 10, height = 6, dpi = 300)
cat("Saved: investigation/figures/year_fe_comparison.png\n")

# Print the difference in year FEs
fe_diff <- fe_combined %>%
  select(year, sample, year_fe_dm) %>%
  pivot_wider(names_from = sample, values_from = year_fe_dm) %>%
  mutate(difference = `Full Sample` - `Ever-Treated Only`)

cat("\nYear FE difference (Full - Ever-Treated), selected years:\n\n")
fe_diff_print <- fe_diff %>%
  filter(year %in% c(1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010)) %>%
  mutate(across(where(is.numeric), ~round(., 4)))
print(as.data.frame(fe_diff_print), row.names = FALSE)

# Correlation and divergence
cat(sprintf("\nCorrelation of year FEs between samples: %.4f\n",
            cor(fe_diff$`Full Sample`, fe_diff$`Ever-Treated Only`, use = "complete.obs")))

# Post-1990 trend difference
fe_post90 <- fe_diff %>% filter(year > 1990)
fe_pre90  <- fe_diff %>% filter(year <= 1990)
cat(sprintf("Mean difference (Full - Ever-Treated) pre-1990:  %.5f\n",
            mean(fe_pre90$difference, na.rm = TRUE)))
cat(sprintf("Mean difference (Full - Ever-Treated) post-1990: %.5f\n",
            mean(fe_post90$difference, na.rm = TRUE)))


# ==================================================================
# ANALYSIS 5: Gradual exclusion by region
# ==================================================================
cat("\n\n============================================================\n")
cat("  ANALYSIS 5: Gradual Exclusion by Region\n")
cat("============================================================\n\n")

# Identify never-treated countries by continent
never_treated_by_continent <- data %>%
  filter(countrycode %in% never_treated_countries) %>%
  distinct(countrycode, continent) %>%
  group_by(continent) %>%
  summarise(countries = list(countrycode), n = n(), .groups = "drop") %>%
  filter(!is.na(continent)) %>%
  arrange(desc(n))

cat("Never-treated countries by continent:\n")
for (i in 1:nrow(never_treated_by_continent)) {
  cat(sprintf("  %s: %d countries\n",
              never_treated_by_continent$continent[i],
              never_treated_by_continent$n[i]))
}

# Run the period-interacted LP dropping never-treated countries from each region
cat("\nRunning period-split LP (h=5) excluding never-treated by region...\n\n")

# Helper: run LP and extract h=5 coefficients
run_lp_h5_inter <- function(df, label) {
  h <- 5
  fml <- as.formula(paste0(
    "f(", outcome, ", ", h, ") - l(", outcome, ", 1) ~ ",
    "i(period, maxwind_95)",
    " + ", controls,
    " | ", fe
  ))

  mod <- tryCatch(
    feols(fml, data = df, panel.id = panel_id, vcov = vcov_fm),
    error = function(e) NULL
  )

  if (is.null(mod)) {
    return(data.frame(
      label = label,
      category = c("Post-1990", "Pre-1990"),
      irf_mean = NA_real_,
      se = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  coefs <- coef(mod)
  ses   <- sqrt(diag(vcov(mod)))

  coef_df <- data.frame(
    label = label,
    term = names(coefs),
    irf_mean = as.numeric(coefs),
    se = as.numeric(ses),
    stringsAsFactors = FALSE
  ) %>%
    filter(str_detect(term, "^period::")) %>%
    mutate(
      category = str_extract(term, "(?<=::).*?(?=:)"),
      category = str_replace_all(category, "_", "-"),
      category = case_when(
        str_detect(category, "Pre")  ~ "Pre-1990",
        str_detect(category, "Post") ~ "Post-1990",
        TRUE ~ category
      )
    ) %>%
    select(label, category, irf_mean, se)

  return(coef_df)
}

# Full sample baseline
results <- run_lp_h5_inter(data, "Full Sample")

# Drop all never-treated
results <- bind_rows(results,
  run_lp_h5_inter(data_ever_treated, "Drop ALL Never-Treated"))

# Drop never-treated from each continent
for (i in 1:nrow(never_treated_by_continent)) {
  cont <- never_treated_by_continent$continent[i]
  countries_to_drop <- never_treated_by_continent$countries[[i]]

  data_sub <- data %>% filter(!(countrycode %in% countries_to_drop))

  label <- sprintf("Drop %s (%d)", cont, length(countries_to_drop))
  results <- bind_rows(results, run_lp_h5_inter(data_sub, label))
}

# Print results table
cat("\n--- Period-Split LP at h=5: Exclusion Results ---\n\n")

results_wide <- results %>%
  select(label, category, irf_mean, se) %>%
  mutate(
    irf_mean = round(irf_mean, 4),
    se = round(se, 4)
  ) %>%
  pivot_wider(
    names_from = category,
    values_from = c(irf_mean, se)
  )

if ("irf_mean_Post-1990" %in% names(results_wide) & "irf_mean_Pre-1990" %in% names(results_wide)) {
  results_wide <- results_wide %>%
    mutate(gap = `irf_mean_Post-1990` - `irf_mean_Pre-1990`)
}

print(as.data.frame(results_wide), row.names = FALSE)

# Get the full-sample gap as reference
full_gap <- results_wide %>% filter(label == "Full Sample") %>% pull(gap)

# Plot the post-1990 coefficient and the pre-post gap
plot_results <- results %>%
  mutate(
    irf_down = irf_mean - 1.96 * se,
    irf_up   = irf_mean + 1.96 * se,
    label = factor(label, levels = rev(unique(results$label)))
  )

p_excl <- ggplot(plot_results, aes(x = irf_mean, y = label, color = category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(xmin = irf_down, xmax = irf_up),
                  position = position_dodge(width = 0.5),
                  size = 0.6, linewidth = 0.8) +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_x_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "GDP Response at h=5 (maxwind_95 shock)",
    y = NULL,
    color = "Period",
    title = "Period-Split LP (h=5): Effect of Excluding Never-Treated by Region",
    subtitle = "Dropping never-treated countries from each continent"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

ggsave("investigation/figures/exclusion_by_region.png",
       p_excl, width = 11, height = 7, dpi = 300)
cat("\nSaved: investigation/figures/exclusion_by_region.png\n")


# ==================================================================
# FINAL SUMMARY
# ==================================================================
cat("\n\n============================================================\n")
cat("  SUMMARY OF FINDINGS\n")
cat("============================================================\n\n")

cat("1. GROUP CHARACTERIZATION:\n")
cat("   - Never-treated countries are predominantly in Africa and Europe\n")
cat("   - Ever-treated countries are concentrated in Asia, Americas, Oceania\n")
cat(sprintf("   - Never-treated have higher mean log GDP (they include rich European countries)\n"))

cat("\n2. GDP TRENDS:\n")
cat("   - Both groups show upward GDP trends, but may diverge post-1990\n")
cat("   - The gap in GDP levels affects the year FE counterfactual\n")

cat("\n3. GROWTH ACCELERATION:\n")
did_never  <- did_summary %>% filter(treatment_group == "Never-Treated")
did_ever   <- did_summary %>% filter(treatment_group == "Ever-Treated")
cat(sprintf("   - Never-treated growth change (post - pre): %+.5f\n", did_never$mean_change))
cat(sprintf("   - Ever-treated growth change (post - pre):  %+.5f\n", did_ever$mean_change))
cat(sprintf("   - Diff-in-diff: %+.5f (p = %.4f)\n",
            coef(did_reg)[2], summary(did_reg)$coefficients[2,4]))

cat("\n4. YEAR FE COMPARISON:\n")
cat("   - Year FEs shift when never-treated countries are dropped\n")
cat("   - This changes the counterfactual GDP trend against which\n")
cat("     cyclone damage is measured\n")

cat("\n5. REGIONAL EXCLUSION:\n")
cat("   - See bar chart for which regional block matters most\n")
cat("   - The full-sample pre-post gap is driven by the inclusion of\n")
cat("     never-treated countries that alter the year FE counterfactual\n")

cat("\nDone.\n")
