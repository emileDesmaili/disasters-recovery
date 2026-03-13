# investigation_never_treated.R
# Robustness check: what happens when we drop "never-treated" countries
# (countries that never experience a p95 cyclone shock) from the estimation sample?
#
# Key question: Does the inclusion of countries with no variation in treatment
# affect the estimated IRFs, particularly the pre/post-1990 divergence?

library(tidyverse)
library(fixest)
library(haven)
source("../emileRegs.R")

setFixest_notes(FALSE)

cat("============================================================\n")
cat("  NEVER-TREATED COUNTRIES INVESTIGATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------------
# 1. Load and construct data (identically to other investigation scripts)
# ------------------------------------------------------------------
pwt  <- read_stata("../raw_data/pwt_clean.dta")
tcs  <- read_stata("../raw_data/ibtracs_clean.dta")
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

# ==================================================================
# ANALYSIS 1: Identify never-treated countries
# ==================================================================
cat("============================================================\n")
cat("  ANALYSIS 1: Identify Never-Treated Countries\n")
cat("============================================================\n\n")

# Countries that EVER have maxwind_95 == 1
ever_treated_countries <- data %>%
  filter(maxwind_95 == 1) %>%
  pull(countrycode) %>%
  unique() %>%
  sort()

# Countries that NEVER have maxwind_95 == 1
all_countries <- unique(data$countrycode)
never_treated_countries <- sort(setdiff(all_countries, ever_treated_countries))

cat(sprintf("Total countries in sample:     %d\n", length(all_countries)))
cat(sprintf("Ever-treated (>= 1 shock):     %d\n", length(ever_treated_countries)))
cat(sprintf("Never-treated (0 shocks):      %d\n", length(never_treated_countries)))

cat(sprintf("\nNever-treated countries (%d):\n", length(never_treated_countries)))
cat(paste("  ", never_treated_countries, collapse = "\n"), "\n")

# Shock counts for ever-treated countries
shock_counts <- data %>%
  filter(maxwind_95 == 1) %>%
  group_by(countrycode) %>%
  summarise(
    n_shocks = n(),
    first_shock = min(year),
    last_shock  = max(year),
    .groups = "drop"
  ) %>%
  arrange(desc(n_shocks))

cat(sprintf("\nEver-treated countries with shock counts (%d):\n", nrow(shock_counts)))
print(as.data.frame(shock_counts), row.names = FALSE)

# Which countries have shocks in which periods?
shock_by_period <- data %>%
  filter(maxwind_95 == 1) %>%
  mutate(period = ifelse(year <= 1990, "Pre-1990", "Post-1990")) %>%
  group_by(countrycode, period) %>%
  summarise(n_shocks = n(), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = n_shocks, values_fill = 0)

countries_pre_only  <- shock_by_period %>%
  filter(`Pre-1990` > 0, `Post-1990` == 0) %>% pull(countrycode)
countries_post_only <- shock_by_period %>%
  filter(`Pre-1990` == 0, `Post-1990` > 0) %>% pull(countrycode)
countries_both      <- shock_by_period %>%
  filter(`Pre-1990` > 0, `Post-1990` > 0) %>% pull(countrycode)

cat(sprintf("\nShock distribution across periods (among ever-treated):\n"))
cat(sprintf("  Shocks in Pre-1990 ONLY:   %d countries\n", length(countries_pre_only)))
cat(sprintf("  Shocks in Post-1990 ONLY:  %d countries\n", length(countries_post_only)))
cat(sprintf("  Shocks in BOTH periods:    %d countries\n", length(countries_both)))

# ------------------------------------------------------------------
# Define subsamples
# ------------------------------------------------------------------
data_ever_treated  <- data %>% filter(countrycode %in% ever_treated_countries)
data_both_periods  <- data %>% filter(countrycode %in% countries_both)

cat(sprintf("\nSample sizes:\n"))
cat(sprintf("  Full sample:       %d obs, %d countries\n",
            nrow(data), n_distinct(data$countrycode)))
cat(sprintf("  Ever-treated:      %d obs, %d countries\n",
            nrow(data_ever_treated), n_distinct(data_ever_treated$countrycode)))
cat(sprintf("  Both-period shocks: %d obs, %d countries\n",
            nrow(data_both_periods), n_distinct(data_both_periods$countrycode)))

# ==================================================================
# ANALYSIS 2: Period-interacted LP -- full vs ever-treated
# ==================================================================
cat("\n============================================================\n")
cat("  ANALYSIS 2: Period-Interacted LP (Full vs Ever-Treated)\n")
cat("============================================================\n\n")

cat("Running period-interacted LP on full sample...\n")
irf_full_inter <- lp_panel_inter(
  data = data, outcome = outcome, main_var = "maxwind_95",
  interact_var = "period", controls = controls,
  horizon = horizon, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
)

cat("Running period-interacted LP on ever-treated sample...\n")
irf_ever_inter <- lp_panel_inter(
  data = data_ever_treated, outcome = outcome, main_var = "maxwind_95",
  interact_var = "period", controls = controls,
  horizon = horizon, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
)

# Clean up category names
clean_categories <- function(df) {
  df %>%
    mutate(
      category = str_replace_all(category, "_", "-"),
      category = case_when(
        str_detect(category, "Pre")  ~ "Pre-1990",
        str_detect(category, "Post") ~ "Post-1990",
        TRUE ~ category
      )
    )
}

irf_full_inter  <- clean_categories(irf_full_inter)
irf_ever_inter  <- clean_categories(irf_ever_inter)

# --- Plot: full vs ever-treated, period-interacted ---
irf_full_inter$sample  <- "Full Sample"
irf_ever_inter$sample  <- "Ever-Treated Only"

plot_data_2 <- bind_rows(irf_full_inter, irf_ever_inter)

p2 <- ggplot(plot_data_2, aes(x = horizon, y = irf_mean, color = category, linetype = sample)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up, fill = category),
              alpha = 0.08, color = NA,
              data = plot_data_2 %>% filter(sample == "Full Sample")) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey40") +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_fill_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_linetype_manual(values = c("Full Sample" = "solid", "Ever-Treated Only" = "dashed")) +
  scale_x_continuous(breaks = 0:10) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "Horizon (years)",
    y = "GDP response to maxwind_95 shock",
    color = "Period", fill = "Period", linetype = "Sample",
    title = "Period-Interacted IRFs: Full Sample vs Ever-Treated Only",
    subtitle = sprintf("Dropping %d never-treated countries (solid = full, dashed = ever-treated only)",
                        length(never_treated_countries))
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

ggsave("figures/irf_drop_never_treated.png", p2, width = 12, height = 6, dpi = 300)
cat("Saved: figures/irf_drop_never_treated.png\n")

# ==================================================================
# ANALYSIS 3: Baseline LP (no period interaction) -- full vs ever-treated
# ==================================================================
cat("\n============================================================\n")
cat("  ANALYSIS 3: Baseline LP (No Period Interaction)\n")
cat("============================================================\n\n")

cat("Running baseline LP on full sample...\n")
irf_full_base <- lp_panel(
  data = data, outcome = outcome, main_var = "maxwind_95",
  controls = controls, horizon = horizon,
  fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
)

cat("Running baseline LP on ever-treated sample...\n")
irf_ever_base <- lp_panel(
  data = data_ever_treated, outcome = outcome, main_var = "maxwind_95",
  controls = controls, horizon = horizon,
  fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
)

irf_full_base$sample <- "Full Sample"
irf_ever_base$sample <- "Ever-Treated Only"

plot_data_3 <- bind_rows(irf_full_base, irf_ever_base)

p3 <- ggplot(plot_data_3, aes(x = horizon, y = irf_mean, linetype = sample)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up),
              alpha = 0.10, fill = "grey50", color = NA,
              data = plot_data_3 %>% filter(sample == "Full Sample")) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up),
              alpha = 0.08, fill = "darkorange", color = NA,
              data = plot_data_3 %>% filter(sample == "Ever-Treated Only")) +
  geom_line(linewidth = 1.1, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey40") +
  scale_linetype_manual(values = c("Full Sample" = "solid", "Ever-Treated Only" = "dashed")) +
  scale_x_continuous(breaks = 0:10) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "Horizon (years)",
    y = "GDP response to maxwind_95 shock",
    linetype = "Sample",
    title = "Baseline IRF: Full Sample vs Ever-Treated Only",
    subtitle = sprintf("No period interaction | Dropping %d never-treated countries",
                        length(never_treated_countries))
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

ggsave("figures/irf_baseline_drop_never_treated.png", p3, width = 10, height = 6, dpi = 300)
cat("Saved: figures/irf_baseline_drop_never_treated.png\n")

# ==================================================================
# ANALYSIS 4: Strictest sample -- only countries with shocks in BOTH periods
# ==================================================================
cat("\n============================================================\n")
cat("  ANALYSIS 4: Strictest Sample -- Shocks in BOTH Periods\n")
cat("============================================================\n\n")

cat(sprintf("Countries with shocks in both periods: %d\n", length(countries_both)))
cat(paste("  ", sort(countries_both), collapse = "\n"), "\n")

cat("\nRunning period-interacted LP on both-period-treated sample...\n")
irf_both_inter <- lp_panel_inter(
  data = data_both_periods, outcome = outcome, main_var = "maxwind_95",
  interact_var = "period", controls = controls,
  horizon = horizon, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
)
irf_both_inter <- clean_categories(irf_both_inter)
irf_both_inter$sample <- "Both-Period Treated"

# Combine all three samples
plot_data_4 <- bind_rows(
  irf_full_inter,
  irf_ever_inter,
  irf_both_inter
) %>%
  mutate(sample = factor(sample, levels = c("Full Sample", "Ever-Treated Only", "Both-Period Treated")))

p4 <- ggplot(plot_data_4, aes(x = horizon, y = irf_mean, color = category, linetype = sample)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up, fill = category),
              alpha = 0.06, color = NA,
              data = plot_data_4 %>% filter(sample == "Full Sample")) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey40") +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_fill_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_linetype_manual(values = c("Full Sample" = "solid",
                                    "Ever-Treated Only" = "dashed",
                                    "Both-Period Treated" = "dotted")) +
  scale_x_continuous(breaks = 0:10) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "Horizon (years)",
    y = "GDP response to maxwind_95 shock",
    color = "Period", fill = "Period", linetype = "Sample",
    title = "Period-Interacted IRFs: Varying Sample Restrictions",
    subtitle = sprintf("Full (%d countries) vs Ever-Treated (%d) vs Both-Period Treated (%d)",
                        n_distinct(data$countrycode),
                        length(ever_treated_countries),
                        length(countries_both))
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top") +
  guides(
    color = guide_legend(order = 1),
    fill  = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )

ggsave("figures/irf_both_period_treated.png", p4, width = 12, height = 6, dpi = 300)
cat("Saved: figures/irf_both_period_treated.png\n")

# ==================================================================
# SUMMARY TABLES
# ==================================================================
cat("\n============================================================\n")
cat("  SUMMARY TABLES: Coefficient Comparison Across Samples\n")
cat("============================================================\n\n")

# Collect all period-interacted results
all_inter <- bind_rows(
  irf_full_inter %>% mutate(sample = "Full Sample"),
  irf_ever_inter %>% mutate(sample = "Ever-Treated Only"),
  irf_both_inter %>% mutate(sample = "Both-Period Treated")
)

# --- h = 5 ---
cat("--- Period-Interacted IRFs at h = 5 ---\n\n")
summary_h5 <- all_inter %>%
  filter(horizon == 5) %>%
  select(sample, category, irf_mean, se, irf_down, irf_up) %>%
  arrange(category, sample)

print(as.data.frame(summary_h5), row.names = FALSE)

# --- h = 10 ---
cat("\n--- Period-Interacted IRFs at h = 10 ---\n\n")
summary_h10 <- all_inter %>%
  filter(horizon == 10) %>%
  select(sample, category, irf_mean, se, irf_down, irf_up) %>%
  arrange(category, sample)

print(as.data.frame(summary_h10), row.names = FALSE)

# --- Baseline (no interaction) at h=5 and h=10 ---
all_base <- bind_rows(
  irf_full_base %>% mutate(sample = "Full Sample"),
  irf_ever_base %>% mutate(sample = "Ever-Treated Only")
)

cat("\n--- Baseline (No Interaction) IRFs at h = 5 ---\n\n")
summary_base_h5 <- all_base %>%
  filter(horizon == 5) %>%
  select(sample, irf_mean, se, irf_down, irf_up)

print(as.data.frame(summary_base_h5), row.names = FALSE)

cat("\n--- Baseline (No Interaction) IRFs at h = 10 ---\n\n")
summary_base_h10 <- all_base %>%
  filter(horizon == 10) %>%
  select(sample, irf_mean, se, irf_down, irf_up)

print(as.data.frame(summary_base_h10), row.names = FALSE)

# ==================================================================
# VERDICT
# ==================================================================
cat("\n============================================================\n")
cat("  VERDICT\n")
cat("============================================================\n\n")

# Compare the pre-post gap across samples
gap_table <- all_inter %>%
  filter(horizon %in% c(5, 10)) %>%
  select(sample, category, horizon, irf_mean) %>%
  pivot_wider(names_from = category, values_from = irf_mean) %>%
  mutate(gap = `Post-1990` - `Pre-1990`)

cat("Pre-Post gap by sample and horizon:\n")
print(as.data.frame(gap_table), row.names = FALSE)

cat("\nInterpretation:\n")
cat("- If dropping never-treated countries substantially changes the IRFs,\n")
cat("  the control group of unexposed countries matters for identification.\n")
cat("- If results are stable, the estimates are robust to this sample restriction.\n")

# Check stability
full_h5_post  <- all_inter %>% filter(sample == "Full Sample", category == "Post-1990", horizon == 5) %>% pull(irf_mean)
ever_h5_post  <- all_inter %>% filter(sample == "Ever-Treated Only", category == "Post-1990", horizon == 5) %>% pull(irf_mean)
both_h5_post  <- all_inter %>% filter(sample == "Both-Period Treated", category == "Post-1990", horizon == 5) %>% pull(irf_mean)

if (length(full_h5_post) > 0 & length(ever_h5_post) > 0) {
  pct_change_ever <- (ever_h5_post - full_h5_post) / abs(full_h5_post) * 100
  cat(sprintf("\nPost-1990 h=5 coefficient change (Full -> Ever-Treated): %.1f%%\n", pct_change_ever))
}
if (length(full_h5_post) > 0 & length(both_h5_post) > 0) {
  pct_change_both <- (both_h5_post - full_h5_post) / abs(full_h5_post) * 100
  cat(sprintf("Post-1990 h=5 coefficient change (Full -> Both-Period):  %.1f%%\n", pct_change_both))
}

cat("\nDone.\n")
