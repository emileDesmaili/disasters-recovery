# investigation_influential.R
# Leave-one-out influence analysis for the period-interacted LP results
# Question: Are a few extreme disaster years or countries driving the post-1990 result?

library(tidyverse)
library(fixest)
library(haven)
source("../emileRegs.R")

setFixest_notes(FALSE)

# ------------------------------------------------------------------
# 1. Load and construct data (exactly as in lp_gdp_shocks.R)
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
outcome      <- "loggdp"
shock_var    <- "maxwind_95"
controls     <- "l(gdp_diff,1:2) + l(maxwind_95,1:2)"
fe           <- "countrycode[year] + countrycode[year2] + year"
panel_id     <- c("countrycode", "year")
vcov_fm      <- DK ~ year
horizon      <- 10
target_h     <- 5

# ------------------------------------------------------------------
# Helper: run the period-interacted LP and extract h=target_h coefficients
# ------------------------------------------------------------------
run_lp_extract <- function(df) {
  irf <- lp_panel_inter(
    data         = df,
    outcome      = outcome,
    main_var     = shock_var,
    interact_var = "period",
    controls     = controls,
    horizon      = horizon,
    fe           = fe,
    panel_id     = panel_id,
    vcov_formula = vcov_fm
  )
  irf %>%
    filter(horizon == target_h) %>%
    select(category, irf_mean, se, irf_down, irf_up)
}

# ------------------------------------------------------------------
# Full-sample baseline
# ------------------------------------------------------------------
cat("Running full-sample baseline estimation...\n")
baseline <- run_lp_extract(data)
cat("Full-sample h=5 coefficients:\n")
print(baseline)

baseline_post <- baseline %>% filter(category == "Post-1990") %>% pull(irf_mean)
baseline_pre  <- baseline %>% filter(category == "Pre-1990")  %>% pull(irf_mean)

# ==================================================================
# PART 1: Leave-one-year-out
# ==================================================================
cat("\n========================================\n")
cat("PART 1: Leave-one-year-out analysis\n")
cat("========================================\n")

years_to_drop <- 1970:2014

loo_year_results <- map_dfr(years_to_drop, function(yr) {
  cat(sprintf("\r  Dropping year %d...", yr))
  df_sub <- data %>% filter(year != yr)
  tryCatch({
    res <- run_lp_extract(df_sub)
    res$dropped_year <- yr
    res
  }, error = function(e) {
    tibble(
      category = c("Pre-1990", "Post-1990"),
      irf_mean = NA_real_, se = NA_real_,
      irf_down = NA_real_, irf_up = NA_real_,
      dropped_year = yr
    )
  })
})
cat("\n")

# Compute influence: absolute change in Post-1990 coefficient when year is dropped
loo_year_influence <- loo_year_results %>%
  filter(category == "Post-1990") %>%
  mutate(
    influence = abs(irf_mean - baseline_post),
    direction = ifelse(irf_mean > baseline_post, "more negative (dampens)", "less negative (amplifies)")
  ) %>%
  arrange(desc(influence))

top5_years <- head(loo_year_influence, 5)

cat("\nTop 5 most influential years (Post-1990 coefficient):\n")
cat(sprintf("  Full-sample Post-1990 h=5 coef: %.4f\n", baseline_post))
for (i in 1:nrow(top5_years)) {
  cat(sprintf("  %d. Drop %d -> coef = %.4f (change = %.4f)\n",
              i, top5_years$dropped_year[i], top5_years$irf_mean[i], top5_years$influence[i]))
}

# Plot: LOO year influence
loo_year_plot_data <- loo_year_results %>%
  filter(!is.na(irf_mean)) %>%
  mutate(category = factor(category, levels = c("Pre-1990", "Post-1990")))

baseline_lines <- tibble(
  category = factor(c("Pre-1990", "Post-1990"), levels = c("Pre-1990", "Post-1990")),
  baseline = c(baseline_pre, baseline_post)
)

p1 <- ggplot(loo_year_plot_data, aes(x = dropped_year, y = irf_mean, color = category)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(data = baseline_lines, aes(yintercept = baseline, color = category),
             linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  labs(
    title = "Leave-One-Year-Out: h=5 IRF Coefficient",
    subtitle = paste0("Shock: maxwind_95 | Dashed lines = full-sample estimates"),
    x = "Year Dropped",
    y = "h=5 IRF Coefficient",
    color = "Period"
  ) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

ggsave("figures/loo_year_influence.png", p1, width = 12, height = 5, dpi = 300)
cat("  Saved: figures/loo_year_influence.png\n")

# ==================================================================
# PART 2: Leave-one-country-out
# ==================================================================
cat("\n========================================\n")
cat("PART 2: Leave-one-country-out analysis\n")
cat("========================================\n")

# Identify countries with at least one maxwind_95 = 1 event
countries_with_events <- data %>%
  filter(maxwind_95 == 1) %>%
  pull(countrycode) %>%
  unique() %>%
  sort()

cat(sprintf("  Countries with at least one p95 shock: %d\n", length(countries_with_events)))

loo_country_results <- map_dfr(countries_with_events, function(cc) {
  cat(sprintf("\r  Dropping country %s...     ", cc))
  df_sub <- data %>% filter(countrycode != cc)
  tryCatch({
    res <- run_lp_extract(df_sub)
    res$dropped_country <- cc
    res
  }, error = function(e) {
    tibble(
      category = c("Pre-1990", "Post-1990"),
      irf_mean = NA_real_, se = NA_real_,
      irf_down = NA_real_, irf_up = NA_real_,
      dropped_country = cc
    )
  })
})
cat("\n")

# Compute influence
loo_country_influence <- loo_country_results %>%
  filter(category == "Post-1990") %>%
  mutate(
    influence = abs(irf_mean - baseline_post)
  ) %>%
  arrange(desc(influence))

top5_countries <- head(loo_country_influence, 5)

cat("\nTop 5 most influential countries (Post-1990 coefficient):\n")
cat(sprintf("  Full-sample Post-1990 h=5 coef: %.4f\n", baseline_post))
for (i in 1:nrow(top5_countries)) {
  cat(sprintf("  %d. Drop %s -> coef = %.4f (change = %.4f)\n",
              i, top5_countries$dropped_country[i], top5_countries$irf_mean[i], top5_countries$influence[i]))
}

# Plot: LOO country influence, sorted by influence on Post-1990
loo_country_plot_data <- loo_country_results %>%
  filter(!is.na(irf_mean)) %>%
  mutate(category = factor(category, levels = c("Pre-1990", "Post-1990")))

# Sort countries by influence on Post-1990
country_order <- loo_country_influence %>%
  arrange(influence) %>%
  pull(dropped_country)

loo_country_plot_data <- loo_country_plot_data %>%
  mutate(dropped_country = factor(dropped_country, levels = country_order))

p2 <- ggplot(loo_country_plot_data, aes(x = dropped_country, y = irf_mean, color = category)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_hline(data = baseline_lines, aes(yintercept = baseline, color = category),
             linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  labs(
    title = "Leave-One-Country-Out: h=5 IRF Coefficient",
    subtitle = paste0("Shock: maxwind_95 | Countries sorted by influence on Post-1990 estimate"),
    x = "Country Dropped",
    y = "h=5 IRF Coefficient",
    color = "Period"
  ) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
  ) +
  coord_flip()

ggsave("figures/loo_country_influence.png", p2, width = 8, height = 10, dpi = 300)
cat("  Saved: figures/loo_country_influence.png\n")

# ==================================================================
# PART 3: Event-year analysis
# ==================================================================
cat("\n========================================\n")
cat("PART 3: Event-year time series\n")
cat("========================================\n")

events_by_year <- data %>%
  filter(year >= 1970, year <= 2014) %>%
  group_by(year) %>%
  summarise(
    n_events = sum(maxwind_95, na.rm = TRUE),
    .groups = "drop"
  )

# Mark the most influential years
influential_years <- top5_years$dropped_year
events_by_year <- events_by_year %>%
  mutate(
    influential = year %in% influential_years,
    label = ifelse(influential, as.character(year), NA_character_)
  )

p3 <- ggplot(events_by_year, aes(x = year, y = n_events)) +
  geom_col(aes(fill = influential), width = 0.7) +
  geom_text(aes(label = label), vjust = -0.5, size = 3.5, fontface = "bold", na.rm = TRUE) +
  geom_vline(xintercept = 1990.5, linetype = "dotted", color = "grey40", linewidth = 0.6) +
  scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick"),
                    labels = c("Other", "Top 5 influential (LOO)"),
                    name = "") +
  labs(
    title = "Number of Countries Experiencing a p95 Max Wind Shock by Year",
    subtitle = "Influential years from leave-one-year-out analysis highlighted",
    x = "Year",
    y = "Number of countries with maxwind_95 = 1"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top") +
  annotate("text", x = 1980, y = max(events_by_year$n_events) * 0.95,
           label = "Pre-1990", color = "steelblue", fontface = "bold", size = 4.5) +
  annotate("text", x = 2002, y = max(events_by_year$n_events) * 0.95,
           label = "Post-1990", color = "firebrick", fontface = "bold", size = 4.5)

ggsave("figures/shock_events_timeline.png", p3, width = 12, height = 5, dpi = 300)
cat("  Saved: figures/shock_events_timeline.png\n")

# ==================================================================
# SUMMARY
# ==================================================================
cat("\n========================================\n")
cat("SUMMARY OF FINDINGS\n")
cat("========================================\n\n")

cat(sprintf("Full-sample h=%d coefficients:\n", target_h))
cat(sprintf("  Pre-1990:  %.4f (SE: %.4f)\n", baseline_pre,
            baseline %>% filter(category == "Pre-1990") %>% pull(se)))
cat(sprintf("  Post-1990: %.4f (SE: %.4f)\n", baseline_post,
            baseline %>% filter(category == "Post-1990") %>% pull(se)))

cat("\n--- Leave-One-Year-Out ---\n")
loo_year_post <- loo_year_results %>% filter(category == "Post-1990", !is.na(irf_mean))
cat(sprintf("  Range of Post-1990 coef when dropping years: [%.4f, %.4f]\n",
            min(loo_year_post$irf_mean), max(loo_year_post$irf_mean)))
cat(sprintf("  SD of Post-1990 coef across LOO runs: %.4f\n", sd(loo_year_post$irf_mean)))
cat(sprintf("  Top 5 influential years: %s\n",
            paste(top5_years$dropped_year, collapse = ", ")))

# Check: does dropping any single year flip the sign?
any_sign_flip_year <- any(sign(loo_year_post$irf_mean) != sign(baseline_post))
cat(sprintf("  Does dropping any single year flip the sign? %s\n",
            ifelse(any_sign_flip_year, "YES", "NO")))

cat("\n--- Leave-One-Country-Out ---\n")
loo_country_post <- loo_country_results %>% filter(category == "Post-1990", !is.na(irf_mean))
cat(sprintf("  Range of Post-1990 coef when dropping countries: [%.4f, %.4f]\n",
            min(loo_country_post$irf_mean), max(loo_country_post$irf_mean)))
cat(sprintf("  SD of Post-1990 coef across LOO runs: %.4f\n", sd(loo_country_post$irf_mean)))
cat(sprintf("  Top 5 influential countries: %s\n",
            paste(top5_countries$dropped_country, collapse = ", ")))

any_sign_flip_country <- any(sign(loo_country_post$irf_mean) != sign(baseline_post))
cat(sprintf("  Does dropping any single country flip the sign? %s\n",
            ifelse(any_sign_flip_country, "YES", "NO")))

cat("\n--- Event Distribution ---\n")
cat(sprintf("  Total p95 shock events 1970-2014: %d\n", sum(events_by_year$n_events)))
cat(sprintf("  Pre-1990 events: %d\n",
            sum(events_by_year$n_events[events_by_year$year <= 1990])))
cat(sprintf("  Post-1990 events: %d\n",
            sum(events_by_year$n_events[events_by_year$year > 1990])))
cat(sprintf("  Mean events/year (Pre-1990): %.1f\n",
            mean(events_by_year$n_events[events_by_year$year <= 1990])))
cat(sprintf("  Mean events/year (Post-1990): %.1f\n",
            mean(events_by_year$n_events[events_by_year$year > 1990])))

cat("\n--- Robustness Assessment ---\n")
# Compute coefficient of variation for the LOO estimates
cv_year <- sd(loo_year_post$irf_mean) / abs(mean(loo_year_post$irf_mean))
cv_country <- sd(loo_country_post$irf_mean) / abs(mean(loo_country_post$irf_mean))
cat(sprintf("  CV of Post-1990 coef (LOO year):    %.3f\n", cv_year))
cat(sprintf("  CV of Post-1990 coef (LOO country): %.3f\n", cv_country))

if (cv_year < 0.1 & cv_country < 0.1) {
  cat("  CONCLUSION: The post-1990 result appears BROADLY BASED.\n")
  cat("  Neither individual years nor countries exert outsized influence.\n")
} else if (cv_year > 0.25 | cv_country > 0.25) {
  cat("  CONCLUSION: The post-1990 result appears FRAGILE.\n")
  cat("  A few observations exert substantial influence on the estimate.\n")
} else {
  cat("  CONCLUSION: The post-1990 result shows MODERATE sensitivity.\n")
  cat("  Some observations are influential but the result is not driven by a single case.\n")
}

cat("\nDone.\n")
