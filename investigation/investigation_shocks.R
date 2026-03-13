# investigation_shocks.R
# Investigate whether cyclone shock characteristics differ pre- vs post-1990,
# potentially explaining divergent local projection IRF estimates.
#
# Author: Investigation script for disasters-recovery project

library(tidyverse)
library(haven)

# ==================================================================
# 1. Load and construct data identically to lp_gdp_shocks.R
# ==================================================================
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

# Global 95th-percentile thresholds (1970-2014)
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

# Standardize continuous measures (unit = 1 SD)
data <- data %>%
  mutate(
    maxwind_sd = maxwind / sd(maxwind[maxwind > 0], na.rm = TRUE),
    energy_sd  = energy  / sd(energy[energy > 0],   na.rm = TRUE),
    nlands_sd  = nlands  / sd(nlands[nlands > 0],   na.rm = TRUE)
  )

cat("\n")
cat("================================================================\n")
cat("  INVESTIGATION: Cyclone Shock Characteristics Pre vs Post 1990\n")
cat("================================================================\n\n")

# ==================================================================
# Basic dataset dimensions
# ==================================================================
cat("--- Dataset Overview ---\n")
cat(sprintf("Total country-years: %d\n", nrow(data)))
cat(sprintf("Year range: %d - %d\n", min(data$year), max(data$year)))
cat(sprintf("Countries: %d\n", n_distinct(data$countrycode)))
cat(sprintf("Pre-1990 country-years:  %d\n", sum(data$period == "Pre-1990")))
cat(sprintf("Post-1990 country-years: %d\n", sum(data$period == "Post-1990")))

# Report the thresholds
cat(sprintf("\nGlobal 95th percentile thresholds (1970-2014):\n"))
cat(sprintf("  maxwind_p95 = %.2f m/s\n", data$maxwind_p95[1]))
cat(sprintf("  energy_p95  = %.2f\n", data$energy_p95[1]))
cat(sprintf("  nlands_p95  = %.2f\n", data$nlands_p95[1]))

# ==================================================================
# 2. Intensity: Mean and median of shock variables (conditional on > 0)
# ==================================================================
cat("\n\n================================================================\n")
cat("  2. SHOCK INTENSITY (conditional on exposure > 0)\n")
cat("================================================================\n\n")

intensity_stats <- data %>%
  filter(maxwind > 0 | energy > 0 | nlands > 0) %>%
  group_by(period) %>%
  summarise(
    n_exposed = n(),
    # Maxwind
    maxwind_mean   = mean(maxwind[maxwind > 0], na.rm = TRUE),
    maxwind_median = median(maxwind[maxwind > 0], na.rm = TRUE),
    maxwind_sd_val = sd(maxwind[maxwind > 0], na.rm = TRUE),
    maxwind_p90    = quantile(maxwind[maxwind > 0], 0.90, na.rm = TRUE),
    maxwind_p95    = quantile(maxwind[maxwind > 0], 0.95, na.rm = TRUE),
    # Energy
    energy_mean    = mean(energy[energy > 0], na.rm = TRUE),
    energy_median  = median(energy[energy > 0], na.rm = TRUE),
    energy_sd_val  = sd(energy[energy > 0], na.rm = TRUE),
    # Nlands
    nlands_mean    = mean(nlands[nlands > 0], na.rm = TRUE),
    nlands_median  = median(nlands[nlands > 0], na.rm = TRUE),
    nlands_sd_val  = sd(nlands[nlands > 0], na.rm = TRUE),
    .groups = "drop"
  )

cat("Variable          | Pre-1990          | Post-1990\n")
cat("------------------+-------------------+-------------------\n")
for (v in c("maxwind", "energy", "nlands")) {
  pre  <- intensity_stats %>% filter(period == "Pre-1990")
  post <- intensity_stats %>% filter(period == "Post-1990")
  cat(sprintf("%-17s | mean=%.2f  med=%.2f | mean=%.2f  med=%.2f\n",
              v,
              pull(pre, paste0(v, "_mean")),
              pull(pre, paste0(v, "_median")),
              pull(post, paste0(v, "_mean")),
              pull(post, paste0(v, "_median"))))
}

cat("\n--- N exposed country-years ---\n")
cat(sprintf("Pre-1990:  %d\n", intensity_stats$n_exposed[intensity_stats$period == "Pre-1990"]))
cat(sprintf("Post-1990: %d\n", intensity_stats$n_exposed[intensity_stats$period == "Post-1990"]))

# Wilcoxon test for intensity differences
cat("\n--- Wilcoxon rank-sum tests (exposed obs only) ---\n")
pre_mw  <- data %>% filter(period == "Pre-1990", maxwind > 0) %>% pull(maxwind)
post_mw <- data %>% filter(period == "Post-1990", maxwind > 0) %>% pull(maxwind)
wt <- wilcox.test(pre_mw, post_mw)
cat(sprintf("maxwind: W = %.0f, p = %.4f\n", wt$statistic, wt$p.value))

pre_en  <- data %>% filter(period == "Pre-1990", energy > 0) %>% pull(energy)
post_en <- data %>% filter(period == "Post-1990", energy > 0) %>% pull(energy)
wt2 <- wilcox.test(pre_en, post_en)
cat(sprintf("energy:  W = %.0f, p = %.4f\n", wt2$statistic, wt2$p.value))

pre_nl  <- data %>% filter(period == "Pre-1990", nlands > 0) %>% pull(nlands)
post_nl <- data %>% filter(period == "Post-1990", nlands > 0) %>% pull(nlands)
wt3 <- wilcox.test(pre_nl, post_nl)
cat(sprintf("nlands:  W = %.0f, p = %.4f\n", wt3$statistic, wt3$p.value))

# ==================================================================
# 3. Frequency: Fraction of country-years with any exposure
# ==================================================================
cat("\n\n================================================================\n")
cat("  3. SHOCK FREQUENCY\n")
cat("================================================================\n\n")

freq_stats <- data %>%
  group_by(period) %>%
  summarise(
    n_total = n(),
    # Any exposure
    n_any_maxwind = sum(maxwind > 0),
    n_any_energy  = sum(energy > 0),
    n_any_nlands  = sum(nlands > 0),
    frac_maxwind  = mean(maxwind > 0),
    frac_energy   = mean(energy > 0),
    frac_nlands   = mean(nlands > 0),
    # 95th percentile shocks
    n_maxwind_95  = sum(maxwind_95 == 1),
    n_energy_95   = sum(energy_95 == 1),
    n_nlands_95   = sum(nlands_95 == 1),
    frac_maxwind_95 = mean(maxwind_95 == 1),
    frac_energy_95  = mean(energy_95 == 1),
    frac_nlands_95  = mean(nlands_95 == 1),
    .groups = "drop"
  )

cat("--- Fraction of country-years with ANY exposure ---\n")
cat(sprintf("  maxwind > 0:  Pre=%.3f (%d/%d)  Post=%.3f (%d/%d)\n",
            freq_stats$frac_maxwind[freq_stats$period=="Pre-1990"],
            freq_stats$n_any_maxwind[freq_stats$period=="Pre-1990"],
            freq_stats$n_total[freq_stats$period=="Pre-1990"],
            freq_stats$frac_maxwind[freq_stats$period=="Post-1990"],
            freq_stats$n_any_maxwind[freq_stats$period=="Post-1990"],
            freq_stats$n_total[freq_stats$period=="Post-1990"]))
cat(sprintf("  energy > 0:   Pre=%.3f (%d/%d)  Post=%.3f (%d/%d)\n",
            freq_stats$frac_energy[freq_stats$period=="Pre-1990"],
            freq_stats$n_any_energy[freq_stats$period=="Pre-1990"],
            freq_stats$n_total[freq_stats$period=="Pre-1990"],
            freq_stats$frac_energy[freq_stats$period=="Post-1990"],
            freq_stats$n_any_energy[freq_stats$period=="Post-1990"],
            freq_stats$n_total[freq_stats$period=="Post-1990"]))
cat(sprintf("  nlands > 0:   Pre=%.3f (%d/%d)  Post=%.3f (%d/%d)\n",
            freq_stats$frac_nlands[freq_stats$period=="Pre-1990"],
            freq_stats$n_any_nlands[freq_stats$period=="Pre-1990"],
            freq_stats$n_total[freq_stats$period=="Pre-1990"],
            freq_stats$frac_nlands[freq_stats$period=="Post-1990"],
            freq_stats$n_any_nlands[freq_stats$period=="Post-1990"],
            freq_stats$n_total[freq_stats$period=="Post-1990"]))

cat("\n--- Fraction of country-years with 95th-percentile shock ---\n")
cat(sprintf("  maxwind_95=1: Pre=%.4f (%d/%d)  Post=%.4f (%d/%d)\n",
            freq_stats$frac_maxwind_95[freq_stats$period=="Pre-1990"],
            freq_stats$n_maxwind_95[freq_stats$period=="Pre-1990"],
            freq_stats$n_total[freq_stats$period=="Pre-1990"],
            freq_stats$frac_maxwind_95[freq_stats$period=="Post-1990"],
            freq_stats$n_maxwind_95[freq_stats$period=="Post-1990"],
            freq_stats$n_total[freq_stats$period=="Post-1990"]))
cat(sprintf("  energy_95=1:  Pre=%.4f (%d/%d)  Post=%.4f (%d/%d)\n",
            freq_stats$frac_energy_95[freq_stats$period=="Pre-1990"],
            freq_stats$n_energy_95[freq_stats$period=="Pre-1990"],
            freq_stats$n_total[freq_stats$period=="Pre-1990"],
            freq_stats$frac_energy_95[freq_stats$period=="Post-1990"],
            freq_stats$n_energy_95[freq_stats$period=="Post-1990"],
            freq_stats$n_total[freq_stats$period=="Post-1990"]))
cat(sprintf("  nlands_95=1:  Pre=%.4f (%d/%d)  Post=%.4f (%d/%d)\n",
            freq_stats$frac_nlands_95[freq_stats$period=="Pre-1990"],
            freq_stats$n_nlands_95[freq_stats$period=="Pre-1990"],
            freq_stats$n_total[freq_stats$period=="Pre-1990"],
            freq_stats$frac_nlands_95[freq_stats$period=="Post-1990"],
            freq_stats$n_nlands_95[freq_stats$period=="Post-1990"],
            freq_stats$n_total[freq_stats$period=="Post-1990"]))

# ==================================================================
# 4. Imbalanced treatment: what fraction of all 95th-pctile events
#    fall in each period?
# ==================================================================
cat("\n\n================================================================\n")
cat("  4. TREATMENT IMBALANCE: Share of p95 events by period\n")
cat("================================================================\n\n")

total_mw95 <- sum(data$maxwind_95)
total_en95 <- sum(data$energy_95)
total_nl95 <- sum(data$nlands_95)

pre_mw95  <- sum(data$maxwind_95[data$period == "Pre-1990"])
post_mw95 <- sum(data$maxwind_95[data$period == "Post-1990"])
pre_en95  <- sum(data$energy_95[data$period == "Pre-1990"])
post_en95 <- sum(data$energy_95[data$period == "Post-1990"])
pre_nl95  <- sum(data$nlands_95[data$period == "Pre-1990"])
post_nl95 <- sum(data$nlands_95[data$period == "Post-1990"])

cat(sprintf("maxwind_95: Total=%d, Pre=%d (%.1f%%), Post=%d (%.1f%%)\n",
            total_mw95, pre_mw95, 100*pre_mw95/total_mw95, post_mw95, 100*post_mw95/total_mw95))
cat(sprintf("energy_95:  Total=%d, Pre=%d (%.1f%%), Post=%d (%.1f%%)\n",
            total_en95, pre_en95, 100*pre_en95/total_en95, post_en95, 100*post_en95/total_en95))
cat(sprintf("nlands_95:  Total=%d, Pre=%d (%.1f%%), Post=%d (%.1f%%)\n",
            total_nl95, pre_nl95, 100*pre_nl95/total_nl95, post_nl95, 100*post_nl95/total_nl95))

# Also report the expected split if shocks were uniform
n_pre  <- sum(data$period == "Pre-1990")
n_post <- sum(data$period == "Post-1990")
cat(sprintf("\nExpected share if uniform: Pre=%.1f%%, Post=%.1f%%\n",
            100*n_pre/(n_pre+n_post), 100*n_post/(n_pre+n_post)))

# ==================================================================
# 5. Serial correlation in shocks: consecutive hits
# ==================================================================
cat("\n\n================================================================\n")
cat("  5. SERIAL CORRELATION: Consecutive shock years\n")
cat("================================================================\n\n")

# For each country, flag years where maxwind > 0 in both year t and year t-1
serial_data <- data %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  mutate(
    lag_maxwind    = lag(maxwind),
    lag_maxwind_95 = lag(maxwind_95),
    lag_energy     = lag(energy),
    lag_nlands     = lag(nlands),
    # Consecutive exposure: hit in year t AND year t-1
    consec_exposed  = (maxwind > 0) & (lag_maxwind > 0),
    consec_p95      = (maxwind_95 == 1) & (lag_maxwind_95 == 1),
    consec_nlands   = (nlands > 0) & (lag_nlands > 0),
    # Current year exposed
    exposed = maxwind > 0
  ) %>%
  ungroup()

consec_stats <- serial_data %>%
  filter(!is.na(lag_maxwind)) %>%
  group_by(period) %>%
  summarise(
    n = n(),
    # Among exposed country-years, what fraction were also exposed the previous year?
    n_exposed = sum(exposed),
    n_consec  = sum(consec_exposed, na.rm = TRUE),
    frac_consec_given_exposed = sum(consec_exposed, na.rm = TRUE) / sum(exposed),
    # Among p95 shocks, fraction preceded by any exposure
    n_p95 = sum(maxwind_95 == 1),
    n_p95_preceded = sum(consec_exposed & maxwind_95 == 1, na.rm = TRUE),
    frac_p95_preceded = ifelse(sum(maxwind_95 == 1) > 0,
                               sum(consec_exposed & maxwind_95 == 1, na.rm = TRUE) / sum(maxwind_95 == 1),
                               NA),
    # Consecutive extreme shocks
    n_consec_p95 = sum(consec_p95, na.rm = TRUE),
    # Consecutive landfalls
    n_consec_nlands = sum(consec_nlands, na.rm = TRUE),
    frac_consec_nlands = sum(consec_nlands, na.rm = TRUE) / sum(nlands > 0),
    .groups = "drop"
  )

cat("--- Consecutive exposure (maxwind > 0 in year t AND t-1) ---\n")
for (p in c("Pre-1990", "Post-1990")) {
  s <- consec_stats %>% filter(period == p)
  cat(sprintf("  %s: %d/%d exposed c-y had prior-year exposure (%.1f%%)\n",
              p, s$n_consec, s$n_exposed, 100*s$frac_consec_given_exposed))
}

cat("\n--- Consecutive landfall exposure (nlands > 0 in year t AND t-1) ---\n")
for (p in c("Pre-1990", "Post-1990")) {
  s <- consec_stats %>% filter(period == p)
  cat(sprintf("  %s: frac of landfall years preceded by landfall = %.1f%%\n",
              p, 100*s$frac_consec_nlands))
}

cat("\n--- P95 shocks preceded by any exposure ---\n")
for (p in c("Pre-1990", "Post-1990")) {
  s <- consec_stats %>% filter(period == p)
  cat(sprintf("  %s: %d/%d p95 shocks preceded by exposure (%.1f%%)\n",
              p, s$n_p95_preceded, s$n_p95, 100*s$frac_p95_preceded))
}

cat("\n--- Consecutive p95 shocks (maxwind_95=1 in year t AND t-1) ---\n")
for (p in c("Pre-1990", "Post-1990")) {
  s <- consec_stats %>% filter(period == p)
  cat(sprintf("  %s: %d events\n", p, s$n_consec_p95))
}

# ==================================================================
# 5b. Average gap between shocks (for exposed countries)
# ==================================================================
cat("\n--- Mean years between shocks (for ever-exposed countries) ---\n")

gap_data <- data %>%
  filter(maxwind > 0) %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  mutate(gap = year - lag(year)) %>%
  filter(!is.na(gap)) %>%
  ungroup() %>%
  mutate(period = ifelse(year <= 1990, "Pre-1990", "Post-1990"))

gap_stats <- gap_data %>%
  group_by(period) %>%
  summarise(
    mean_gap   = mean(gap),
    median_gap = median(gap),
    n_gaps     = n(),
    .groups = "drop"
  )

for (p in c("Pre-1990", "Post-1990")) {
  s <- gap_stats %>% filter(period == p)
  cat(sprintf("  %s: mean gap = %.2f years, median = %.1f years (n=%d inter-shock intervals)\n",
              p, s$mean_gap, s$median_gap, s$n_gaps))
}

# ==================================================================
# 5c. Autocorrelation of shock series within countries
# ==================================================================
cat("\n--- Within-country autocorrelation of maxwind ---\n")

autocorr <- data %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  filter(n() >= 10) %>%
  summarise(
    ar1 = cor(maxwind[-1], maxwind[-n()], use = "pairwise.complete.obs"),
    n_years = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(ar1))

cat(sprintf("  Mean AR(1) across %d countries: %.4f\n", nrow(autocorr), mean(autocorr$ar1)))
cat(sprintf("  Median AR(1): %.4f\n", median(autocorr$ar1)))

# Now split by period
autocorr_period <- data %>%
  arrange(countrycode, year) %>%
  group_by(countrycode, period) %>%
  filter(n() >= 5) %>%
  summarise(
    ar1 = cor(maxwind[-1], maxwind[-n()], use = "pairwise.complete.obs"),
    n_years = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(ar1))

for (p in c("Pre-1990", "Post-1990")) {
  ac <- autocorr_period %>% filter(period == p)
  cat(sprintf("  %s: mean AR(1) = %.4f across %d countries\n",
              p, mean(ac$ar1), nrow(ac)))
}

# ==================================================================
# 6. Country composition: which countries are being treated?
# ==================================================================
cat("\n\n================================================================\n")
cat("  6. COUNTRY COMPOSITION OF TREATED OBSERVATIONS\n")
cat("================================================================\n\n")

country_composition <- data %>%
  filter(maxwind_95 == 1) %>%
  group_by(period, countrycode) %>%
  summarise(n_shocks = n(), .groups = "drop")

cat("--- Top 10 countries by p95 shock count (Pre-1990) ---\n")
pre_top <- country_composition %>% filter(period == "Pre-1990") %>% arrange(desc(n_shocks)) %>% head(10)
for (i in 1:nrow(pre_top)) {
  cat(sprintf("  %s: %d\n", pre_top$countrycode[i], pre_top$n_shocks[i]))
}

cat("\n--- Top 10 countries by p95 shock count (Post-1990) ---\n")
post_top <- country_composition %>% filter(period == "Post-1990") %>% arrange(desc(n_shocks)) %>% head(10)
for (i in 1:nrow(post_top)) {
  cat(sprintf("  %s: %d\n", post_top$countrycode[i], post_top$n_shocks[i]))
}

# Countries unique to one period
pre_countries  <- country_composition %>% filter(period == "Pre-1990") %>% pull(countrycode)
post_countries <- country_composition %>% filter(period == "Post-1990") %>% pull(countrycode)
cat(sprintf("\nCountries with p95 shock ONLY pre-1990:  %d (%s)\n",
            length(setdiff(pre_countries, post_countries)),
            paste(setdiff(pre_countries, post_countries), collapse=", ")))
cat(sprintf("Countries with p95 shock ONLY post-1990: %d (%s)\n",
            length(setdiff(post_countries, pre_countries)),
            paste(setdiff(post_countries, pre_countries), collapse=", ")))
cat(sprintf("Countries with p95 shock in BOTH periods: %d\n",
            length(intersect(pre_countries, post_countries))))

# ==================================================================
# 7. Summary table
# ==================================================================
cat("\n\n================================================================\n")
cat("  7. SUMMARY TABLE\n")
cat("================================================================\n\n")

summary_table <- tibble(
  Metric = c(
    "Country-years",
    "Frac with any cyclone exposure (maxwind>0)",
    "Frac with p95 maxwind shock",
    "N p95 maxwind events",
    "Mean maxwind (cond. on >0, m/s)",
    "Median maxwind (cond. on >0, m/s)",
    "Mean energy (cond. on >0)",
    "Mean nlands (cond. on >0)",
    "Frac of all p95 events in this period",
    "Consec. exposure rate (given exposed)",
    "P95 shocks preceded by prior-yr exposure",
    "Mean inter-shock gap (years)",
    "Mean within-country AR(1) of maxwind"
  ),
  `Pre-1990` = c(
    sprintf("%d", freq_stats$n_total[freq_stats$period=="Pre-1990"]),
    sprintf("%.3f", freq_stats$frac_maxwind[freq_stats$period=="Pre-1990"]),
    sprintf("%.4f", freq_stats$frac_maxwind_95[freq_stats$period=="Pre-1990"]),
    sprintf("%d", freq_stats$n_maxwind_95[freq_stats$period=="Pre-1990"]),
    sprintf("%.2f", intensity_stats$maxwind_mean[intensity_stats$period=="Pre-1990"]),
    sprintf("%.2f", intensity_stats$maxwind_median[intensity_stats$period=="Pre-1990"]),
    sprintf("%.2f", intensity_stats$energy_mean[intensity_stats$period=="Pre-1990"]),
    sprintf("%.2f", intensity_stats$nlands_mean[intensity_stats$period=="Pre-1990"]),
    sprintf("%.1f%%", 100*pre_mw95/total_mw95),
    sprintf("%.1f%%", 100*consec_stats$frac_consec_given_exposed[consec_stats$period=="Pre-1990"]),
    sprintf("%.1f%%", 100*consec_stats$frac_p95_preceded[consec_stats$period=="Pre-1990"]),
    sprintf("%.2f", gap_stats$mean_gap[gap_stats$period=="Pre-1990"]),
    sprintf("%.4f", mean(autocorr_period$ar1[autocorr_period$period=="Pre-1990"]))
  ),
  `Post-1990` = c(
    sprintf("%d", freq_stats$n_total[freq_stats$period=="Post-1990"]),
    sprintf("%.3f", freq_stats$frac_maxwind[freq_stats$period=="Post-1990"]),
    sprintf("%.4f", freq_stats$frac_maxwind_95[freq_stats$period=="Post-1990"]),
    sprintf("%d", freq_stats$n_maxwind_95[freq_stats$period=="Post-1990"]),
    sprintf("%.2f", intensity_stats$maxwind_mean[intensity_stats$period=="Post-1990"]),
    sprintf("%.2f", intensity_stats$maxwind_median[intensity_stats$period=="Post-1990"]),
    sprintf("%.2f", intensity_stats$energy_mean[intensity_stats$period=="Post-1990"]),
    sprintf("%.2f", intensity_stats$nlands_mean[intensity_stats$period=="Post-1990"]),
    sprintf("%.1f%%", 100*post_mw95/total_mw95),
    sprintf("%.1f%%", 100*consec_stats$frac_consec_given_exposed[consec_stats$period=="Post-1990"]),
    sprintf("%.1f%%", 100*consec_stats$frac_p95_preceded[consec_stats$period=="Post-1990"]),
    sprintf("%.2f", gap_stats$mean_gap[gap_stats$period=="Post-1990"]),
    sprintf("%.4f", mean(autocorr_period$ar1[autocorr_period$period=="Post-1990"]))
  )
)

# Print nicely
max_metric_width <- max(nchar(summary_table$Metric))
max_pre_width    <- max(nchar(summary_table$`Pre-1990`))
max_post_width   <- max(nchar(summary_table$`Post-1990`))

header <- sprintf("%-*s | %-*s | %-*s",
                  max_metric_width, "Metric",
                  max_pre_width, "Pre-1990",
                  max_post_width, "Post-1990")
cat(header, "\n")
cat(paste(rep("-", nchar(header)), collapse=""), "\n")
for (i in 1:nrow(summary_table)) {
  cat(sprintf("%-*s | %-*s | %-*s\n",
              max_metric_width, summary_table$Metric[i],
              max_pre_width, summary_table$`Pre-1990`[i],
              max_post_width, summary_table$`Post-1990`[i]))
}

cat("\n\n================================================================\n")
cat("  INTERPRETATION\n")
cat("================================================================\n\n")
cat("See console output above for detailed findings.\n")
cat("Key questions answered:\n")
cat("  1. Do shock intensity distributions differ across periods?\n")
cat("  2. Are shocks more frequent post-1990?\n")
cat("  3. Is the global 95th-pctile threshold balanced across periods?\n")
cat("  4. Are consecutive shocks more common post-1990?\n")
cat("\nDone.\n")
