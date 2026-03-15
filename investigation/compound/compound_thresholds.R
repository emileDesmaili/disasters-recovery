# compound_thresholds.R
# Explore how varying the treatment definition (shock threshold)
# affects the compound vs isolated shock analysis
#
# The p95 threshold (60 knots) is very low — barely tropical storm.
# This means large cyclone-prone countries like AUS, CHN, PHL are
# "treated" nearly every year, making almost all their shocks compound.
# Higher thresholds should produce more isolated shocks even in
# these countries, providing better identifying variation.

library(tidyverse)
library(fixest)
library(haven)
library(countrycode)

source("../../emileRegs.R")
setFixest_notes(FALSE)

# ------------------------------------------------------------------
# 1. Load and construct data
# ------------------------------------------------------------------
pwt  <- read_stata("../../raw_data/pwt_clean.dta")
tcs  <- read_stata("../../raw_data/ibtracs_clean.dta")
data <- pwt %>% left_join(tcs, by = c("year", "countrycode"))

knot_to_ms <- 0.514444

data <- data %>%
  mutate(
    maxwind  = replace_na(max_ann_wind_i_nob, 0) * knot_to_ms,
    energy   = replace_na(sum_ann_energy_i_nob, 0),
    loggdp   = ln_real_gdp_usd_pc,
    gdp_diff = growth_real_gdp_usd_pc,
    year2    = year^2,
    period   = ifelse(year <= 1990, "Pre-1990", "Post-1990"),
    region   = countrycode(countrycode, "iso3c", "region"),
    region   = ifelse(is.na(region), "Other", region)
  ) %>%
  filter(year >= 1969)

# ------------------------------------------------------------------
# 2. Define multiple thresholds
# ------------------------------------------------------------------
# Percentile-based (computed among ALL country-years including zeros)
p95  <- quantile(data$maxwind[data$year >= 1970 & data$year <= 2014], 0.95, na.rm = TRUE)
p975 <- quantile(data$maxwind[data$year >= 1970 & data$year <= 2014], 0.975, na.rm = TRUE)
p99  <- quantile(data$maxwind[data$year >= 1970 & data$year <= 2014], 0.99, na.rm = TRUE)

# Category-based (Saffir-Simpson scale, in m/s)
cat1_thresh <- 64 * knot_to_ms   # 33 m/s — Category 1 hurricane
cat2_thresh <- 83 * knot_to_ms   # 43 m/s — Category 2
cat3_thresh <- 96 * knot_to_ms   # 49 m/s — Category 3 (major)

thresholds <- list(
  "p95 (60 kts)"    = p95,
  "p97.5 (72 kts)"  = p975,
  "p99 (90 kts)"    = p99,
  "Cat 1+ (64 kts)" = cat1_thresh,
  "Cat 2+ (83 kts)" = cat2_thresh,
  "Cat 3+ (96 kts)" = cat3_thresh
)

cat("========== THRESHOLD VALUES ==========\n")
for (nm in names(thresholds)) {
  cat(sprintf("  %-20s: %.1f m/s = %.0f knots\n",
              nm, thresholds[[nm]], thresholds[[nm]] / knot_to_ms))
}

# ------------------------------------------------------------------
# 3. Characterize each threshold
# ------------------------------------------------------------------
cat("\n========== THRESHOLD DIAGNOSTICS ==========\n")
cat(sprintf("%-20s %8s %8s %8s %8s %8s %8s\n",
            "Threshold", "Shocks", "Treated", "Compound", "Isolated", "%Comp", "Both"))
cat(paste(rep("-", 85), collapse = ""), "\n")

threshold_stats <- list()

for (nm in names(thresholds)) {
  thresh <- thresholds[[nm]]

  d <- data %>%
    mutate(shock = as.integer(thresh > 0 & maxwind >= thresh))

  treated <- d %>% filter(shock == 1) %>% pull(countrycode) %>% unique()

  d_tr <- d %>%
    filter(countrycode %in% treated) %>%
    arrange(countrycode, year) %>%
    group_by(countrycode) %>%
    mutate(
      shock_lag1 = lag(shock, 1, default = 0L),
      compound   = as.integer(shock == 1 & shock_lag1 == 1),
      isolated   = as.integer(shock == 1 & shock_lag1 == 0)
    ) %>%
    ungroup()

  n_shocks   <- sum(d_tr$shock)
  n_compound <- sum(d_tr$compound)
  n_isolated <- sum(d_tr$isolated)

  # Countries with both types
  country_types <- d_tr %>%
    filter(shock == 1) %>%
    group_by(countrycode) %>%
    summarise(has_comp = any(compound == 1), has_isol = any(isolated == 1))
  n_both <- sum(country_types$has_comp & country_types$has_isol)

  cat(sprintf("%-20s %8d %8d %8d %8d %7.1f%% %8d\n",
              nm, n_shocks, length(treated), n_compound, n_isolated,
              n_compound / n_shocks * 100, n_both))

  threshold_stats[[nm]] <- list(
    name = nm, thresh = thresh, n_shocks = n_shocks,
    n_treated = length(treated), n_compound = n_compound,
    n_isolated = n_isolated, n_both = n_both,
    pct_compound = n_compound / n_shocks * 100
  )
}

# ------------------------------------------------------------------
# 4. Australia specifically across thresholds
# ------------------------------------------------------------------
cat("\n========== AUSTRALIA ACROSS THRESHOLDS ==========\n")
for (nm in names(thresholds)) {
  thresh <- thresholds[[nm]]
  aus <- data %>% filter(countrycode == "AUS")
  n_above <- sum(aus$maxwind >= thresh)
  cat(sprintf("  %-20s: %d / %d years above threshold\n", nm, n_above, nrow(aus)))
}

# ------------------------------------------------------------------
# 5. Estimate compound vs isolated IRFs at each threshold
# ------------------------------------------------------------------
cat("\n========== ESTIMATING IRFs ==========\n")

outcome  <- "loggdp"
horizon  <- 10
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

all_irfs <- list()

for (nm in names(thresholds)) {
  thresh <- thresholds[[nm]]
  cat(sprintf("Estimating: %s\n", nm))

  # Create shock and compound variables
  d <- data %>%
    mutate(shock = as.integer(thresh > 0 & maxwind >= thresh))

  treated <- d %>% filter(shock == 1) %>% pull(countrycode) %>% unique()

  d_tr <- d %>%
    filter(countrycode %in% treated) %>%
    arrange(countrycode, year) %>%
    group_by(countrycode) %>%
    mutate(
      shock_lag1     = lag(shock, 1, default = 0L),
      had_recent     = as.factor(shock_lag1)
    ) %>%
    ungroup()

  # Skip if too few shocks or no within-country variation
  n_comp <- sum(d_tr$shock == 1 & d_tr$shock_lag1 == 1)
  n_isol <- sum(d_tr$shock == 1 & d_tr$shock_lag1 == 0)

  if (n_comp < 5 || n_isol < 5) {
    cat(sprintf("  Skipping %s: too few compound (%d) or isolated (%d)\n",
                nm, n_comp, n_isol))
    next
  }

  controls <- "l(gdp_diff,1:2) + l(shock,1:2)"

  irf <- tryCatch({
    lp_panel_inter(
      data = d_tr, outcome = outcome, main_var = "shock",
      interact_var = "had_recent",
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
        ),
        threshold = nm
      )
  }, error = function(e) {
    cat(sprintf("  Error for %s: %s\n", nm, e$message))
    NULL
  })

  if (!is.null(irf)) {
    all_irfs[[nm]] <- irf
  }
}

irf_all <- bind_rows(all_irfs)

# ------------------------------------------------------------------
# 6. Print h=5 results across thresholds
# ------------------------------------------------------------------
cat("\n========== COMPOUND vs ISOLATED AT h=5 ACROSS THRESHOLDS ==========\n")
cat(sprintf("%-20s %10s %10s %10s\n", "Threshold", "Isolated", "Compound", "Diff (pp)"))
cat(paste(rep("-", 55), collapse = ""), "\n")

for (nm in names(thresholds)) {
  h5 <- irf_all %>% filter(threshold == nm, horizon == 5)
  if (nrow(h5) == 0) next
  isol <- h5 %>% filter(category == "Isolated") %>% pull(irf_mean)
  comp <- h5 %>% filter(category == "Compound") %>% pull(irf_mean)
  if (length(isol) > 0 && length(comp) > 0)
    cat(sprintf("%-20s %9.1f%% %9.1f%% %9.2f\n",
                nm, isol * 100, comp * 100, (comp - isol) * 100))
}

# ------------------------------------------------------------------
# 7. Figures
# ------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# 7a. IRF comparison across thresholds
irf_all <- irf_all %>%
  mutate(threshold = factor(threshold, levels = names(thresholds)))

p_thresh <- ggplot(irf_all, aes(x = horizon, y = irf_mean,
                                 color = category, fill = category)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.2) +
  facet_wrap(~threshold, ncol = 3) +
  scale_color_manual(values = c("Isolated" = "steelblue", "Compound" = "firebrick")) +
  scale_fill_manual(values = c("Isolated" = "steelblue", "Compound" = "firebrick")) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Compound vs Isolated Shocks Across Treatment Thresholds",
    subtitle = "shock × 1(shock in t-1) | Ever-treated sample, year FE + country trends",
    color = NULL, fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40", size = 9),
    legend.position = "bottom"
  )

ggsave("figures/irf_compound_thresholds.png", p_thresh,
       width = 13, height = 8, dpi = 300)
cat("\nSaved: figures/irf_compound_thresholds.png\n")

# 7b. Threshold diagnostics summary plot
stats_df <- bind_rows(lapply(threshold_stats, function(x) {
  data.frame(threshold = x$name, n_shocks = x$n_shocks,
             n_treated = x$n_treated, n_compound = x$n_compound,
             n_isolated = x$n_isolated, pct_compound = x$pct_compound,
             n_both = x$n_both)
})) %>%
  mutate(threshold = factor(threshold, levels = names(thresholds)))

p_stats <- stats_df %>%
  select(threshold, Compound = n_compound, Isolated = n_isolated) %>%
  pivot_longer(-threshold, names_to = "type", values_to = "count") %>%
  ggplot(aes(x = threshold, y = count, fill = type)) +
  geom_col(position = "stack", alpha = 0.8) +
  scale_fill_manual(values = c("Compound" = "firebrick", "Isolated" = "steelblue")) +
  labs(x = NULL, y = "Number of shock-years",
       title = "Shock Composition Across Thresholds",
       subtitle = "Higher thresholds → fewer shocks, lower compound share",
       fill = NULL) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave("figures/threshold_composition.png", p_stats,
       width = 8, height = 5, dpi = 300)
cat("Saved: figures/threshold_composition.png\n")

# 7c. h=5 gap across thresholds
gap_df <- irf_all %>%
  filter(horizon == 5) %>%
  select(threshold, category, irf_mean, se) %>%
  pivot_wider(names_from = category, values_from = c(irf_mean, se)) %>%
  mutate(
    gap = irf_mean_Compound - irf_mean_Isolated,
    gap_se = sqrt(se_Compound^2 + se_Isolated^2)  # conservative (ignores covariance)
  )

p_gap <- ggplot(gap_df, aes(x = threshold, y = gap)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_pointrange(aes(ymin = gap - 1.96 * gap_se, ymax = gap + 1.96 * gap_se),
                  size = 0.6, color = "grey30") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Compound - Isolated effect at h=5",
       title = "Compound-Isolated Gap Across Thresholds",
       subtitle = "Positive = compound less damaging | 95% CI (conservative)") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave("figures/threshold_gap.png", p_gap, width = 8, height = 5, dpi = 300)
cat("Saved: figures/threshold_gap.png\n")

# ------------------------------------------------------------------
# 8. Country-level detail at selected thresholds
# ------------------------------------------------------------------
cat("\n========== TOP COUNTRIES BY THRESHOLD ==========\n")
for (nm in c("p95 (60 kts)", "Cat 1+ (64 kts)", "Cat 2+ (83 kts)", "Cat 3+ (96 kts)")) {
  thresh <- thresholds[[nm]]
  d <- data %>%
    mutate(shock = as.integer(thresh > 0 & maxwind >= thresh))
  country_tab <- d %>%
    filter(shock == 1) %>%
    group_by(countrycode) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    head(10)
  cat(sprintf("\n--- %s ---\n", nm))
  cat(sprintf("  %-5s  %s\n", "ISO", "Shocks"))
  for (i in 1:nrow(country_tab))
    cat(sprintf("  %-5s  %d\n", country_tab$countrycode[i], country_tab$n[i]))
}

saveRDS(list(irf_all = irf_all, threshold_stats = threshold_stats,
             gap_df = gap_df),
        "threshold_results.rds")
cat("\nSaved: threshold_results.rds\n")
cat("Done.\n")
