# investigation_placebo_scramble.R
# Shock scrambling placebo test for the period-interacted LP
#
# Logic: If the post-pre gap in the period-interacted LP is driven by actual
# cyclone damage (real shock-GDP relationship), then randomly permuting the
# shock indicator among treated country-years should eliminate the gap
# (distribution centered on zero). If the gap persists under scrambling, it
# means the gap is an artifact of the estimation framework (year FE structure,
# compositional differences between treated and never-treated countries).
#
# Design:
#   1. Keep full sample (all 182 countries)
#   2. Among the 40 treated countries, randomly permute maxwind_95 across
#      their country-years (preserving which countries are treated vs
#      never-treated, and total number of shocks)
#   3. Run period-interacted LP at h=5 on scrambled data
#   4. Record post-pre gap
#   5. Repeat 500 times
#   6. Plot distribution of scrambled gaps vs actual gap

library(tidyverse)
library(fixest)
library(haven)

source("emileRegs.R")

setFixest_notes(FALSE)
set.seed(42)

cat("============================================================\n")
cat("  PLACEBO TEST: SHOCK SCRAMBLING\n")
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

# ------------------------------------------------------------------
# 2. LP settings
# ------------------------------------------------------------------
outcome  <- "loggdp"
horizon  <- 5
controls <- "l(gdp_diff,1:2) + l(maxwind_95,1:2)"
fe       <- "countrycode[year] + countrycode[year2] + year"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

# ------------------------------------------------------------------
# 3. Identify treated countries
# ------------------------------------------------------------------
ever_treated <- data %>%
  filter(maxwind_95 == 1) %>%
  pull(countrycode) %>%
  unique()

treated_idx <- which(data$countrycode %in% ever_treated)

cat(sprintf("Total countries:       %d\n", n_distinct(data$countrycode)))
cat(sprintf("Ever-treated:          %d\n", length(ever_treated)))
cat(sprintf("Never-treated:         %d\n",
            n_distinct(data$countrycode) - length(ever_treated)))
cat(sprintf("Treated country-years: %d\n", length(treated_idx)))
cat(sprintf("Total shocks:          %d\n", sum(data$maxwind_95)))
cat("\n")

# ------------------------------------------------------------------
# 4. Helper: extract post-pre gap at h=5 from lp_panel_inter output
# ------------------------------------------------------------------
extract_gap <- function(irf_df) {
  # Clean category names
  irf_df <- irf_df %>%
    mutate(
      category = str_replace_all(category, "_", "-"),
      category = case_when(
        str_detect(category, "Pre")  ~ "Pre-1990",
        str_detect(category, "Post") ~ "Post-1990",
        TRUE ~ category
      )
    )

  h5 <- irf_df %>% filter(horizon == 5)
  post <- h5 %>% filter(category == "Post-1990") %>% pull(irf_mean)
  pre  <- h5 %>% filter(category == "Pre-1990")  %>% pull(irf_mean)

  if (length(post) == 1 & length(pre) == 1) {
    return(post - pre)
  } else {
    return(NA_real_)
  }
}

# ------------------------------------------------------------------
# 5. Compute actual gap
# ------------------------------------------------------------------
cat("Computing actual gap from real data...\n")
irf_actual <- lp_panel_inter(
  data = data, outcome = outcome, main_var = "maxwind_95",
  interact_var = "period", controls = controls,
  horizon = horizon, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
)

actual_gap <- extract_gap(irf_actual)
cat(sprintf("Actual post-pre gap at h=5: %.4f (%.2f pp)\n\n", actual_gap, actual_gap * 100))

# ------------------------------------------------------------------
# 6. Scrambling simulation
# ------------------------------------------------------------------
n_sims <- 500
scrambled_gaps <- numeric(n_sims)

cat(sprintf("Running %d scrambling simulations...\n", n_sims))
t_start <- Sys.time()

for (s in 1:n_sims) {
  if (s %% 50 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    cat(sprintf("  Sim %d/%d  (%.1f min elapsed)\n", s, n_sims, elapsed))
  }

  # Copy data
  d_sim <- data

  # Permute maxwind_95 among treated country-years
  d_sim$maxwind_95[treated_idx] <- sample(d_sim$maxwind_95[treated_idx])

  # Run LP (lags of maxwind_95 computed internally by fixest from the panel)
  result <- tryCatch({
    irf_sim <- lp_panel_inter(
      data = d_sim, outcome = outcome, main_var = "maxwind_95",
      interact_var = "period", controls = controls,
      horizon = horizon, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
    )
    extract_gap(irf_sim)
  }, error = function(e) {
    NA_real_
  })

  scrambled_gaps[s] <- result
}

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("\nDone. Total time: %.1f minutes\n", elapsed_total))

# ------------------------------------------------------------------
# 7. Summary statistics
# ------------------------------------------------------------------
valid_gaps <- scrambled_gaps[!is.na(scrambled_gaps)]
n_valid <- length(valid_gaps)
n_failed <- n_sims - n_valid

# Two-sided p-value: fraction of scrambled gaps at least as extreme as actual
p_value <- mean(abs(valid_gaps) >= abs(actual_gap))

cat("\n============================================================\n")
cat("  RESULTS\n")
cat("============================================================\n\n")
cat(sprintf("Simulations completed: %d / %d  (%d failed)\n", n_valid, n_sims, n_failed))
cat(sprintf("Actual gap (post-pre at h=5): %.4f (%.2f pp)\n", actual_gap, actual_gap * 100))
cat(sprintf("\nScrambled gap distribution:\n"))
cat(sprintf("  Mean:   %.4f (%.2f pp)\n", mean(valid_gaps), mean(valid_gaps) * 100))
cat(sprintf("  Median: %.4f (%.2f pp)\n", median(valid_gaps), median(valid_gaps) * 100))
cat(sprintf("  SD:     %.4f (%.2f pp)\n", sd(valid_gaps), sd(valid_gaps) * 100))
cat(sprintf("  Min:    %.4f (%.2f pp)\n", min(valid_gaps), min(valid_gaps) * 100))
cat(sprintf("  Max:    %.4f (%.2f pp)\n", max(valid_gaps), max(valid_gaps) * 100))
cat(sprintf("\nTwo-sided p-value: %.4f\n", p_value))

if (p_value < 0.05) {
  cat("=> The actual gap is statistically unusual under shock scrambling.\n")
  cat("   The gap appears driven by the ACTUAL shock timing, not by the\n")
  cat("   estimation framework or compositional differences.\n")
} else {
  cat("=> The actual gap is NOT unusual under shock scrambling.\n")
  cat("   The gap may be an artifact of the estimation framework\n")
  cat("   (year FE dominated by never-treated countries, differential trends).\n")
}

# ------------------------------------------------------------------
# 8. Plot
# ------------------------------------------------------------------
dir.create("investigation/figures", showWarnings = FALSE, recursive = TRUE)

plot_df <- data.frame(gap = valid_gaps)

p <- ggplot(plot_df, aes(x = gap)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "grey70", color = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = actual_gap, color = "firebrick",
             linewidth = 1.2, linetype = "solid") +
  annotate("text",
           x = actual_gap, y = Inf, vjust = 2, hjust = -0.1,
           label = sprintf("Actual gap = %.2f pp\np-value = %.3f",
                           actual_gap * 100, p_value),
           color = "firebrick", fontface = "bold", size = 4.5) +
  labs(
    x = "Post-1990 minus Pre-1990 gap at h = 5",
    y = "Density",
    title = "Placebo Test: Shock Scrambling Among Treated Countries",
    subtitle = sprintf(
      "%d simulations | Shocks permuted within %d treated countries | %d never-treated unchanged",
      n_valid, length(ever_treated),
      n_distinct(data$countrycode) - length(ever_treated)
    )
  ) +
  scale_x_continuous(labels = function(x) paste0(round(x * 100, 1), " pp")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40")
  )

ggsave("investigation/figures/placebo_shock_scramble.png", p,
       width = 9, height = 6, dpi = 300)
cat("\nFigure saved: investigation/figures/placebo_shock_scramble.png\n")
cat("\nDone.\n")
