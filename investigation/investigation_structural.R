# investigation_structural.R
# Test whether the post-1990 divergence in cyclone-GDP IRFs is explained
# by structural economic changes, using variables available in PWT data.
#
# Three analyses:
#   1. Capital intensity interaction (rkna / real_gdp_usd)
#   2. Continuous time trend interaction (year - 1990)
#   3. Human capital interaction (hc index)

library(tidyverse)
library(fixest)
library(haven)
library(patchwork)

source("../emileRegs.R")

# ------------------------------------------------------------------
# 1. Load and construct data (same pipeline as other investigation scripts)
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

# Global 95th-percentile thresholds
data <- data %>%
  mutate(
    maxwind_p95 = quantile(maxwind[year >= 1970 & year <= 2014], 0.95, na.rm = TRUE)
  ) %>%
  mutate(
    maxwind_95 = as.integer(maxwind_p95 > 0 & maxwind >= maxwind_p95)
  )

# Standardize continuous measures
data <- data %>%
  mutate(
    maxwind_sd = maxwind / sd(maxwind[maxwind > 0], na.rm = TRUE)
  )

# ------------------------------------------------------------------
# 2. LP settings (matching lp_gdp_shocks.R)
# ------------------------------------------------------------------
outcome  <- "loggdp"
horizon  <- 10
fe       <- "countrycode[year] + countrycode[year2] + year"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

make_controls <- function(shock_var) {
  paste0("l(gdp_diff,1:2) + l(", shock_var, ",1:2)")
}

# ------------------------------------------------------------------
# 3. State-dependent LP function (from lp_gdp_shocks.R lines 213-254)
#    Adapted to accept custom evaluation points and labels
# ------------------------------------------------------------------
lp_state_dep <- function(data, outcome, main_var, state_var,
                         controls = NULL, horizon = 10,
                         fe = "countrycode + year",
                         panel_id = c("countrycode", "year"),
                         vcov_formula = DK ~ year,
                         eval_quantiles = c(0.10, 0.50, 0.90),
                         eval_labels = NULL,
                         eval_values = NULL) {
  rhs_controls <- if (!is.null(controls)) paste0(" + ", controls) else ""
  irf_list <- list()

  # Determine evaluation points
  if (is.null(eval_values)) {
    state_vals <- quantile(data[[state_var]], eval_quantiles, na.rm = TRUE)
  } else {
    state_vals <- eval_values
  }

  if (is.null(eval_labels)) {
    eval_labels <- paste0("p", eval_quantiles * 100)
  }

  for (h in 0:horizon) {
    fml <- as.formula(paste0(
      "f(", outcome, ", ", h, ") - l(", outcome, ", 1) ~ ",
      main_var, " + ", main_var, ":", state_var,
      rhs_controls, " | ", fe
    ))
    mod <- feols(fml, data = data, panel.id = panel_id, vcov = vcov_formula)

    # Extract coefficients
    beta_main  <- coef(mod)[main_var]
    inter_name <- grep(paste0(main_var, ":", state_var), names(coef(mod)), value = TRUE)
    if (length(inter_name) == 0)
      inter_name <- grep(paste0(state_var, ":", main_var), names(coef(mod)), value = TRUE)
    beta_inter <- coef(mod)[inter_name]

    # Evaluate marginal effect at each state value
    for (j in seq_along(state_vals)) {
      me <- beta_main + beta_inter * state_vals[j]
      # Delta method SE
      grad <- c(1, state_vals[j])
      se <- sqrt(t(grad) %*% vcov(mod)[c(main_var, inter_name),
                                        c(main_var, inter_name)] %*% grad)
      irf_list[[length(irf_list) + 1]] <- data.frame(
        horizon   = h,
        irf_mean  = me,
        se        = as.numeric(se),
        irf_down  = me - 1.96 * as.numeric(se),
        irf_up    = me + 1.96 * as.numeric(se),
        quantile  = eval_labels[j],
        state_val = state_vals[j]
      )
    }
  }
  bind_rows(irf_list)
}

# ------------------------------------------------------------------
# Plotting theme (consistent with other investigation scripts)
# ------------------------------------------------------------------
irf_theme <- theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40")
  )


# ==================================================================
# ANALYSIS 1: Capital Intensity Interaction
# ==================================================================
cat("\n============================================================\n")
cat("  ANALYSIS 1: Capital Intensity (rkna / real_gdp_usd)\n")
cat("============================================================\n\n")

# Construct capital/GDP ratio and demean
data <- data %>%
  mutate(
    cap_gdp = rkna / real_gdp_usd,
    cap_gdp_dm = cap_gdp - mean(cap_gdp, na.rm = TRUE)
  )

cat(sprintf("Capital/GDP ratio summary:\n"))
cat(sprintf("  Mean:   %.3f\n", mean(data$cap_gdp, na.rm = TRUE)))
cat(sprintf("  Median: %.3f\n", median(data$cap_gdp, na.rm = TRUE)))
cat(sprintf("  SD:     %.3f\n", sd(data$cap_gdp, na.rm = TRUE)))
cat(sprintf("  p10:    %.3f\n", quantile(data$cap_gdp, 0.10, na.rm = TRUE)))
cat(sprintf("  p50:    %.3f\n", quantile(data$cap_gdp, 0.50, na.rm = TRUE)))
cat(sprintf("  p90:    %.3f\n", quantile(data$cap_gdp, 0.90, na.rm = TRUE)))
cat(sprintf("  N non-missing: %d\n", sum(!is.na(data$cap_gdp))))

# Run state-dependent LP: maxwind_sd interacted with cap_gdp_dm
irfs_capital <- lp_state_dep(
  data = data,
  outcome = outcome,
  main_var = "maxwind_sd",
  state_var = "cap_gdp_dm",
  controls = make_controls("maxwind_sd"),
  horizon = horizon,
  fe = fe,
  panel_id = panel_id,
  vcov_formula = vcov_fm
)

irfs_capital <- irfs_capital %>%
  mutate(quantile = factor(quantile, levels = c("p10", "p50", "p90"),
                           labels = c("Low capital intensity (p10)",
                                      "Median (p50)",
                                      "High capital intensity (p90)")))

# Print summary table
cat("\n--- IRF at h=5 and h=10 by capital intensity ---\n")
cap_summary <- irfs_capital %>%
  filter(horizon %in% c(5, 10)) %>%
  select(horizon, quantile, irf_mean, se, irf_down, irf_up, state_val) %>%
  mutate(across(c(irf_mean, se, irf_down, irf_up), ~ round(. * 100, 3)))
print(as.data.frame(cap_summary), row.names = FALSE)

# Plot
p_capital <- ggplot(irfs_capital, aes(x = horizon, y = irf_mean,
                                      color = quantile, fill = quantile)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Low capital intensity (p10)" = "steelblue",
                                "Median (p50)" = "grey40",
                                "High capital intensity (p90)" = "firebrick")) +
  scale_fill_manual(values = c("Low capital intensity (p10)" = "steelblue",
                               "Median (p50)" = "grey40",
                               "High capital intensity (p90)" = "firebrick")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "Horizon (years)",
    y = "GDP response (per 1 SD wind)",
    color = "Capital/GDP ratio",
    fill = "Capital/GDP ratio",
    title = "Cyclone IRF by Capital Intensity",
    subtitle = "Marginal effect of maxwind (1 SD) at p10/p50/p90 of capital/GDP"
  ) +
  irf_theme +
  theme(legend.position = "top")

ggsave("figures/irf_capital_intensity.png", p_capital,
       width = 8, height = 6, dpi = 300)
cat("\nFigure saved to figures/irf_capital_intensity.png\n")


# ==================================================================
# ANALYSIS 2: Continuous Time Trend Interaction
# ==================================================================
cat("\n\n============================================================\n")
cat("  ANALYSIS 2: Continuous Time Trend (year - 1990)\n")
cat("============================================================\n\n")

# Construct time trend centered at 1990
data <- data %>%
  mutate(time_trend = year - 1990)

cat(sprintf("Time trend range: %d to %d (0 = 1990)\n",
            min(data$time_trend), max(data$time_trend)))

# Evaluation points: -15 (approx 1975), 0 (1990), +15 (approx 2005)
eval_vals <- c(-15, 0, 15)
eval_labs <- c("~1975 (t=-15)", "1990 (t=0)", "~2005 (t=+15)")

irfs_trend <- lp_state_dep(
  data = data,
  outcome = outcome,
  main_var = "maxwind_95",
  state_var = "time_trend",
  controls = make_controls("maxwind_95"),
  horizon = horizon,
  fe = fe,
  panel_id = panel_id,
  vcov_formula = vcov_fm,
  eval_values = eval_vals,
  eval_labels = eval_labs
)

irfs_trend <- irfs_trend %>%
  mutate(quantile = factor(quantile,
                           levels = c("~1975 (t=-15)", "1990 (t=0)", "~2005 (t=+15)")))

# Print summary table
cat("\n--- IRF at h=5 and h=10 by time trend ---\n")
trend_summary <- irfs_trend %>%
  filter(horizon %in% c(5, 10)) %>%
  select(horizon, quantile, irf_mean, se, irf_down, irf_up, state_val) %>%
  mutate(across(c(irf_mean, se, irf_down, irf_up), ~ round(. * 100, 3)))
print(as.data.frame(trend_summary), row.names = FALSE)

# Plot
p_trend <- ggplot(irfs_trend, aes(x = horizon, y = irf_mean,
                                   color = quantile, fill = quantile)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("~1975 (t=-15)" = "steelblue",
                                "1990 (t=0)" = "grey40",
                                "~2005 (t=+15)" = "firebrick")) +
  scale_fill_manual(values = c("~1975 (t=-15)" = "steelblue",
                               "1990 (t=0)" = "grey40",
                               "~2005 (t=+15)" = "firebrick")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "Horizon (years)",
    y = "GDP response to maxwind p95 shock",
    color = "Evaluated at",
    fill = "Evaluated at",
    title = "Cyclone IRF with Continuous Time Trend Interaction",
    subtitle = "maxwind_95 x (year - 1990); marginal effect at 1975, 1990, 2005"
  ) +
  irf_theme +
  theme(legend.position = "top")

ggsave("figures/irf_time_trend.png", p_trend,
       width = 8, height = 6, dpi = 300)
cat("\nFigure saved to figures/irf_time_trend.png\n")


# ==================================================================
# ANALYSIS 3: Human Capital Interaction
# ==================================================================
cat("\n\n============================================================\n")
cat("  ANALYSIS 3: Human Capital (hc index)\n")
cat("============================================================\n\n")

# Demean human capital
data <- data %>%
  mutate(hc_dm = hc - mean(hc, na.rm = TRUE))

cat(sprintf("Human capital index (hc) summary:\n"))
cat(sprintf("  Mean:   %.3f\n", mean(data$hc, na.rm = TRUE)))
cat(sprintf("  Median: %.3f\n", median(data$hc, na.rm = TRUE)))
cat(sprintf("  SD:     %.3f\n", sd(data$hc, na.rm = TRUE)))
cat(sprintf("  p10:    %.3f\n", quantile(data$hc, 0.10, na.rm = TRUE)))
cat(sprintf("  p50:    %.3f\n", quantile(data$hc, 0.50, na.rm = TRUE)))
cat(sprintf("  p90:    %.3f\n", quantile(data$hc, 0.90, na.rm = TRUE)))
cat(sprintf("  N non-missing: %d\n", sum(!is.na(data$hc))))

# Run state-dependent LP: maxwind_95 interacted with hc_dm
irfs_hc <- lp_state_dep(
  data = data,
  outcome = outcome,
  main_var = "maxwind_95",
  state_var = "hc_dm",
  controls = make_controls("maxwind_95"),
  horizon = horizon,
  fe = fe,
  panel_id = panel_id,
  vcov_formula = vcov_fm
)

irfs_hc <- irfs_hc %>%
  mutate(quantile = factor(quantile, levels = c("p10", "p50", "p90"),
                           labels = c("Low human capital (p10)",
                                      "Median (p50)",
                                      "High human capital (p90)")))

# Print summary table
cat("\n--- IRF at h=5 and h=10 by human capital ---\n")
hc_summary <- irfs_hc %>%
  filter(horizon %in% c(5, 10)) %>%
  select(horizon, quantile, irf_mean, se, irf_down, irf_up, state_val) %>%
  mutate(across(c(irf_mean, se, irf_down, irf_up), ~ round(. * 100, 3)))
print(as.data.frame(hc_summary), row.names = FALSE)

# Plot
p_hc <- ggplot(irfs_hc, aes(x = horizon, y = irf_mean,
                              color = quantile, fill = quantile)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Low human capital (p10)" = "steelblue",
                                "Median (p50)" = "grey40",
                                "High human capital (p90)" = "firebrick")) +
  scale_fill_manual(values = c("Low human capital (p10)" = "steelblue",
                               "Median (p50)" = "grey40",
                               "High human capital (p90)" = "firebrick")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "Horizon (years)",
    y = "GDP response to maxwind p95 shock",
    color = "Human capital index",
    fill = "Human capital index",
    title = "Cyclone IRF by Human Capital Level",
    subtitle = "Marginal effect of maxwind p95 shock at p10/p50/p90 of hc index"
  ) +
  irf_theme +
  theme(legend.position = "top")

ggsave("figures/irf_human_capital.png", p_hc,
       width = 8, height = 6, dpi = 300)
cat("\nFigure saved to figures/irf_human_capital.png\n")


# ==================================================================
# SUMMARY
# ==================================================================
cat("\n\n============================================================\n")
cat("  SUMMARY: Structural Change Explanations\n")
cat("============================================================\n\n")

# Compare the spread of IRFs across the three dimensions
for (analysis_name in c("Capital Intensity", "Time Trend", "Human Capital")) {
  df <- switch(analysis_name,
               "Capital Intensity" = irfs_capital,
               "Time Trend"        = irfs_trend,
               "Human Capital"     = irfs_hc)

  h10 <- df %>% filter(horizon == 10)
  spread <- max(h10$irf_mean) - min(h10$irf_mean)
  cat(sprintf("%s:\n", analysis_name))
  for (i in 1:nrow(h10)) {
    cat(sprintf("  %s: h=10 IRF = %.4f (%.2f%%)\n",
                h10$quantile[i], h10$irf_mean[i], h10$irf_mean[i] * 100))
  }
  cat(sprintf("  Spread (max - min at h=10): %.4f (%.2f pp)\n\n",
              spread, spread * 100))
}

cat("Figures saved:\n")
cat("  - figures/irf_capital_intensity.png\n")
cat("  - figures/irf_time_trend.png\n")
cat("  - figures/irf_human_capital.png\n")
cat("\nDone.\n")
