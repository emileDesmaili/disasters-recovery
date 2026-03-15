# matching_summary.R
# Run all matching methods and produce comparison figures
# This is the main entry point: run this script to produce all results

library(tidyverse)
library(fixest)
library(haven)
library(MatchIt)
library(countrycode)

source("../../emileRegs.R")
setFixest_notes(FALSE)

# ------------------------------------------------------------------
# 1. Run all matching scripts
# ------------------------------------------------------------------
cat("=== Step 1: Diagnostics ===\n")
source("matching_diagnostics.R")

cat("\n=== Step 2: Propensity Score Matching ===\n")
source("matching_psm.R")

cat("\n=== Step 3: Coarsened Exact Matching ===\n")
source("matching_cem.R")

cat("\n=== Step 4: Mahalanobis Distance Matching ===\n")
source("matching_mahal.R")

# ------------------------------------------------------------------
# 2. Load all results
# ------------------------------------------------------------------
res_psm   <- readRDS("results_psm.rds")
res_cem   <- readRDS("results_cem.rds")
res_mahal <- readRDS("results_mahal.rds")
data      <- readRDS("panel_data.rds")

# ------------------------------------------------------------------
# 3. Full-sample baselines for comparison
# ------------------------------------------------------------------
outcome  <- "loggdp"
horizon  <- 10
controls <- "l(gdp_diff,1:2) + l(maxwind_95,1:2)"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

# Baseline: full sample + year FE
irf_baseline <- lp_panel_inter(
  data = data, outcome = outcome, main_var = "maxwind_95",
  interact_var = "period", controls = controls,
  horizon = horizon,
  fe = "countrycode[year] + countrycode[year2] + year",
  panel_id = panel_id, vcov_formula = vcov_fm
) %>%
  mutate(
    category = case_when(
      str_detect(category, "Pre")  ~ "Pre-1990",
      str_detect(category, "Post") ~ "Post-1990",
      TRUE ~ category
    ),
    method = "Full Sample"
  )

# Region x year FE (the "solution")
data_reg <- data %>%
  mutate(
    region = countrycode(countrycode, "iso3c", "region"),
    region = ifelse(is.na(region), "Other", region)
  )

irf_regionfe <- lp_panel_inter(
  data = data_reg, outcome = outcome, main_var = "maxwind_95",
  interact_var = "period", controls = controls,
  horizon = horizon,
  fe = "countrycode[year] + countrycode[year2] + region^year",
  panel_id = panel_id, vcov_formula = vcov_fm
) %>%
  mutate(
    category = case_when(
      str_detect(category, "Pre")  ~ "Pre-1990",
      str_detect(category, "Post") ~ "Post-1990",
      TRUE ~ category
    ),
    method = "Region×Year FE"
  )

# Extract matched IRFs with year FE only (the interesting comparison)
irf_psm   <- res_psm$irf   %>% filter(spec == "PSM + Year FE") %>% mutate(method = "PSM")
irf_cem   <- res_cem$irf   %>% filter(spec == "CEM + Year FE") %>% mutate(method = "CEM")
irf_mahal <- res_mahal$irf %>% filter(spec == "Mahalanobis + Year FE") %>% mutate(method = "Mahalanobis")

# ------------------------------------------------------------------
# 4. Combined IRF comparison plot
# ------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

all_irfs <- bind_rows(irf_baseline, irf_psm, irf_cem, irf_mahal, irf_regionfe) %>%
  mutate(method = factor(method,
                         levels = c("Full Sample", "PSM", "CEM",
                                    "Mahalanobis", "Region×Year FE")))

p_combined <- ggplot(all_irfs, aes(x = horizon, y = irf_mean,
                                    color = category, fill = category)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.3) +
  facet_wrap(~method, ncol = 5) +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_fill_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "IRF Comparison Across Matching Methods",
    subtitle = "All matching methods use year FE; rightmost panel uses region×year FE (no matching)",
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

ggsave("figures/irf_matching_comparison.png", p_combined,
       width = 16, height = 5, dpi = 300)
cat("\nSaved: figures/irf_matching_comparison.png\n")

# ------------------------------------------------------------------
# 5. Post-pre gap summary at h=5
# ------------------------------------------------------------------
cat("\n========== POST-PRE GAP AT h=5 ==========\n")
cat(sprintf("%-25s %10s %10s %10s\n", "Method", "Pre-1990", "Post-1990", "Gap (pp)"))
cat(paste(rep("-", 58), collapse = ""), "\n")

for (m in levels(all_irfs$method)) {
  h5 <- all_irfs %>% filter(method == m, horizon == 5)
  pre  <- h5 %>% filter(category == "Pre-1990")  %>% pull(irf_mean)
  post <- h5 %>% filter(category == "Post-1990") %>% pull(irf_mean)
  if (length(pre) > 0 && length(post) > 0) {
    cat(sprintf("%-25s %9.1f%% %9.1f%% %9.2f\n",
                m, pre * 100, post * 100, (post - pre) * 100))
  }
}

# ------------------------------------------------------------------
# 6. Balance improvement summary
# ------------------------------------------------------------------
cat("\n========== BALANCE IMPROVEMENT (|SMD|) ==========\n")

country_chars <- readRDS("country_chars_pre1990.rds")

compute_smd <- function(df, treated_var = "ever_treated") {
  vars <- c("mean_loggdp", "mean_growth", "mean_hc", "mean_pop")
  sapply(vars, function(v) {
    t_vals <- df[[v]][df[[treated_var]]]
    c_vals <- df[[v]][!df[[treated_var]]]
    pooled_sd <- sqrt((sd(t_vals, na.rm = TRUE)^2 + sd(c_vals, na.rm = TRUE)^2) / 2)
    abs(mean(t_vals, na.rm = TRUE) - mean(c_vals, na.rm = TRUE)) / pooled_sd
  })
}

# Before matching
smd_before <- compute_smd(country_chars)

# After PSM
psm_ids <- res_psm$matched_ids
psm_chars <- country_chars %>%
  filter(countrycode %in% psm_ids$countrycode)
smd_psm <- compute_smd(psm_chars)

# After CEM
cem_ids <- res_cem$matched_ids
cem_chars <- country_chars %>%
  filter(countrycode %in% cem_ids$countrycode)
smd_cem <- compute_smd(cem_chars)

# After Mahalanobis
mahal_ids <- res_mahal$matched_ids
mahal_chars <- country_chars %>%
  filter(countrycode %in% mahal_ids$countrycode)
smd_mahal <- compute_smd(mahal_chars)

var_labels <- c("Log GDP pc", "GDP Growth", "Human Capital", "Log Pop")
cat(sprintf("%-20s %10s %10s %10s %12s\n",
            "Variable", "Before", "PSM", "CEM", "Mahalanobis"))
cat(paste(rep("-", 65), collapse = ""), "\n")
for (i in seq_along(var_labels)) {
  cat(sprintf("%-20s %10.3f %10.3f %10.3f %12.3f\n",
              var_labels[i],
              smd_before[i], smd_psm[i], smd_cem[i], smd_mahal[i]))
}

# ------------------------------------------------------------------
# 7. Balance improvement plot
# ------------------------------------------------------------------
balance_df <- data.frame(
  variable = rep(var_labels, 4),
  smd = c(smd_before, smd_psm, smd_cem, smd_mahal),
  method = rep(c("Before", "PSM", "CEM", "Mahalanobis"), each = length(var_labels))
) %>%
  mutate(
    method = factor(method, levels = c("Before", "PSM", "CEM", "Mahalanobis")),
    variable = factor(variable, levels = rev(var_labels))
  )

p_balance <- ggplot(balance_df, aes(x = smd, y = variable, color = method)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "firebrick", alpha = 0.5) +
  scale_color_manual(values = c("Before" = "grey50", "PSM" = "steelblue",
                                 "CEM" = "forestgreen", "Mahalanobis" = "darkorange")) +
  labs(
    x = "Absolute Standardized Mean Difference",
    y = NULL,
    title = "Covariate Balance Before and After Matching",
    subtitle = "Dashed lines: 0.10 (good) and 0.25 (poor) thresholds",
    color = "Method"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "bottom"
  )

ggsave("figures/balance_improvement.png", p_balance,
       width = 8, height = 5, dpi = 300)
cat("\nSaved: figures/balance_improvement.png\n")

cat("\n=== All matching analysis complete ===\n")
