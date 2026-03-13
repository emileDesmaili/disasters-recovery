# lp_gdp_shocks.R
# Local projection analysis of natural disaster shocks on GDP
# Produces four figures: baseline and period-split IRFs, binary and continuous

library(tidyverse)
library(fixest)
library(haven)
source("emileRegs.R")

# ------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------
pwt <- read_stata("raw_data/pwt_clean.dta")
tcs <- read_stata("raw_data/ibtracs_clean.dta")
data <- pwt %>% left_join(tcs, by = c("year", "countrycode"))

# ------------------------------------------------------------------
# 2. Variable construction
# ------------------------------------------------------------------
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

# ------------------------------------------------------------------
# 3. LP settings
# ------------------------------------------------------------------
outcome  <- "loggdp"
horizon  <- 10
make_controls <- function(shock_var) {
  paste0("l(gdp_diff,1:2) + l(", shock_var, ",1:2)")
}
fe       <- "countrycode[year] + countrycode[year2] + year"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

shocks       <- c("maxwind_95", "energy_95", "nlands_95")
shocks_cont  <- c("maxwind_sd", "energy_sd", "nlands_sd")
shock_labels <- c("Max Wind Speed", "Cyclone Energy", "Landfalls")

# ------------------------------------------------------------------
# 4. Plotting helpers
# ------------------------------------------------------------------
irf_theme <- theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

plot_irf <- function(df, ylab = "GDP response", color_var = NULL) {
  if (is.null(color_var)) {
    p <- ggplot(df, aes(x = horizon, y = irf_mean)) +
      geom_ribbon(aes(ymin = irf_down, ymax = irf_up),
                  alpha = 0.2, fill = "grey50") +
      geom_line(linewidth = 1.1, color = "grey30")
  } else {
    df[[color_var]] <- factor(df[[color_var]],
                              levels = c("Pre-1990", "Post-1990"))
    p <- ggplot(df, aes(x = horizon, y = irf_mean,
                        color = .data[[color_var]],
                        fill  = .data[[color_var]])) +
      geom_ribbon(aes(ymin = irf_down, ymax = irf_up),
                  alpha = 0.15, color = NA) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("Pre-1990" = "steelblue",
                                    "Post-1990" = "firebrick")) +
      scale_fill_manual(values = c("Pre-1990" = "steelblue",
                                   "Post-1990" = "firebrick"))
  }

  p +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
    facet_wrap(~shock, ncol = 3) +
    labs(x = "Horizon (years)", y = ylab, color = "Period", fill = "Period") +
    scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
    irf_theme +
    theme(legend.position = if (!is.null(color_var)) "top" else "none")
}

# ------------------------------------------------------------------
# 5. Estimate and plot: binary shocks
# ------------------------------------------------------------------
irfs_baseline <- map_dfr(seq_along(shocks), ~ {
  lp_panel(data = data, outcome = outcome, main_var = shocks[.x],
           controls = make_controls(shocks[.x]), horizon = horizon, fe = fe,
           panel_id = panel_id, vcov_formula = vcov_fm) %>%
    mutate(shock = shock_labels[.x])
})

irfs_period <- map_dfr(seq_along(shocks), ~ {
  lp_panel_inter(data = data, outcome = outcome, main_var = shocks[.x],
                 interact_var = "period", controls = make_controls(shocks[.x]), horizon = horizon,
                 fe = fe, panel_id = panel_id, vcov_formula = vcov_fm) %>%
    mutate(shock = shock_labels[.x])
})

ggsave("figures/irf_all_shocks_95.png",
       plot_irf(irfs_baseline), width = 10, height = 4, dpi = 300)
ggsave("figures/irf_all_shocks_95_by_period.png",
       plot_irf(irfs_period, color_var = "category"), width = 10, height = 4, dpi = 300)

# ------------------------------------------------------------------
# 6. Estimate and plot: continuous shocks
# ------------------------------------------------------------------
irfs_cont <- map_dfr(seq_along(shocks_cont), ~ {
  lp_panel(data = data, outcome = outcome, main_var = shocks_cont[.x],
           controls = make_controls(shocks_cont[.x]), horizon = horizon, fe = fe,
           panel_id = panel_id, vcov_formula = vcov_fm) %>%
    mutate(shock = shock_labels[.x])
})

irfs_cont_period <- map_dfr(seq_along(shocks_cont), ~ {
  lp_panel_inter(data = data, outcome = outcome, main_var = shocks_cont[.x],
                 interact_var = "period", controls = make_controls(shocks_cont[.x]),
                 horizon = horizon, fe = fe, panel_id = panel_id,
                 vcov_formula = vcov_fm) %>%
    mutate(shock = shock_labels[.x])
})

ggsave("figures/irf_all_shocks_continuous.png",
       plot_irf(irfs_cont, ylab = "GDP response (per 1 SD)"),
       width = 10, height = 4, dpi = 300)
ggsave("figures/irf_all_shocks_continuous_by_period.png",
       plot_irf(irfs_cont_period, ylab = "GDP response (per 1 SD)",
                color_var = "category"),
       width = 10, height = 4, dpi = 300)

# ==================================================================
# 7. Nonlinear effects (following Nath, Ramey & Klenow 2025)
# ==================================================================

# --- 7a. Intensity bins: separate IRFs by storm severity ---
data <- data %>%
  mutate(
    intensity_bin = case_when(
      maxwind >= 50 * knot_to_ms ~ "Cat 3+",
      maxwind >= 33 * knot_to_ms ~ "Cat 1-2",
      maxwind >= 18 * knot_to_ms ~ "Trop. Storm",
      TRUE ~ "None"
    ),
    intensity_bin = factor(intensity_bin,
                           levels = c("None", "Trop. Storm", "Cat 1-2", "Cat 3+"))
  )

irfs_bins <- lp_panel_inter(
  data = data, outcome = outcome, main_var = "maxwind_sd",
  interact_var = "intensity_bin",
  controls = make_controls("maxwind_sd"),
  horizon = horizon, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
) %>%
  filter(category != "None") %>%
  mutate(category = factor(category,
                           levels = c("Trop. Storm", "Cat 1-2", "Cat 3+")))

p_bins <- ggplot(irfs_bins, aes(x = horizon, y = irf_mean,
                                color = category, fill = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Trop. Storm" = "goldenrod",
                                "Cat 1-2" = "steelblue",
                                "Cat 3+" = "firebrick")) +
  scale_fill_manual(values = c("Trop. Storm" = "goldenrod",
                               "Cat 1-2" = "steelblue",
                               "Cat 3+" = "firebrick")) +
  labs(x = "Horizon (years)", y = "GDP response (per 1 SD wind)",
       color = "Intensity", fill = "Intensity") +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  irf_theme + theme(legend.position = "top")

ggsave("figures/irf_intensity_bins.png", p_bins, width = 7, height = 5, dpi = 300)

# --- 7b. State-dependent model (Nath et al. style) ---
# Interact continuous shock with country's mean historical exposure
data <- data %>%
  group_by(countrycode) %>%
  mutate(mean_exposure = mean(maxwind[year >= 1970 & year <= 2014], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(mean_exposure_dm = mean_exposure - mean(mean_exposure))  # demean

# LP with shock + shock * mean_exposure
lp_state_dep <- function(data, outcome, main_var, state_var,
                         controls = NULL, horizon = 10,
                         fe = "countrycode + year",
                         panel_id = c("countrycode", "year"),
                         vcov_formula = DK ~ year) {
  rhs_controls <- if (!is.null(controls)) paste0(" + ", controls) else ""
  irf_list <- vector("list", horizon + 1)

  for (h in 0:horizon) {
    fml <- as.formula(paste0(
      "f(", outcome, ", ", h, ") - l(", outcome, ", 1) ~ ",
      main_var, " + ", main_var, ":", state_var,
      rhs_controls, " | ", fe
    ))
    mod <- feols(fml, data = data, panel.id = panel_id, vcov = vcov_formula)

    # Evaluate marginal effect at different state values
    beta_main  <- coef(mod)[main_var]
    inter_name <- grep(paste0(main_var, ":", state_var), names(coef(mod)), value = TRUE)
    if (length(inter_name) == 0)
      inter_name <- grep(paste0(state_var, ":", main_var), names(coef(mod)), value = TRUE)
    beta_inter <- coef(mod)[inter_name]

    # Evaluate at 10th, 50th, 90th percentile of state_var
    state_vals <- quantile(data[[state_var]], c(0.10, 0.50, 0.90), na.rm = TRUE)
    for (j in seq_along(state_vals)) {
      me <- beta_main + beta_inter * state_vals[j]
      # Delta method SE
      grad <- c(1, state_vals[j])
      se <- sqrt(t(grad) %*% vcov(mod)[c(main_var, inter_name),
                                        c(main_var, inter_name)] %*% grad)
      irf_list[[length(irf_list) + 1]] <- data.frame(
        horizon = h, irf_mean = me, se = as.numeric(se),
        irf_down = me - 1.96 * as.numeric(se),
        irf_up = me + 1.96 * as.numeric(se),
        quantile = paste0("p", c(10, 50, 90)[j]),
        state_val = state_vals[j]
      )
    }
  }
  bind_rows(irf_list)
}

irfs_state <- lp_state_dep(
  data = data, outcome = outcome, main_var = "maxwind_sd",
  state_var = "mean_exposure_dm", controls = make_controls("maxwind_sd"),
  horizon = horizon, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
) %>%
  mutate(quantile = factor(quantile, levels = c("p10", "p50", "p90"),
                           labels = c("Low exposure (p10)",
                                      "Median exposure (p50)",
                                      "High exposure (p90)")))

p_state <- ggplot(irfs_state, aes(x = horizon, y = irf_mean,
                                  color = quantile, fill = quantile)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Low exposure (p10)" = "steelblue",
                                "Median exposure (p50)" = "grey40",
                                "High exposure (p90)" = "firebrick")) +
  scale_fill_manual(values = c("Low exposure (p10)" = "steelblue",
                               "Median exposure (p50)" = "grey40",
                               "High exposure (p90)" = "firebrick")) +
  labs(x = "Horizon (years)", y = "GDP response (per 1 SD wind)",
       color = "Country mean exposure", fill = "Country mean exposure") +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  irf_theme + theme(legend.position = "top")

ggsave("figures/irf_state_dependent.png", p_state, width = 7, height = 5, dpi = 300)

cat("Done. Figures saved to figures/\n")
