# investigation_breakpoint.R
# Breakpoint robustness analysis: Is 1990 special, or does the pre/post
# divergence in cyclone-shock IRFs hold for other break years?
#
# Strategy:
#   1. Loop over candidate break years {1980, 1985, 1990, 1995, 2000}
#   2. For each, create a period split and run lp_panel_inter()
#   3. Extract h=5 coefficients for "Pre" and "Post" periods
#   4. Plot how the IRF gap varies across break years
#   5. Also run a decade-interacted specification to test gradual vs discrete

library(tidyverse)
library(fixest)
library(haven)

source("../emileRegs.R")

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
    year2    = year^2
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

# ------------------------------------------------------------------
# 2. LP settings (matching lp_gdp_shocks.R)
# ------------------------------------------------------------------
outcome  <- "loggdp"
horizon  <- 10
shock    <- "maxwind_95"
controls <- "l(gdp_diff,1:2) + l(maxwind_95,1:2)"
fe       <- "countrycode[year] + countrycode[year2] + year"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

# ------------------------------------------------------------------
# 3. Loop over candidate break years
# ------------------------------------------------------------------
break_years <- c(1980, 1985, 1990, 1995, 2000)

results_list <- list()

for (by in break_years) {

  cat("=== Break year:", by, "===\n")

  # Create period variable for this break year
  data$period <- ifelse(data$year <= by,
                        paste0("Pre-", by),
                        paste0("Post-", by))

  # Run LP with interaction
  irfs <- lp_panel_inter(
    data        = data,
    outcome     = outcome,
    main_var    = shock,
    interact_var = "period",
    controls    = controls,
    horizon     = horizon,
    fe          = fe,
    panel_id    = panel_id,
    vcov_formula = vcov_fm
  )

  # Tag with break year
  irfs$break_year <- by

  # Clean up category labels to generic "Pre" / "Post"
  irfs <- irfs %>%
    mutate(
      era = case_when(
        str_detect(category, "^Pre")  ~ "Pre",
        str_detect(category, "^Post") ~ "Post",
        TRUE ~ category
      )
    )

  results_list[[as.character(by)]] <- irfs
}

all_irfs <- bind_rows(results_list)

# ------------------------------------------------------------------
# 4. Extract h=5 coefficients and print summary table
# ------------------------------------------------------------------
h5 <- all_irfs %>%
  filter(horizon == 5) %>%
  select(break_year, era, irf_mean, se, irf_down, irf_up)

cat("\n\n========== HORIZON-5 IRF BY BREAK YEAR ==========\n")
h5_wide <- h5 %>%
  select(break_year, era, irf_mean, se) %>%
  pivot_wider(
    names_from  = era,
    values_from = c(irf_mean, se),
    names_glue  = "{era}_{.value}"
  ) %>%
  mutate(
    gap       = Post_irf_mean - Pre_irf_mean,
    gap_pct   = round(gap * 100, 3)
  )

print(h5_wide, width = 200)

# Also print the full table
cat("\nDetailed h=5 table:\n")
h5_print <- h5 %>%
  mutate(across(where(is.numeric) & !matches("break_year"),
                ~ round(. * 100, 3))) %>%
  arrange(break_year, era)
print(as.data.frame(h5_print))

# ------------------------------------------------------------------
# 5. Figure: pre/post IRF at h=5 across break years
# ------------------------------------------------------------------
plot_data <- h5 %>%
  mutate(
    irf_pct  = irf_mean * 100,
    down_pct = irf_down * 100,
    up_pct   = irf_up   * 100,
    era      = factor(era, levels = c("Pre", "Post"))
  )

p1 <- ggplot(plot_data, aes(x = break_year, y = irf_pct,
                              color = era, shape = era)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_pointrange(aes(ymin = down_pct, ymax = up_pct),
                  size = 0.8, linewidth = 0.7,
                  position = position_dodge(width = 1.5)) +
  scale_color_manual(values = c("Pre" = "steelblue", "Post" = "firebrick")) +
  labs(
    x     = "Break Year",
    y     = "IRF at h = 5 (%)",
    color = "Period", shape = "Period",
    title = "h=5 response by break year"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        plot.title    = element_text(face = "bold"))

# ------------------------------------------------------------------
# 6. Full IRF paths for all break years (supplementary panel)
# ------------------------------------------------------------------
all_irfs_plot <- all_irfs %>%
  mutate(
    irf_pct  = irf_mean * 100,
    down_pct = irf_down * 100,
    up_pct   = irf_up   * 100,
    era      = factor(era, levels = c("Pre", "Post")),
    break_label = paste0("Break: ", break_year)
  )

p2 <- ggplot(all_irfs_plot, aes(x = horizon, y = irf_pct,
                                 color = era, fill = era)) +
  geom_ribbon(aes(ymin = down_pct, ymax = up_pct),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~break_label, nrow = 1) +
  scale_color_manual(values = c("Pre" = "steelblue", "Post" = "firebrick")) +
  scale_fill_manual(values  = c("Pre" = "steelblue", "Post" = "firebrick")) +
  labs(
    x     = "Horizon (years)",
    y     = "GDP response (%)",
    color = "Period", fill = "Period",
    title = "Full IRF paths by break year"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

# ------------------------------------------------------------------
# 7. Decade-interacted specification
# ------------------------------------------------------------------
cat("\n\n========== DECADE-INTERACTED SPECIFICATION ==========\n")

data$decade <- paste0(floor(data$year / 10) * 10, "s")
cat("Decade counts:\n")
print(table(data$decade))

irfs_decade <- lp_panel_inter(
  data         = data,
  outcome      = outcome,
  main_var     = shock,
  interact_var = "decade",
  controls     = controls,
  horizon      = horizon,
  fe           = fe,
  panel_id     = panel_id,
  vcov_formula = vcov_fm
)

# Print h=5 by decade
h5_decade <- irfs_decade %>%
  filter(horizon == 5) %>%
  select(category, irf_mean, se, irf_down, irf_up) %>%
  mutate(across(c(irf_mean, se, irf_down, irf_up), ~ round(. * 100, 3)))

cat("\nHorizon-5 IRF by decade (%):\n")
print(as.data.frame(h5_decade))

# Plot decade IRFs
decade_plot_data <- irfs_decade %>%
  mutate(
    irf_pct  = irf_mean * 100,
    down_pct = irf_down * 100,
    up_pct   = irf_up   * 100,
    decade   = factor(category)
  )

p3 <- ggplot(decade_plot_data, aes(x = horizon, y = irf_pct)) +
  geom_ribbon(aes(ymin = down_pct, ymax = up_pct),
              fill = "steelblue", alpha = 0.15) +
  geom_line(color = "steelblue", linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~decade, nrow = 1) +
  labs(
    x     = "Horizon (years)",
    y     = "GDP response (%)",
    title = "Decade-interacted IRFs"
  ) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

# Bar chart of h=5 by decade
p4 <- ggplot(h5_decade, aes(x = category, y = irf_mean)) +
  geom_col(fill = "grey60", width = 0.6) +
  geom_errorbar(aes(ymin = irf_down, ymax = irf_up), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x     = "Decade",
    y     = "h=5 GDP response (%)",
    title = "h=5 response by decade (decade-interacted spec.)"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

# ------------------------------------------------------------------
# 8. Combine and save
# ------------------------------------------------------------------
library(patchwork)

combined <- (p1 + p4) / (p2) / (p3) +
  plot_layout(heights = c(1, 0.8, 0.8))

ggsave("figures/breakpoint_robustness.png", combined,
       width = 14, height = 14, dpi = 300)

cat("\n\nFigure saved to figures/breakpoint_robustness.png\n")
cat("Done.\n")
