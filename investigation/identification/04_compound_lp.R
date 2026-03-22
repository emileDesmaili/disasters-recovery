# =============================================================================
# Compound LP: State-Dependent GDP Response × Time Since Last Cyclone
# =============================================================================
# Model: f(loggdp,h) - l(loggdp,1) ~ w + w:time_since + controls | FE
#   where time_since = years elapsed since the previous cyclone event
#   (measured at the time of the current observation, not the current strike)
#
# Marginal effect ME(w | time_since = c) = β₁ + β₂·c, delta-method SE
#
# Five shock variables, all interacted with years_since_any:
#   maxwind       — raw max wind speed (m/s)
#   maxwind_resid — residualized maxwind (country linear trend + region×year removed)
#   wind_cat1plus — Cat 1+ indicator (maxwind ≥ 33 m/s)
#   wind_cat2plus — Cat 2+ indicator (maxwind ≥ 43 m/s)
#   wind_cat3plus — Cat 3+ indicator (maxwind ≥ 50 m/s)
#
# State variable: years_since_any — years since last maxwind > 0 (any TC passage)
#
# Evaluated at c ∈ {1, 3, 5} years:
#   1 year  → compounding: country was hit very recently
#   3 years → medium gap
#   5 years → isolated: long recovery window since last event
#
# FE: countrycode[year] + countrycode[year2] + region^year
# SE: Driscoll–Kraay (vcov = DK ~ year)
# =============================================================================

library(tidyverse)
library(fixest)
library(patchwork)
source("../../emileRegs.R")
setFixest_notes(FALSE)
setFixest_nthreads(1)

# ── Helper: years since last positive event -----------------------------------
# For each observation i, returns year[i] - year[last j < i with x[j] > 0].
# Returns NA if no prior positive event exists (first-ever strike in the panel).
# Note: at strike year t, records the gap since the *previous* strike — not 0.

time_since_last <- function(x, y) {
  result  <- rep(NA_real_, length(x))
  last_yr <- NA_real_
  for (i in seq_along(x)) {
    if (!is.na(last_yr))             result[i] <- y[i] - last_yr
    if (!is.na(x[i]) && x[i] > 0)   last_yr   <- y[i]
  }
  result
}

# ── Data -----------------------------------------------------------------------

data <- readRDS("../direct_damages/panel_damages.rds") %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  mutate(
    years_since_any = time_since_last(maxwind, year)
  ) %>%
  ungroup()

cat("years_since_any quantiles (all non-NA obs):\n")
print(round(quantile(data$years_since_any, c(.10,.25,.50,.75,.90), na.rm = TRUE), 1))

# ── Settings ------------------------------------------------------------------

FE      <- "countrycode[year] + countrycode[year2] + region^year"
PANEL   <- c("countrycode", "year")
DK      <- DK ~ year

# ── Residualise maxwind on country linear trend + region×year FE --------------
# Removes predictable cyclone-season variation and common regional trends.
# The residual captures idiosyncratic within-country wind shocks.
# Note: years_since_any is still defined on original maxwind > 0 (event timing).

resid_mod <- feols(maxwind ~ 1 | countrycode[year] + region^year,
                   data = data, panel.id = PANEL)
data$maxwind_resid <- residuals(resid_mod)

cat("\nmaxwind_resid summary:\n")
print(round(summary(data$maxwind_resid), 3))
cat("Correlation(maxwind, maxwind_resid):",
    round(cor(data$maxwind, data$maxwind_resid, use = "complete.obs"), 3), "\n\n")

HORIZON <- 10

eval_vals   <- c(1, 3, 5)
eval_labels <- c("1 yr (compounding)", "3 yrs", "5 yrs (isolated)")

# Red = compounding (recent prior event), gold = medium, blue = isolated
pal_since <- c(
  "1 yr (compounding)" = "#d7191c",
  "3 yrs"              = "#fdae61",
  "5 yrs (isolated)"   = "#2c7bb6"
)

myblue  <- "#2c7bb6"
myred   <- "#d7191c"
mygold  <- "#fdae61"

theme_nature <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 1,
                                   colour = "grey10", margin = margin(b = 4)),
      plot.subtitle = element_text(colour = "grey45", size = base_size - 2,
                                   margin = margin(b = 8)),
      axis.title    = element_text(colour = "grey15", size = base_size),
      axis.text     = element_text(colour = "grey25", size = base_size - 2),
      axis.line     = element_line(colour = "grey40", linewidth = 0.45),
      axis.ticks    = element_line(colour = "grey50", linewidth = 0.45),
      legend.position  = "top",
      legend.text      = element_text(colour = "grey20", size = base_size - 1),
      legend.key.size  = unit(1.3, "lines"),
      legend.title     = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(face = "bold", size = base_size - 2),
      plot.margin      = margin(10, 14, 8, 10)
    )
}

irf_plot <- function(df, title = NULL, sub = NULL) {
  ggplot(df, aes(x = horizon, colour = quantile, fill = quantile)) +
    geom_ribbon(aes(ymin = irf_down * 100, ymax = irf_up * 100),
                alpha = 0.10, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.5) +
    geom_line(aes(y = irf_mean * 100), linewidth = 1.6) +
    scale_x_continuous(breaks = seq(0, HORIZON, 2)) +
    scale_colour_manual(values = pal_since) +
    scale_fill_manual(values   = pal_since) +
    labs(title = title, subtitle = sub,
         x = "Years after strike", y = "Cumulative GDP growth (%)") +
    theme_nature() +
    guides(colour = guide_legend(nrow = 1))
}

# ── Shock variables to run ----------------------------------------------------

shock_vars <- list(
  maxwind = list(
    label    = "Max wind (raw, m/s)",
    controls = "l(maxwind, 1:2) + l(gdp_diff, 1:2)"
  ),
  maxwind_resid = list(
    label    = "Max wind (residualised, m/s)",
    controls = "l(maxwind_resid, 1:2) + l(gdp_diff, 1:2)"
  ),
  wind_cat1plus = list(
    label    = "Cat\u00a01+ indicator (\u226533\u00a0m/s)",
    controls = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2)"
  ),
  wind_cat2plus = list(
    label    = "Cat\u00a02+ indicator (\u226543\u00a0m/s)",
    controls = "l(wind_cat2plus, 1:2) + l(gdp_diff, 1:2)"
  ),
  wind_cat3plus = list(
    label    = "Cat\u00a03+ indicator (\u226550\u00a0m/s)",
    controls = "l(wind_cat3plus, 1:2) + l(gdp_diff, 1:2)"
  )
)

# ── Run lp_state_dep ----------------------------------------------------------

irf_all <- map_dfr(names(shock_vars), function(sv) {
  info <- shock_vars[[sv]]
  cat("Running state-dep LP, shock =", sv, "\n")
  tryCatch(
    lp_state_dep(
      data         = data,
      outcome      = "loggdp",
      main_var     = sv,
      state_var    = "years_since_any",
      controls     = info$controls,
      horizon      = HORIZON,
      fe           = FE,
      panel_id     = PANEL,
      vcov_formula = DK,
      eval_values  = eval_vals,
      eval_labels  = eval_labels
    ) %>% mutate(shock = sv, shock_label = info$label),
    error = function(e) { message("  Skip ", sv, ": ", e$message); NULL }
  )
}) %>%
  mutate(
    quantile    = factor(quantile, levels = eval_labels),
    shock       = factor(shock,    levels = names(shock_vars))
  )

# ── Summary table -------------------------------------------------------------

cat("\nResults at h = 5:\n")
irf_all %>%
  filter(horizon == 5) %>%
  select(shock, quantile, irf_mean, se, irf_down, irf_up) %>%
  mutate(across(c(irf_mean, se, irf_down, irf_up), ~ round(. * 100, 2))) %>%
  as.data.frame() %>%
  print()

# ── Individual figures --------------------------------------------------------

plots <- map(names(shock_vars), function(sv) {
  df <- irf_all %>% filter(shock == sv)
  if (nrow(df) == 0) return(NULL)
  irf_plot(
    df,
    title = paste0("GDP\u00d7Time-Since-Last: shock = ", shock_vars[[sv]]$label),
    sub   = "ME(shock | years since last TC) evaluated at \u03c4 \u2208 {1, 3, 5} yrs"
  )
})
names(plots) <- names(shock_vars)

iwalk(plots, function(p, sv) {
  if (is.null(p)) return(invisible(NULL))
  ggsave(paste0("figures/irf_compound_", sv, ".png"), p,
         width = 8, height = 6, dpi = 300)
  message("Saved irf_compound_", sv, ".png")
})

# ── Combined figure -----------------------------------------------------------

valid_plots <- compact(plots)
if (length(valid_plots) >= 2) {
  p_combined <- wrap_plots(valid_plots, ncol = 1) +
    plot_annotation(
      title = "GDP Response \u00d7 Time Since Last Cyclone: Five Shock Measures",
      theme = theme(
        plot.title = element_text(face = "bold", size = 16,
                                  colour = "grey10", hjust = 0.5,
                                  margin = margin(b = 6))
      )
    )

  ggsave("figures/irf_compound_all_shocks.png", p_combined,
         width = 8, height = 30, dpi = 300)
  message("Saved irf_compound_all_shocks.png")
}

# ── Save results --------------------------------------------------------------

saveRDS(irf_all, "results_compound_lp.rds")
message("Saved results_compound_lp.rds")
cat("\nDone.\n")
