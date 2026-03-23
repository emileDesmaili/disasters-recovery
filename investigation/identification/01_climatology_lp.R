# =============================================================================
# Identification: State-Dependent LP × Country Wind Climatology
# =============================================================================
# Specification:
#   Δ_h log GDP_{it} = β1·w_{it} + β2·(w_{it} × W̄_i) + controls + FE + ε
#
# where W̄_i = mean maxwind in STRIKE YEARS for country i (conditional
# climatology: average intensity experienced when a storm hits; NA for never-struck).
#
# The marginal effect ME(w | W̄ = c) = β1 + β2·c is evaluated at
# c ∈ {10, 40, 70} m/s to span low-, mid-, and high-climatology countries.
#
# Storm variables tested: maxwind, wind_cat1plus, wind_cat2plus,
#                         wind_cat3plus, maxwind_sqkm
#
# FE: countrycode[year] + countrycode[year2] + region^year
#     (identical to all other investigation scripts)
#
# Key identification note:
#   W̄_i is time-invariant → absorbed by country FE as a main effect, but
#   the interaction w_{it}×W̄_i is time-varying and therefore identified
#   within-country by variation in wind across years.
# =============================================================================

library(tidyverse)
library(fixest)
library(patchwork)
library(latex2exp)
source("../../emileRegs.R")
setFixest_notes(FALSE)

# ── Data ----------------------------------------------------------------------

data <- readRDS("../direct_damages/panel_damages.rds") %>%
  arrange(countrycode, year)

# Country climatology: mean maxwind across ALL years (including zeros)
data <- data %>%
  group_by(countrycode) %>%
  mutate(
    W_bar = mean(maxwind, na.rm = TRUE)   # m/s, time-invariant
  ) %>%
  ungroup()

cat("W\u0304 (unconditional climatology) quantiles across countries:\n")
print(round(quantile(
  distinct(data, countrycode, W_bar)$W_bar,
  c(0.10, 0.25, 0.50, 0.75, 0.90, 0.95), na.rm = TRUE
), 1))
cat("n countries with W̄ defined:", sum(!is.na(distinct(data, countrycode, W_bar)$W_bar)), "\n\n")

# ── Settings ------------------------------------------------------------------

FE      <- "countrycode[year] + countrycode[year2] + region^year"
PANEL   <- c("countrycode", "year")
DK      <- DK ~ year
HORIZON <- 10

eval_vals   <- c(20, 40)
eval_labels <- c("20 m/s", "40 m/s")

pal_clim <- c(
  "20 m/s" = "#2c7bb6",
  "40 m/s" = "#d7191c"
)

storm_vars <- list(
  maxwind      = list(label = "Max wind speed (m/s)",                    controls = "l(maxwind,      1:2) + l(gdp_diff, 1:2)"),
  wind_cat1plus = list(label = "Cat 1+ binary ($\\geq$33 m/s)",          controls = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2)"),
  wind_cat2plus = list(label = "Cat 2+ binary ($\\geq$43 m/s)",          controls = "l(wind_cat2plus, 1:2) + l(gdp_diff, 1:2)"),
  wind_cat3plus = list(label = "Cat 3+ binary ($\\geq$50 m/s)",          controls = "l(wind_cat3plus, 1:2) + l(gdp_diff, 1:2)"),
  maxwind_sqkm  = list(label = "Max wind $\\times$ area (m/s $\\cdot$ km$^2$)", controls = "l(maxwind_sqkm,  1:2) + l(gdp_diff, 1:2)")
)

theme_nature <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 1,
                                   colour = "grey10", margin = margin(b = 4)),
      plot.subtitle = element_text(colour = "grey45", size = base_size - 2,
                                   margin = margin(b = 8)),
      axis.title    = element_text(colour = "grey15", size = base_size),
      axis.text     = element_text(colour = "grey25", size = base_size - 1),
      axis.line     = element_line(colour = "grey40", linewidth = 0.45),
      axis.ticks    = element_line(colour = "grey50", linewidth = 0.45),
      legend.position  = "top",
      legend.text      = element_text(colour = "grey20", size = base_size - 1),
      legend.key.size  = unit(1.3, "lines"),
      legend.title     = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(12, 16, 10, 12)
    )
}

# ── Run state-dependent LP for each shock variable ----------------------------

irf_all <- map_dfr(names(storm_vars), function(sv) {
  info <- storm_vars[[sv]]
  cat("Running lp_state_dep for:", sv, "\n")
  tryCatch(
    lp_state_dep(
      data         = data,
      outcome      = "loggdp",
      main_var     = sv,
      state_var    = "W_bar",
      controls     = info$controls,
      horizon      = HORIZON,
      fe           = FE,
      panel_id     = PANEL,
      vcov_formula = DK,
      eval_values  = eval_vals,
      eval_labels  = eval_labels
    ) %>% mutate(shock = sv, shock_label = info$label),
    error = function(e) {
      message("  Skip ", sv, ": ", e$message)
      NULL
    }
  )
}) %>%
  mutate(
    quantile = factor(quantile, levels = eval_labels),
    shock    = factor(shock,    levels = names(storm_vars))
  )

cat("\nResults at h = 5:\n")
irf_all %>%
  filter(horizon == 5) %>%
  select(shock, quantile, irf_mean, se, irf_down, irf_up) %>%
  mutate(across(c(irf_mean, se, irf_down, irf_up), ~round(. * 100, 2))) %>%
  as.data.frame() %>%
  print()

# ── Plots: one per shock variable, lines by climatology ----------------------

plot_irf <- function(df, title, subtitle = NULL) {
  p <- ggplot(df, aes(x = horizon, colour = quantile, fill = quantile)) +
    geom_ribbon(aes(ymin = irf_down * 100, ymax = irf_up * 100),
                alpha = 0.10, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.5) +
    geom_line(aes(y = irf_mean * 100), linewidth = 1.4) +
    scale_x_continuous(breaks = 0:HORIZON) +
    scale_colour_manual(values = pal_clim, name = NULL) +
    scale_fill_manual(values   = pal_clim, name = NULL) +
    labs(title = title, subtitle = subtitle,
         x = "Years after strike", y = "Cumulative GDP growth (%)") +
    theme_nature() +
    guides(colour = guide_legend(nrow = 1))
  p
}

plots <- map(names(storm_vars), function(sv) {
  df <- irf_all %>% filter(shock == sv)
  if (nrow(df) == 0) return(NULL)
  plot_irf(
    df,
    title    = TeX(paste0("GDP Response: ", storm_vars[[sv]]$label)),
    subtitle = TeX("ME(wind | $\\bar{W}$) evaluated at $\\bar{W}$ $\\in$ \\{20, 40\\} m/s")
  )
})
names(plots) <- names(storm_vars)

# Save individual figures
iwalk(plots, function(p, sv) {
  if (is.null(p)) return(invisible(NULL))
  ggsave(paste0("figures/irf_clim_", sv, ".png"), p,
         width = 7.5, height = 7, dpi = 300)
  message("Saved irf_clim_", sv, ".png")
})

# ── Combined figure: all shocks in one patchwork panel -----------------------

valid_plots <- compact(plots)
if (length(valid_plots) >= 2) {
  p_combined <- wrap_plots(valid_plots, ncol = 3)

  ggsave("figures/irf_climatology_all.png", p_combined,
         width = 10, height = 6, dpi = 300)
  message("Saved irf_climatology_all.png")
}

# ── Save results --------------------------------------------------------------

saveRDS(irf_all, "results_climatology_lp.rds")
message("Saved results_climatology_lp.rds")
cat("\nDone.\n")
