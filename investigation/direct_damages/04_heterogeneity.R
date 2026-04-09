# =============================================================================
# Direct Damages Investigation — Heterogeneous GDP Responses
# =============================================================================
# Two approaches, both using shared helper functions from emileRegs.R:
#
#  1. CONTINUOUS INTERACTION (state-dependent LP)
#     lp_state_dep(): feols(Δ_h loggdp ~ wind + wind:ln_damage_lagGDP + ctrls | FE)
#     Marginal effect evaluated at p25 / p75 of ln_damage_lagGDP.
#
#  2. DISCRETE INTERACTION (categorical vulnerability binary)
#     lp_panel_inter(): feols(Δ_h loggdp ~ i(damage_binary, wind) + ctrls | FE)
#     Gives separate slope for Low / High vulnerability directly.
# =============================================================================

library(tidyverse)
library(fixest)
source("../../emileRegs.R")
setFixest_notes(FALSE)

data <- readRDS("panel_damages.rds")

myblue   <- "#1f78b4"
myred    <- "#e31a1c"
myorange <- "#f28e2b"
mygreen  <- "#2ca25f"
mypurple <- "#756bb1"

FE      <- "countrycode[year] + countrycode[year2] + region^year"
PANEL   <- c("countrycode", "year")
DK      <- DK ~ year
HORIZON <- 8

theme_nature <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      colour = "grey10", margin = margin(b = 5)),
      axis.title       = element_text(colour = "grey15", size = base_size),
      axis.text        = element_text(colour = "grey25", size = base_size - 1),
      axis.line        = element_line(colour = "grey40", linewidth = 0.45),
      axis.ticks       = element_line(colour = "grey50", linewidth = 0.45),
      legend.position  = "top",
      legend.text      = element_text(colour = "grey20", size = base_size - 1),
      legend.key.size  = unit(1.3, "lines"),
      legend.title     = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(12, 16, 10, 12)
    )
}

# ── 1. Continuous interaction: lp_state_dep() --------------------------------
# ME(wind | damage=d) = β_wind + β_inter × d, evaluated at p25 / p75.

irf_state <- lp_state_dep(
  data           = data,
  outcome        = "loggdp",
  main_var       = "wind_cat1plus",
  state_var      = "ln_damage_lagGDP",
  controls       = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2)",
  horizon        = HORIZON,
  fe             = FE,
  panel_id       = PANEL,
  vcov_formula   = DK,
  eval_quantiles = c(0.25, 0.75),
  eval_labels    = c("Low damage (p25)", "High damage (p75)")
)

if (nrow(irf_state) > 0) {
  pal_state <- c("Low damage (p25)" = myblue, "High damage (p75)" = myred)

  irf_state <- irf_state %>%
    mutate(quantile = factor(quantile, levels = names(pal_state)))

  p_state <- ggplot(irf_state,
                    aes(x = horizon, colour = quantile, fill = quantile)) +
    geom_ribbon(aes(ymin = irf_down * 100, ymax = irf_up * 100),
                alpha = 0.10, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60",
               linewidth = 0.5) +
    geom_line(aes(y = irf_mean * 100), linewidth = 1.4) +
    scale_x_continuous(breaks = 0:HORIZON) +
    scale_colour_manual(values = pal_state, name = NULL) +
    scale_fill_manual(values   = pal_state, name = NULL) +
    labs(
      title = "GDP Response to Cat 1+ by log(damage/GDP[t-1])",
      x     = "Years after strike",
      y     = "Cumulative GDP growth (%)"
    ) +
    theme_nature() +
    guides(colour = guide_legend(nrow = 1))

  ggsave("figures/irf_gdp_state_dep_damage.png", p_state,
         width = 7.5, height = 7, dpi = 300)
  message("Saved irf_gdp_state_dep_damage.png")
}

# ── 2. Discrete interaction: lp_panel_inter() --------------------------------
# i(damage_binary, wind_cat1plus) gives a separate wind slope for each
# vulnerability level — Low and High — from a single pooled regression.

pal_vuln <- c("Low vulnerability" = myblue, "High vulnerability" = myred)

irf_by_vuln <- lp_panel_inter(
  data         = data,
  outcome      = "loggdp",
  main_var     = "wind_cat1plus",
  interact_var = "damage_binary",
  controls     = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2)",
  horizon      = HORIZON,
  fe           = FE,
  panel_id     = PANEL,
  vcov_formula = DK
) %>%
  rename(group = category) %>%
  mutate(group = factor(group, levels = names(pal_vuln)))

if (nrow(irf_by_vuln) > 0) {
  p_vuln <- ggplot(irf_by_vuln,
                   aes(x = horizon, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = irf_down * 100, ymax = irf_up * 100),
                alpha = 0.12, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60",
               linewidth = 0.5) +
    geom_line(aes(y = irf_mean * 100), linewidth = 1.4) +
    scale_x_continuous(breaks = 0:HORIZON) +
    scale_colour_manual(values = pal_vuln, name = NULL) +
    scale_fill_manual(values   = pal_vuln, name = NULL) +
    labs(
      title = "GDP Response by Damage Vulnerability (Cat 1+)",
      x     = "Years after strike",
      y     = "Cumulative GDP growth (%)"
    ) +
    theme_nature() +
    guides(colour = guide_legend(nrow = 1))

  ggsave("figures/irf_gdp_by_damage_vuln.png", p_vuln,
         width = 7.5, height = 7, dpi = 300)
  message("Saved irf_gdp_by_damage_vuln.png")
}

saveRDS(list(
  irf_state   = irf_state,
  irf_by_vuln = irf_by_vuln
), "results_heterogeneity.rds")
message("Saved results_heterogeneity.rds")
