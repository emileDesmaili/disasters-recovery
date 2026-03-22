# =============================================================================
# Direct Damages Investigation — LP: GDP Response by Storm Category
# =============================================================================
# LP IRFs of log GDP on categorical wind shocks (Cat 1+/2+/3+).
#
# Controls: l(shock, 1:2) + l(gdp_diff, 1:2)
# FE:       countrycode[year] + countrycode[year2] + region^year
# Shock:    wind_cat1plus / wind_cat2plus / wind_cat3plus (m/s binary, no per-sqkm)
#
# Figures:
#   figures/irf_gdp_by_category.png
# =============================================================================

library(tidyverse)
library(fixest)
source("../../emileRegs.R")
setFixest_notes(FALSE)

data <- readRDS("panel_damages.rds")

myblue   <- "#1f78b4"
myred    <- "#e31a1c"
myorange <- "#f28e2b"

FE      <- "countrycode[year] + countrycode[year2] + region^year"
PANEL   <- c("countrycode", "year")
DK      <- DK ~ year
HORIZON <- 8

theme_nature <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      colour = "grey10", margin = margin(b = 5)),
      plot.subtitle    = element_text(colour = "grey45", size = base_size - 1,
                                      margin = margin(b = 10)),
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

# ── GDP LP by categorical shock ---------------------------------------------

cat_shocks <- c(
  wind_cat1plus = "Cat 1+",
  wind_cat2plus = "Cat 2+",
  wind_cat3plus = "Cat 3+"
)

irf_gdp_cat <- map2_dfr(names(cat_shocks), cat_shocks, function(sh, lab) {
  tryCatch(
    lp_panel(
      data         = data,
      outcome      = "loggdp",
      main_var     = sh,
      controls     = paste0("l(", sh, ", 1:2) + l(gdp_diff, 1:2)"),
      horizon      = HORIZON,
      fe           = FE,
      panel_id     = PANEL,
      vcov_formula = DK
    ) %>% mutate(shock = lab),
    error = function(e) { message("Skip ", sh, ": ", e$message); NULL }
  )
}) %>%
  mutate(shock = factor(shock, levels = cat_shocks))

pal_cat <- setNames(c(myblue, myorange, myred), cat_shocks)

p_gdp_cat <- ggplot(irf_gdp_cat, aes(x = horizon, colour = shock, fill = shock)) +
  geom_ribbon(aes(ymin = irf_down * 100, ymax = irf_up * 100),
              alpha = 0.12, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60",
             linewidth = 0.5) +
  geom_line(aes(y = irf_mean * 100), linewidth = 1.4) +
  scale_x_continuous(breaks = 0:HORIZON) +
  scale_colour_manual(values = pal_cat, name = NULL) +
  scale_fill_manual(values   = pal_cat, name = NULL) +
  labs(
    title = "GDP Response by Cyclone Intensity",
    x     = "Years after strike",
    y     = "Cumulative GDP growth (%)"
  ) +
  theme_nature() +
  guides(colour = guide_legend(nrow = 1))

ggsave("figures/irf_gdp_by_category.png", p_gdp_cat,
       width = 7.5, height = 7, dpi = 300)
message("Saved irf_gdp_by_category.png")

saveRDS(list(irf_gdp_cat = irf_gdp_cat), "results_comparison.rds")
message("Saved results_comparison.rds")
