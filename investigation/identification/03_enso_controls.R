# =============================================================================
# ENSO Controls: Cat1+ IRF — Baseline vs Adding Niño-3.4, Temp, Precip
# =============================================================================
# 1. Baseline IRF: cat1+ → GDP (standard spec, full sample)
# 2. Robustness: add nino34, t, p as controls — does the cyclone effect change?
# 3. Heterogeneity: same comparison by Low / High climatology sub-sample
# 4. Scatter: x = maxwind, y = nino34 / t / p, faceted by climate variable
# =============================================================================

library(tidyverse)
library(fixest)
library(patchwork)
library(latex2exp)
source("../../emileRegs.R")
setFixest_notes(FALSE)
setFixest_nthreads(1)

# ── Data -----------------------------------------------------------------------

data <- readRDS("../direct_damages/panel_damages.rds") %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  mutate(
    W_bar      = mean(maxwind, na.rm = TRUE),
    nino34_std = (nino34 - mean(nino34, na.rm = TRUE)) /
                   sd(nino34, na.rm = TRUE),
    t_std      = (t - mean(t,      na.rm = TRUE)) / sd(t,      na.rm = TRUE),
    p_std      = (p - mean(p,      na.rm = TRUE)) / sd(p,      na.rm = TRUE)
  ) %>%
  ungroup()

# Split among ever-struck countries (W_bar > 0); never-struck get NA
data$clim2 <- ifelse(data$W_bar <= 0,  NA_character_,
                     ifelse(data$W_bar >= 20, "High climatology", "Low climatology"))
cat("clim2 distribution (threshold W_bar >= 20 m/s):\n")
print(table(data$clim2, useNA = "always"))

FE      <- "countrycode[year] + countrycode[year2] + region^year"
PANEL   <- c("countrycode", "year")
DK      <- DK ~ year
HORIZON <- 10

myblue  <- "#2c7bb6"
myred   <- "#d7191c"
mygold  <- "#fdae61"
mygreen <- "#1a9641"

pal_clim2 <- c("Low climatology" = myblue, "High climatology" = myred)

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

irf_plot <- function(df, pal, title = NULL, sub = NULL) {
  ggplot(df, aes(x = horizon, colour = spec, fill = spec)) +
    geom_ribbon(aes(ymin = irf_down * 100, ymax = irf_up * 100),
                alpha = 0.10, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.5) +
    geom_line(aes(y = irf_mean * 100), linewidth = 1.6) +
    scale_x_continuous(breaks = seq(0, HORIZON, 2)) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values   = pal) +
    labs(title = title, subtitle = sub,
         x = "Years after strike", y = "Cumulative GDP growth (%)") +
    theme_nature() +
    guides(colour = guide_legend(nrow = 1))
}

# =============================================================================
# 1. Figure: Scatter — x = maxwind, y = nino34 / t / p, faceted
# =============================================================================

scatter_data <- data %>%
  filter(maxwind > 0) %>%   # only strike-years for clarity
  select(countrycode, year, maxwind, nino34, t, p) %>%
  pivot_longer(c(nino34, t, p), names_to = "var", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(
    var_label = recode(var,
      nino34 = "Ni\u00f1o-3.4",
      t      = "Temperature anomaly",
      p      = "Precipitation anomaly"
    ),
    var_label = factor(var_label,
      levels = c("Ni\u00f1o-3.4", "Temperature anomaly", "Precipitation anomaly"))
  )

# Correlation annotations
cor_ann <- scatter_data %>%
  group_by(var_label) %>%
  summarise(
    r   = cor(maxwind, value, use = "complete.obs"),
    p   = cor.test(maxwind, value)$p.value,
    sig = case_when(p < 0.01 ~ "***", p < 0.05 ~ "**", p < 0.10 ~ "*", TRUE ~ ""),
    x   = max(maxwind, na.rm = TRUE) * 0.6,
    y   = max(value,   na.rm = TRUE) * 0.9,
    lbl = paste0("r = ", round(r, 3), sig),
    .groups = "drop"
  )

p_scatter <- ggplot(scatter_data, aes(x = maxwind, y = value)) +
  geom_point(alpha = 0.15, size = 0.9, shape = 19, colour = myblue) +
  geom_smooth(method = "loess", se = TRUE, colour = myred,
              fill = myred, alpha = 0.15, linewidth = 1.2) +
  geom_text(data = cor_ann,
            aes(x = x, y = y, label = lbl),
            colour = "grey20", size = 3.8, fontface = "bold") +
  facet_wrap(~ var_label, scales = "free_y", nrow = 1) +
  labs(
    title    = "Max Wind Speed vs Climate Anomalies (Strike-Years Only)",
    subtitle = "Pearson r annotated; *** p<0.01, ** p<0.05, * p<0.1",
    x = "Max wind speed (m/s)", y = "Climate anomaly"
  ) +
  theme_nature(base_size = 13)

ggsave("figures/clim_enso_scatter.png", p_scatter,
       width = 12, height = 7, dpi = 300)
message("Saved clim_enso_scatter.png")

# =============================================================================
# 2. Baseline vs climate-controlled IRF (full sample)
# =============================================================================

specs <- list(
  "Baseline"      = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2)",
  "+ Ni\u00f1o-3.4" = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2) + nino34_std",
  "+ Temp"        = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2) + t_std",
  "+ Precip"      = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2) + p_std",
  "+ All climate" = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2) + nino34_std + t_std + p_std"
)

pal_specs <- c(
  "Baseline"        = "grey30",
  "+ Ni\u00f1o-3.4" = myblue,
  "+ Temp"          = mygold,
  "+ Precip"        = mygreen,
  "+ All climate"   = myred
)

irf_specs <- map_dfr(names(specs), function(nm) {
  cat("Running spec:", nm, "\n")
  tryCatch(
    lp_panel(
      data         = data,
      outcome      = "loggdp",
      main_var     = "wind_cat1plus",
      controls     = specs[[nm]],
      horizon      = HORIZON,
      fe           = FE,
      panel_id     = PANEL,
      vcov_formula = DK
    ) %>% mutate(spec = nm),
    error = function(e) { message("Skip ", nm, ": ", e$message); NULL }
  )
}) %>%
  mutate(spec = factor(spec, levels = names(specs)))

cat("\nBaseline vs climate controls at h = 5:\n")
irf_specs %>% filter(horizon == 5) %>%
  select(spec, irf_mean, se) %>%
  mutate(across(c(irf_mean, se), ~round(. * 100, 3))) %>%
  print()

p_ctrl <- irf_plot(irf_specs, pal_specs,
  title = "Cat\u00a01+ IRF: Baseline vs Climate Controls",
  sub   = "Controls added one at a time: Ni\u00f1o-3.4, temperature, precipitation (standardised)")

ggsave("figures/irf_enso_controls.png", p_ctrl,
       width = 9, height = 7, dpi = 300)
message("Saved irf_enso_controls.png")

saveRDS(irf_specs, "results_enso_controls.rds")
message("Saved results_enso_controls.rds")
cat("\nDone.\n")
