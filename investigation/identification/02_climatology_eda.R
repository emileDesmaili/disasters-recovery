# =============================================================================
# Identification EDA: High-Climatology Countries, Wind Time Series,
#                     Autocorrelation, Compound Events, and ENSO Correlates
# =============================================================================

library(tidyverse)
library(fixest)
library(patchwork)
source("../../emileRegs.R")
setFixest_notes(FALSE)

# ── Data -----------------------------------------------------------------------
# Arrange by countrycode, year BEFORE all group-by mutates so lag() is correct

data <- readRDS("../direct_damages/panel_damages.rds") %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  mutate(
    W_bar       = mean(maxwind, na.rm = TRUE),
    struck      = as.integer(maxwind > 0),
    struck_lag1 = lag(struck, 1)              # safe: data sorted by year
  ) %>%
  ungroup()

# Country-level summary — one row per country, guaranteed
clim_cty <- data %>%
  group_by(countrycode) %>%
  summarise(
    W_bar     = first(W_bar),
    continent = first(continent),
    region    = first(region),
    .groups   = "drop"
  ) %>%
  filter(!is.na(W_bar)) %>%
  arrange(desc(W_bar)) %>%
  mutate(
    rank     = row_number(),
    # Binary split: High = W_bar >= 33 m/s (Cat1+ threshold), Low = 0 < W_bar < 33
    clim2    = if_else(W_bar <= 0,  NA_character_,
               if_else(W_bar >= 20, "High climatology", "Low climatology")),
    clim2    = factor(clim2, levels = c("Low climatology", "High climatology")),
    clim_grp = clim2   # alias used in downstream compound-event plots
  )

cat("Top 20 countries by climatology:\n")
print(clim_cty %>% select(rank, countrycode, W_bar, continent) %>% head(20), n = 20)

cat("\nW_bar quantiles:\n")
print(round(quantile(clim_cty$W_bar, c(.10,.25,.50,.75,.90,.95), na.rm=TRUE), 2))

# ── Theme & colours -----------------------------------------------------------

myblue  <- "#2c7bb6"
myred   <- "#d7191c"
mygold  <- "#fdae61"
mygreen <- "#1a9641"

pal_clim2 <- c("High climatology" = myred, "Low climatology" = myblue)
pal_grp   <- pal_clim2   # alias

theme_nature <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      colour = "grey10", margin = margin(b = 4)),
      plot.subtitle    = element_text(colour = "grey45", size = base_size - 2,
                                      margin = margin(b = 8)),
      axis.title       = element_text(colour = "grey15", size = base_size),
      axis.text        = element_text(colour = "grey25", size = base_size - 2),
      axis.line        = element_line(colour = "grey40", linewidth = 0.45),
      axis.ticks       = element_line(colour = "grey50", linewidth = 0.45),
      legend.position  = "top",
      legend.text      = element_text(colour = "grey20", size = base_size - 1),
      legend.key.size  = unit(1.2, "lines"),
      legend.title     = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(face = "bold", size = base_size - 2),
      plot.margin      = margin(10, 14, 8, 10)
    )
}

# =============================================================================
# Figure 1: Top-20 countries by climatology
# =============================================================================

top20 <- clim_cty %>% head(20) %>%
  mutate(cc = fct_reorder(countrycode, W_bar))

p_top20 <- ggplot(top20, aes(x = cc, y = W_bar, fill = continent)) +
  geom_col(alpha = 0.85, width = 0.75) +
  geom_hline(yintercept = quantile(clim_cty$W_bar, 0.90, na.rm = TRUE),
             linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  annotate("text", x = 1.2, y = quantile(clim_cty$W_bar, 0.90, na.rm = TRUE) + 0.3,
           label = "p90", hjust = 0, colour = "grey40", size = 3.8) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Top-20 Countries by Wind Climatology",
    subtitle = expression(bar(W)[i]~"= mean maxwind across all years (m/s)"),
    x = NULL, y = expression(bar(W)[i]~"(m/s)")
  ) +
  theme_nature() +
  theme(legend.position = "right")

ggsave("figures/clim_top20.png", p_top20, width = 8, height = 7, dpi = 300)
message("Saved clim_top20.png")

# =============================================================================
# Figure 2: Wind time series for top-6 countries
# =============================================================================

top6 <- clim_cty %>% head(6) %>% pull(countrycode)

ts_data <- data %>%
  filter(countrycode %in% top6) %>%
  mutate(
    label = paste0(countrycode, " (", round(W_bar, 1), " m/s)"),
    label = fct_reorder(label, -W_bar)
  )

p_ts <- ggplot(ts_data, aes(x = year, y = maxwind)) +
  geom_col(aes(fill = maxwind > 0), width = 0.9, alpha = 0.75) +
  geom_hline(aes(yintercept = W_bar), linetype = "dashed",
             colour = myred, linewidth = 0.6) +
  scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = myblue), guide = "none") +
  scale_x_continuous(breaks = c(1970, 1985, 2000, 2014)) +
  facet_wrap(~ label, ncol = 3, scales = "free_y") +
  labs(
    title    = "Annual Max Wind Speed: Top-6 High-Climatology Countries",
    subtitle = expression("Bars = annual max wind (m/s); dashed red = country climatology"~bar(W)[i]),
    x = NULL, y = "Max wind (m/s)"
  ) +
  theme_nature(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

ggsave("figures/clim_timeseries.png", p_ts, width = 13, height = 7, dpi = 300)
message("Saved clim_timeseries.png")

# =============================================================================
# Figure 3: Autocorrelation of strike indicator — overall across all countries
# =============================================================================

data_acf <- data %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  mutate(
    s_l1 = lag(struck, 1),
    s_l2 = lag(struck, 2),
    s_l3 = lag(struck, 3),
    s_l4 = lag(struck, 4),
    s_l5 = lag(struck, 5)
  ) %>%
  summarise(
    n_obs = n(),
    acf0  = 1,
    acf1  = if (n() > 10) cor(struck, s_l1, use = "complete.obs") else NA_real_,
    acf2  = if (n() > 10) cor(struck, s_l2, use = "complete.obs") else NA_real_,
    acf3  = if (n() > 10) cor(struck, s_l3, use = "complete.obs") else NA_real_,
    acf4  = if (n() > 10) cor(struck, s_l4, use = "complete.obs") else NA_real_,
    acf5  = if (n() > 10) cor(struck, s_l5, use = "complete.obs") else NA_real_,
    .groups = "drop"
  )

acf_long <- data_acf %>%
  pivot_longer(starts_with("acf"), names_to = "lag",
               names_prefix = "acf", values_to = "acf") %>%
  mutate(lag = as.integer(lag))

acf_summary <- acf_long %>%
  group_by(lag) %>%
  summarise(
    mean_acf = mean(acf, na.rm = TRUE),
    se_acf   = sd(acf, na.rm = TRUE) / sqrt(sum(!is.na(acf))),
    .groups  = "drop"
  )

cat("\nMean ACF by lag (all countries):\n")
print(acf_summary)

p_acf <- ggplot(acf_summary, aes(x = lag, y = mean_acf)) +
  geom_ribbon(aes(ymin = mean_acf - 1.96 * se_acf,
                  ymax = mean_acf + 1.96 * se_acf),
              alpha = 0.15, fill = myblue, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.5) +
  geom_line(colour = myblue, linewidth = 1.6) +
  scale_x_continuous(breaks = 0:5) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    title    = "Autocorrelation of Annual Strike Indicator",
    subtitle = "Mean country-level ACF(struck, lag k); ribbon = \u00b11.96 SE across countries",
    x = "Lag (years)", y = "Mean ACF"
  ) +
  theme_nature()

ggsave("figures/clim_acf.png", p_acf, width = 8, height = 6, dpi = 300)
message("Saved clim_acf.png")

# =============================================================================
# Figure 4: Compound events — P(struck | struck last year) by W_bar decile
# =============================================================================

compound <- data %>%
  filter(!is.na(struck_lag1)) %>%
  left_join(clim_cty %>% select(countrycode, clim_grp), by = "countrycode") %>%
  filter(!is.na(clim_grp)) %>%
  mutate(W_decile = ntile(W_bar, 10))

p_compound_cond <- compound %>%
  group_by(W_decile) %>%
  summarise(
    W_bar_mean    = mean(W_bar, na.rm = TRUE),
    p_given_hit   = mean(struck[struck_lag1 == 1], na.rm = TRUE),
    p_given_nohit = mean(struck[struck_lag1 == 0], na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  pivot_longer(c(p_given_hit, p_given_nohit),
               names_to = "condition", values_to = "prob") %>%
  mutate(condition = recode(condition,
    p_given_hit   = "P(struck | struck t\u22121)",
    p_given_nohit = "P(struck | not struck t\u22121)"
  ))

p_cond <- ggplot(p_compound_cond,
                 aes(x = W_bar_mean, y = prob, colour = condition, shape = condition)) +
  geom_smooth(aes(fill = condition), method = "loess", span = 1.2,
              alpha = 0.10, linewidth = 1.2, se = TRUE) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(myred, myblue)) +
  scale_fill_manual(values   = c(myred, myblue)) +
  scale_shape_manual(values  = c(19, 17)) +
  labs(
    title    = "Compound Event Risk by Climatology",
    subtitle = "Conditional strike probability by \u1e84 decile",
    x        = expression(bar(W)[i]~"decile mean (m/s)"),
    y        = "Probability of strike"
  ) +
  theme_nature()

ggsave("figures/clim_compound_cond.png", p_cond, width = 8, height = 6, dpi = 300)
message("Saved clim_compound_cond.png")

# =============================================================================
# Figure 5: Back-to-back strike rate by clim group
# =============================================================================

runs <- data %>%
  left_join(clim_cty %>% select(countrycode, clim_grp), by = "countrycode") %>%
  filter(!is.na(clim_grp), !is.na(struck_lag1)) %>%
  group_by(countrycode, clim_grp) %>%
  summarise(
    n_struck       = sum(struck, na.rm = TRUE),
    n_backtoback   = sum(struck == 1 & struck_lag1 == 1, na.rm = TRUE),
    pct_backtoback = 100 * n_backtoback / max(n_struck, 1),
    .groups        = "drop"
  )

cat("\nBack-to-back strike rates by climatology group:\n")
runs %>%
  group_by(clim_grp) %>%
  summarise(median_pct = median(pct_backtoback, na.rm = TRUE),
            mean_pct   = mean(pct_backtoback, na.rm = TRUE),
            n          = n()) %>% print()

p_runs <- ggplot(runs %>% filter(n_struck > 0),
                 aes(x = clim_grp, y = pct_backtoback, fill = clim_grp)) +
  geom_violin(alpha = 0.35, colour = NA) +
  geom_boxplot(width = 0.22, outlier.shape = NA,
               colour = "grey30", linewidth = 0.7) +
  geom_jitter(aes(colour = clim_grp), width = 0.07, size = 1.6,
              alpha = 0.55, shape = 19) +
  scale_fill_manual(values   = pal_grp, guide = "none") +
  scale_colour_manual(values = pal_grp, guide = "none") +
  labs(
    title    = "Back-to-Back Strike Rate by Climatology Group",
    subtitle = "% of strike-years following another strike-year",
    x = NULL, y = "Back-to-back rate (%)"
  ) +
  theme_nature()

ggsave("figures/clim_backtoback.png", p_runs, width = 7.5, height = 6, dpi = 300)
message("Saved clim_backtoback.png")

# =============================================================================
# Figure 6: Combined compound panel
# =============================================================================

p_combined_compound <- (p_acf | p_cond) / p_runs +
  plot_annotation(
    title    = "Wind Climatology and Compound Cyclone Exposure",
    subtitle = paste0(
      "Do high-climatology countries face more serially correlated storms?\n",
      "Top: ACF of annual strike indicator; Bottom-left: conditional strike probs; ",
      "Bottom-right: back-to-back rates"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 16, colour = "grey10", hjust = 0.5),
      plot.subtitle = element_text(colour = "grey40", size = 12, hjust = 0.5,
                                   margin = margin(b = 8))
    )
  )

ggsave("figures/clim_compound_combined.png", p_combined_compound,
       width = 14, height = 12, dpi = 300)
message("Saved clim_compound_combined.png")

# =============================================================================
# Figure 7: Average inter-strike interval vs average max wind speed
#           (ever-treated countries only)
# =============================================================================

inter_strike <- data %>%
  filter(W_bar > 0) %>%                    # ever-treated countries only
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  summarise(
    W_bar      = first(W_bar),
    continent  = first(continent),
    n_years    = n(),
    n_strikes  = sum(maxwind > 0, na.rm = TRUE),
    # mean gap between consecutive strike years
    mean_gap   = {
      strike_yrs <- year[maxwind > 0]
      if (length(strike_yrs) >= 2) mean(diff(sort(strike_yrs))) else NA_real_
    },
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_gap))

p_gap <- ggplot(inter_strike, aes(x = W_bar, y = mean_gap)) +
  geom_point(aes(colour = continent), alpha = 0.65, size = 2.2, shape = 19) +
  geom_smooth(method = "loess", se = TRUE,
              colour = myred, fill = myred, alpha = 0.15, linewidth = 1.4) +
  scale_colour_brewer(palette = "Set2") +
  labs(
    title    = "Inter-Strike Interval vs Country Wind Climatology",
    subtitle = "Ever-treated countries only; LOESS fit with 95% CI",
    x        = expression(bar(W)[i]~"= mean maxwind across all years (m/s)"),
    y        = "Mean years between consecutive strikes",
    colour   = NULL
  ) +
  theme_nature() +
  theme(legend.position = "right")

ggsave("figures/clim_gap_vs_wbar.png", p_gap, width = 9, height = 6, dpi = 300)
message("Saved clim_gap_vs_wbar.png")

# =============================================================================
# Figure 8: ENSO / T / P distributions — strike vs non-strike years
#           Ever-treated countries only (W_bar > 0)
# =============================================================================

enso_box <- data %>%
  filter(W_bar > 0) %>%
  mutate(strike = if_else(maxwind > 0, "Strike year", "No strike")) %>%
  select(strike, nino34, t, p) %>%
  pivot_longer(c(nino34, t, p), names_to = "var", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(
    var_label = recode(var,
      nino34 = "Ni\u00f1o-3.4",
      t      = "Temperature",
      p      = "Precipitation"
    ),
    var_label = factor(var_label,
      levels = c("Ni\u00f1o-3.4", "Temperature", "Precipitation")),
    strike = factor(strike, levels = c("No strike", "Strike year"))
  )

# Welch t-test annotations
ttest_ann <- enso_box %>%
  group_by(var_label) %>%
  summarise(
    p_val = tryCatch(
      t.test(value[strike == "Strike year"], value[strike == "No strike"])$p.value,
      error = function(e) NA_real_
    ),
    y_pos = quantile(value, 0.97, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sig = case_when(p_val < 0.01 ~ "***", p_val < 0.05 ~ "**",
                    p_val < 0.10 ~ "*",  TRUE ~ "ns"),
    label = paste0("p = ", formatC(p_val, format = "f", digits = 3), "  ", sig)
  )

p_box <- ggplot(enso_box, aes(x = strike, y = value, fill = strike)) +
  geom_boxplot(width = 0.45, outlier.shape = NA,
               linewidth = 0.55, alpha = 0.75) +
  geom_text(data = ttest_ann,
            aes(x = 1.5, y = y_pos, label = label, fill = NULL),
            inherit.aes = FALSE, size = 3.8, colour = "grey25", fontface = "bold") +
  facet_wrap(~ var_label, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("No strike" = "grey70", "Strike year" = myblue),
                    guide = "none") +
  labs(
    title    = "Climate Anomalies: Strike vs Non-Strike Years",
    subtitle = "Ever-treated countries only (W\u0304 > 0); Welch t-test p-value annotated",
    x = NULL, y = "Anomaly"
  ) +
  theme_nature(base_size = 13)

ggsave("figures/clim_enso_boxplot.png", p_box, width = 11, height = 5, dpi = 300)
message("Saved clim_enso_boxplot.png")

cat("\nDone.\n")
