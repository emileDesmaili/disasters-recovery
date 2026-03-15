# compound_figures.R
# Produce all figures for the compound shock analysis slides
# Run after compound_diagnostics.R, compound_estimation.R, compound_robustness.R

library(tidyverse)

# ------------------------------------------------------------------
# 1. Load all results
# ------------------------------------------------------------------
if (!file.exists("estimation_results.rds")) {
  source("compound_estimation.R")
}
if (!file.exists("robustness_results.rds")) {
  source("compound_robustness.R")
}

est_res <- readRDS("estimation_results.rds")
rob_res <- readRDS("robustness_results.rds")

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# Common theme
# ------------------------------------------------------------------
irf_theme <- theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "bottom"
  )

compound_colors <- c("Isolated" = "steelblue",
                      "Compound" = "firebrick",
                      "Isolated (no shock in t-1)" = "steelblue",
                      "Compound (shock in t-1)" = "firebrick",
                      "Difference" = "grey30")

# ------------------------------------------------------------------
# Figure 1: Main interaction result — Year FE vs Region×Year FE
# ------------------------------------------------------------------
plot_spec4 <- bind_rows(
  est_res$spec4_interact_yearfe,
  est_res$spec4r_interact_regionfe
) %>%
  mutate(spec = factor(spec, levels = c("Interaction: Treated + Year FE",
                                         "Interaction: Treated + Region×Year FE")))

p1 <- ggplot(plot_spec4, aes(x = horizon, y = irf_mean,
                              color = category, fill = category)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~spec, ncol = 2) +
  scale_color_manual(values = compound_colors) +
  scale_fill_manual(values = compound_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 0:10) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Compound vs Isolated Shocks: Interaction Approach",
    subtitle = "Shock × 1(shock in t-1) | Ever-treated countries only",
    color = NULL, fill = NULL
  ) +
  irf_theme

ggsave("figures/irf_compound_main.png", p1, width = 13, height = 5.5, dpi = 300)
cat("Saved: figures/irf_compound_main.png\n")

# ------------------------------------------------------------------
# Figure 2: Separate shock indicators (three FE specs)
# ------------------------------------------------------------------
plot_separate <- bind_rows(
  est_res$spec1_separate_yearfe,
  est_res$spec2_separate_regionfe,
  est_res$spec3_fullsample_regionfe
) %>%
  filter(term != "Difference") %>%
  mutate(
    category = case_when(
      term == "isolated_strict" ~ "Isolated",
      term == "compound_strict" ~ "Compound",
      TRUE ~ term
    ),
    spec = factor(spec, levels = c("Treated Only + Year FE",
                                    "Treated Only + Region×Year FE",
                                    "Full Sample + Region×Year FE"))
  )

p2 <- ggplot(plot_separate, aes(x = horizon, y = irf_mean,
                                 color = category, fill = category)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~spec, ncol = 3) +
  scale_color_manual(values = compound_colors) +
  scale_fill_manual(values = compound_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Separate Compound and Isolated Shock Indicators",
    subtitle = "Two binary regressors entered simultaneously | 95% CI (Driscoll-Kraay)",
    color = NULL, fill = NULL
  ) +
  irf_theme

ggsave("figures/irf_compound_separate.png", p2, width = 14, height = 5.5, dpi = 300)
cat("Saved: figures/irf_compound_separate.png\n")

# ------------------------------------------------------------------
# Figure 3: Difference (compound - isolated) across specs
# ------------------------------------------------------------------
plot_diff <- bind_rows(
  est_res$spec1_separate_yearfe,
  est_res$spec2_separate_regionfe,
  est_res$spec3_fullsample_regionfe
) %>%
  filter(term == "Difference") %>%
  mutate(spec = factor(spec, levels = c("Treated Only + Year FE",
                                         "Treated Only + Region×Year FE",
                                         "Full Sample + Region×Year FE")))

p3 <- ggplot(plot_diff, aes(x = horizon, y = irf_mean)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.2, fill = "grey50") +
  geom_line(linewidth = 1, color = "grey30") +
  geom_point(size = 1.5, color = "grey30") +
  facet_wrap(~spec, ncol = 3) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(
    x = "Horizon (years)",
    y = "Compound minus Isolated effect",
    title = "Difference: Compound - Isolated Shock Effects",
    subtitle = "Positive = compound shocks are less damaging; Negative = compound more damaging",
    color = NULL
  ) +
  irf_theme

ggsave("figures/irf_compound_difference.png", p3, width = 14, height = 5, dpi = 300)
cat("Saved: figures/irf_compound_difference.png\n")

# ------------------------------------------------------------------
# Figure 4: Cumulative exposure (0, 1, 2+ recent shocks)
# ------------------------------------------------------------------
cum_colors <- c("0 recent" = "steelblue", "1 recent" = "goldenrod",
                "2+ recent" = "firebrick")

p4 <- ggplot(rob_res$r2_cumulative, aes(x = horizon, y = irf_mean,
                                          color = category, fill = category)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = cum_colors) +
  scale_fill_manual(values = cum_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 0:10) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Effect by Recent Shock Exposure",
    subtitle = "Shock × (number of shocks in t-1 to t-3) | Ever-treated, year FE",
    color = "Shocks in past\n3 years", fill = "Shocks in past\n3 years"
  ) +
  irf_theme

ggsave("figures/irf_cumulative_exposure.png", p4, width = 8, height = 5.5, dpi = 300)
cat("Saved: figures/irf_cumulative_exposure.png\n")

# ------------------------------------------------------------------
# Figure 5: Drop always-compound countries
# ------------------------------------------------------------------
p5 <- ggplot(rob_res$r3_drop_always, aes(x = horizon, y = irf_mean,
                                           color = category, fill = category)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = compound_colors) +
  scale_fill_manual(values = compound_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 0:10) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Dropping Always-Compound Countries",
    subtitle = "Excluding countries with >90% compound shocks and 10+ total shocks",
    color = NULL, fill = NULL
  ) +
  irf_theme

ggsave("figures/irf_drop_always_compound.png", p5, width = 8, height = 5.5, dpi = 300)
cat("Saved: figures/irf_drop_always_compound.png\n")

# ------------------------------------------------------------------
# Figure 6: Compound × Period triple interaction
# ------------------------------------------------------------------
triple_colors <- c(
  "Isolated × Pre-1990"  = "steelblue",
  "Isolated × Post-1990" = "dodgerblue",
  "Compound × Pre-1990"  = "firebrick",
  "Compound × Post-1990" = "tomato"
)

# Clean up category names from the interaction output
r4_clean <- rob_res$r4_triple %>%
  mutate(
    category = str_replace_all(category, "_", " ") %>%
      str_replace("×", "×") %>%
      str_trim()
  )

# Try to map the categories to the expected names
r4_clean <- r4_clean %>%
  mutate(category = case_when(
    str_detect(category, "Isolated") & str_detect(category, "Pre")  ~ "Isolated × Pre-1990",
    str_detect(category, "Isolated") & str_detect(category, "Post") ~ "Isolated × Post-1990",
    str_detect(category, "Compound") & str_detect(category, "Pre")  ~ "Compound × Pre-1990",
    str_detect(category, "Compound") & str_detect(category, "Post") ~ "Compound × Post-1990",
    TRUE ~ category
  ))

p6 <- ggplot(r4_clean, aes(x = horizon, y = irf_mean,
                             color = category, fill = category)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = triple_colors) +
  scale_fill_manual(values = triple_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 0:10) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Compound × Period Interaction",
    subtitle = "Four-way split: compound/isolated × pre/post-1990",
    color = NULL, fill = NULL
  ) +
  irf_theme

ggsave("figures/irf_compound_period.png", p6, width = 9, height = 5.5, dpi = 300)
cat("Saved: figures/irf_compound_period.png\n")

# ------------------------------------------------------------------
# Figure 7: Continuous intensity × recent shock
# ------------------------------------------------------------------
p7 <- ggplot(rob_res$r1_continuous, aes(x = horizon, y = irf_mean,
                                          color = category, fill = category)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = compound_colors) +
  scale_fill_manual(values = compound_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 0:10) +
  labs(
    x = "Horizon (years)",
    y = "GDP response (per 1 SD wind speed)",
    title = "Continuous Wind Speed × Recent Shock",
    subtitle = "Continuous intensity measure interacted with 1(shock in t-1)",
    color = NULL, fill = NULL
  ) +
  irf_theme

ggsave("figures/irf_continuous_compound.png", p7, width = 8, height = 5.5, dpi = 300)
cat("Saved: figures/irf_continuous_compound.png\n")

# ------------------------------------------------------------------
# Figure 8: No-trends specification
# ------------------------------------------------------------------
p8 <- ggplot(rob_res$r5_no_trends, aes(x = horizon, y = irf_mean,
                                          color = category, fill = category)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = compound_colors) +
  scale_fill_manual(values = compound_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 0:10) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Country + Year FE Only (No Trends)",
    subtitle = "Simpler FE: more power but potentially biased by differential trends",
    color = NULL, fill = NULL
  ) +
  irf_theme

ggsave("figures/irf_no_trends.png", p8, width = 8, height = 5.5, dpi = 300)
cat("Saved: figures/irf_no_trends.png\n")

# ------------------------------------------------------------------
# Figure 9: Broad compound definition comparison
# ------------------------------------------------------------------
broad_colors <- c("Isolated (no shock in t-1 or t-2)" = "steelblue",
                   "Compound (shock in t-1 or t-2)" = "firebrick")

p9 <- ggplot(est_res$spec5_broad_yearfe, aes(x = horizon, y = irf_mean,
                                               color = category, fill = category)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = broad_colors) +
  scale_fill_manual(values = broad_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 0:10) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Broader Compound Definition: Shock in t-1 or t-2",
    subtitle = "2-year lookback window | Ever-treated, year FE",
    color = NULL, fill = NULL
  ) +
  irf_theme

ggsave("figures/irf_broad_compound.png", p9, width = 8, height = 5.5, dpi = 300)
cat("Saved: figures/irf_broad_compound.png\n")

cat("\nAll figures generated.\n")
