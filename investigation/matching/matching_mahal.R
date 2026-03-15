# matching_mahal.R
# Mahalanobis Distance Matching with exact matching on continent
# then re-run the period-interacted LP

library(tidyverse)
library(fixest)
library(haven)
library(MatchIt)
library(countrycode)

source("../../emileRegs.R")
setFixest_notes(FALSE)

# ------------------------------------------------------------------
# 1. Load pre-computed data
# ------------------------------------------------------------------
if (!file.exists("country_chars_pre1990.rds") || !file.exists("panel_data.rds")) {
  cat("Running matching_diagnostics.R first...\n")
  source("matching_diagnostics.R")
}

country_chars <- readRDS("country_chars_pre1990.rds")
data <- readRDS("panel_data.rds")

# ------------------------------------------------------------------
# 2. Mahalanobis distance matching with exact on continent
# ------------------------------------------------------------------
match_df <- country_chars %>%
  mutate(treated = as.integer(ever_treated)) %>%
  select(countrycode, treated, mean_loggdp, mean_growth, mean_hc,
         mean_pop, continent)

cat("Fitting Mahalanobis distance matching...\n")

# Check continent distribution — exact matching only works if both groups
# have members in each continent
cat("\nContinent distribution:\n")
print(table(match_df$continent, match_df$treated))

# Only exact match on continents that have both treated and control
continent_both <- match_df %>%
  group_by(continent) %>%
  summarise(has_treated = any(treated == 1),
            has_control = any(treated == 0)) %>%
  filter(has_treated & has_control) %>%
  pull(continent)

cat(sprintf("\nContinents with both treated and control: %s\n",
            paste(continent_both, collapse = ", ")))

# Filter to continents with both groups for exact matching
match_df_filtered <- match_df %>%
  filter(continent %in% continent_both)

m_mahal <- matchit(
  treated ~ mean_loggdp + mean_growth + mean_hc + mean_pop,
  data = match_df_filtered,
  method = "nearest",
  distance = "mahalanobis",
  exact = ~continent,
  ratio = 2,
  replace = FALSE
)

cat("\n========== MAHALANOBIS MATCHING SUMMARY ==========\n")
print(summary(m_mahal))

matched_ids <- match.data(m_mahal)
matched_countries <- matched_ids$countrycode

cat(sprintf("\nMatched sample: %d countries (%d treated, %d control)\n",
            length(matched_countries),
            sum(matched_ids$treated == 1),
            sum(matched_ids$treated == 0)))

# ------------------------------------------------------------------
# 3. Re-run period-interacted LP on matched sample
# ------------------------------------------------------------------
data_matched <- data %>%
  filter(countrycode %in% matched_countries) %>%
  mutate(
    region = countrycode(countrycode, "iso3c", "region"),
    region = ifelse(is.na(region), "Other", region)
  )

outcome  <- "loggdp"
horizon  <- 10
controls <- "l(gdp_diff,1:2) + l(maxwind_95,1:2)"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

irf_mahal_yearfe <- lp_panel_inter(
  data = data_matched, outcome = outcome, main_var = "maxwind_95",
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
    spec = "Mahalanobis + Year FE"
  )

irf_mahal_regionfe <- lp_panel_inter(
  data = data_matched, outcome = outcome, main_var = "maxwind_95",
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
    spec = "Mahalanobis + Region×Year FE"
  )

# Full-sample baseline
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
    spec = "Full Sample + Year FE"
  )

# ------------------------------------------------------------------
# 4. Plot comparison
# ------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

plot_df <- bind_rows(irf_baseline, irf_mahal_yearfe, irf_mahal_regionfe) %>%
  mutate(spec = factor(spec, levels = c("Full Sample + Year FE",
                                         "Mahalanobis + Year FE",
                                         "Mahalanobis + Region×Year FE")))

p <- ggplot(plot_df, aes(x = horizon, y = irf_mean,
                          color = category, fill = category)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~spec, ncol = 3) +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_fill_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 0:10) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Mahalanobis Distance Matching: Period-Interacted IRFs",
    subtitle = "Mahalanobis on pre-1990 covariates, exact match on continent (ratio=2)",
    color = NULL, fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "bottom"
  )

ggsave("figures/irf_mahal.png", p, width = 14, height = 5.5, dpi = 300)
cat("\nSaved: figures/irf_mahal.png\n")

# Print h=5 gaps
cat("\nPost-Pre gaps at h=5:\n")
for (s in unique(plot_df$spec)) {
  h5 <- plot_df %>% filter(spec == s, horizon == 5)
  post <- h5 %>% filter(category == "Post-1990") %>% pull(irf_mean)
  pre  <- h5 %>% filter(category == "Pre-1990")  %>% pull(irf_mean)
  cat(sprintf("  %s: %.2f pp\n", s, (post - pre) * 100))
}

# Save results
saveRDS(list(irf = bind_rows(irf_mahal_yearfe, irf_mahal_regionfe),
             match_obj = m_mahal,
             matched_ids = matched_ids),
        "results_mahal.rds")
cat("Saved: results_mahal.rds\n")
cat("Done.\n")
