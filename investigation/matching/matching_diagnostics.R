# matching_diagnostics.R
# Characterize treated vs never-treated on pre-treatment covariates
# Motivates why matching is needed for control group construction

library(tidyverse)
library(haven)
library(countrycode)

# ------------------------------------------------------------------
# 1. Load and construct data
# ------------------------------------------------------------------
pwt  <- read_stata("../../raw_data/pwt_clean.dta")
tcs  <- read_stata("../../raw_data/ibtracs_clean.dta")
data <- pwt %>% left_join(tcs, by = c("year", "countrycode"))

knot_to_ms <- 0.514444

data <- data %>%
  mutate(
    maxwind  = replace_na(max_ann_wind_i_nob, 0) * knot_to_ms,
    loggdp   = ln_real_gdp_usd_pc,
    gdp_diff = growth_real_gdp_usd_pc,
    year2    = year^2,
    period   = ifelse(year <= 1990, "Pre-1990", "Post-1990")
  ) %>%
  filter(year >= 1969)

data <- data %>%
  mutate(
    maxwind_p95 = quantile(maxwind[year >= 1970 & year <= 2014], 0.95, na.rm = TRUE)
  ) %>%
  mutate(
    maxwind_95 = as.integer(maxwind_p95 > 0 & maxwind >= maxwind_p95)
  )

# Classify treated vs never-treated
treated_countries <- data %>%
  filter(maxwind_95 == 1) %>%
  pull(countrycode) %>%
  unique()

data <- data %>%
  mutate(
    ever_treated = countrycode %in% treated_countries,
    region = countrycode(countrycode, "iso3c", "region"),
    region = ifelse(is.na(region), "Other", region),
    continent = countrycode(countrycode, "iso3c", "continent"),
    continent = ifelse(is.na(continent), "Other", continent)
  )

cat(sprintf("Treated countries: %d\n", length(treated_countries)))
cat(sprintf("Never-treated countries: %d\n",
            n_distinct(data$countrycode) - length(treated_countries)))

# ------------------------------------------------------------------
# 2. Compute pre-treatment country-level means (pre-1990)
# ------------------------------------------------------------------
country_chars <- data %>%
  filter(year <= 1990) %>%
  group_by(countrycode, ever_treated, region, continent) %>%
  summarise(
    mean_loggdp   = mean(loggdp, na.rm = TRUE),
    mean_growth   = mean(gdp_diff, na.rm = TRUE),
    mean_hc       = mean(hc, na.rm = TRUE),
    mean_rkna     = mean(log(rkna + 1), na.rm = TRUE),
    mean_pop      = mean(log(pop + 1), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(complete.cases(.))

cat(sprintf("\nCountries with complete pre-1990 data: %d\n", nrow(country_chars)))
cat(sprintf("  Treated: %d\n", sum(country_chars$ever_treated)))
cat(sprintf("  Never-treated: %d\n", sum(!country_chars$ever_treated)))

# ------------------------------------------------------------------
# 3. Balance table: treated vs never-treated
# ------------------------------------------------------------------
balance <- country_chars %>%
  group_by(ever_treated) %>%
  summarise(
    across(starts_with("mean_"),
           list(avg = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}")
  ) %>%
  mutate(group = ifelse(ever_treated, "Treated", "Never-Treated"))

cat("\n========== BALANCE TABLE (Pre-1990 Means) ==========\n")
cat(sprintf("%-25s %15s %15s %12s\n",
            "Variable", "Treated", "Never-Treated", "Std. Diff"))
cat(paste(rep("-", 70), collapse = ""), "\n")

vars <- c("mean_loggdp", "mean_growth", "mean_hc", "mean_rkna", "mean_pop")
var_labels <- c("Log GDP per capita", "GDP Growth", "Human Capital (hc)",
                "Log Capital Stock", "Log Population")

std_diffs <- numeric(length(vars))
for (i in seq_along(vars)) {
  t_avg <- balance[[paste0(vars[i], "_avg")]][balance$ever_treated]
  t_sd  <- balance[[paste0(vars[i], "_sd")]][balance$ever_treated]
  c_avg <- balance[[paste0(vars[i], "_avg")]][!balance$ever_treated]
  c_sd  <- balance[[paste0(vars[i], "_sd")]][!balance$ever_treated]
  pooled_sd <- sqrt((t_sd^2 + c_sd^2) / 2)
  std_diffs[i] <- (t_avg - c_avg) / pooled_sd

  cat(sprintf("%-25s %7.3f (%4.2f) %7.3f (%4.2f)  %8.3f\n",
              var_labels[i], t_avg, t_sd, c_avg, c_sd, std_diffs[i]))
}

# ------------------------------------------------------------------
# 4. Region distribution
# ------------------------------------------------------------------
cat("\n========== REGION DISTRIBUTION ==========\n")
region_tab <- country_chars %>%
  count(ever_treated, continent) %>%
  pivot_wider(names_from = ever_treated, values_from = n, values_fill = 0) %>%
  rename(NeverTreated = `FALSE`, Treated = `TRUE`)
print(region_tab)

# ------------------------------------------------------------------
# 5. Distribution plots
# ------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

plot_df <- country_chars %>%
  mutate(group = ifelse(ever_treated, "Treated (40)", "Never-Treated (142)")) %>%
  select(countrycode, group, mean_loggdp, mean_growth, mean_hc, mean_pop) %>%
  pivot_longer(cols = starts_with("mean_"),
               names_to = "variable", values_to = "value") %>%
  mutate(variable = case_when(
    variable == "mean_loggdp" ~ "Log GDP per capita",
    variable == "mean_growth" ~ "GDP Growth",
    variable == "mean_hc"     ~ "Human Capital Index",
    variable == "mean_pop"    ~ "Log Population"
  ))

p_dist <- ggplot(plot_df, aes(x = value, fill = group, color = group)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  facet_wrap(~variable, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("Treated (40)" = "firebrick",
                                "Never-Treated (142)" = "steelblue")) +
  scale_color_manual(values = c("Treated (40)" = "firebrick",
                                 "Never-Treated (142)" = "steelblue")) +
  labs(x = NULL, y = "Density",
       title = "Pre-Treatment Covariate Distributions: Treated vs Never-Treated",
       subtitle = "Pre-1990 country means",
       fill = NULL, color = NULL) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "bottom"
  )

ggsave("figures/covariate_imbalance.png", p_dist,
       width = 10, height = 7, dpi = 300)
cat("\nSaved: figures/covariate_imbalance.png\n")

# ------------------------------------------------------------------
# 6. Standardized difference bar chart
# ------------------------------------------------------------------
smd_df <- data.frame(
  variable = var_labels,
  smd = abs(std_diffs)
) %>%
  mutate(variable = fct_reorder(variable, smd))

p_smd <- ggplot(smd_df, aes(x = smd, y = variable)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "firebrick") +
  annotate("text", x = 0.12, y = 0.5, label = "Good balance", size = 3, hjust = 0) +
  annotate("text", x = 0.27, y = 0.5, label = "Poor balance", size = 3, hjust = 0,
           color = "firebrick") +
  labs(x = "Absolute Standardized Mean Difference",
       y = NULL,
       title = "Covariate Imbalance: Treated vs Never-Treated",
       subtitle = "Pre-1990 means | Dashed lines: conventional thresholds") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40")
  )

ggsave("figures/smd_before_matching.png", p_smd,
       width = 8, height = 4.5, dpi = 300)
cat("Saved: figures/smd_before_matching.png\n")

# Save country characteristics for other scripts
saveRDS(country_chars, "country_chars_pre1990.rds")
saveRDS(data, "panel_data.rds")
cat("\nSaved: country_chars_pre1990.rds, panel_data.rds\n")
cat("Done.\n")
