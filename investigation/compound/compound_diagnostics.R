# compound_diagnostics.R
# Characterize compound vs isolated shocks among ever-treated countries
# Understand who gets compound shocks, when, and covariate differences

library(tidyverse)
library(haven)
library(countrycode)

# ------------------------------------------------------------------
# 1. Load and construct data (ever-treated only)
# ------------------------------------------------------------------
pwt  <- read_stata("../../raw_data/pwt_clean.dta")
tcs  <- read_stata("../../raw_data/ibtracs_clean.dta")
data <- pwt %>% left_join(tcs, by = c("year", "countrycode"))

knot_to_ms <- 0.514444

data <- data %>%
  mutate(
    maxwind  = replace_na(max_ann_wind_i_nob, 0) * knot_to_ms,
    energy   = replace_na(sum_ann_energy_i_nob, 0),
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

# Identify treated countries
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
    continent = ifelse(is.na(continent), "Other", continent),
    country_name = countrycode(countrycode, "iso3c", "country.name")
  )

# Restrict to ever-treated
d_treated <- data %>%
  filter(ever_treated) %>%
  arrange(countrycode, year)

cat(sprintf("Ever-treated countries: %d\n", length(treated_countries)))
cat(sprintf("Total obs (treated panel): %d\n", nrow(d_treated)))

# ------------------------------------------------------------------
# 2. Classify compound vs isolated shocks
# ------------------------------------------------------------------
# Definition: A shock in year t is "compound" if there was also a shock
# in t-1. This is the backward-looking definition — it asks whether the
# economy was already recovering from a recent shock when hit again.
d_treated <- d_treated %>%
  group_by(countrycode) %>%
  mutate(
    shock_lag1  = lag(maxwind_95, 1, default = 0L),
    shock_lag2  = lag(maxwind_95, 2, default = 0L),
    shock_lead1 = lead(maxwind_95, 1, default = 0L),
    # Compound: shock in t AND shock in t-1
    compound_strict = as.integer(maxwind_95 == 1 & shock_lag1 == 1),
    # Compound broad: shock in t AND shock in t-1 or t-2
    compound_broad  = as.integer(maxwind_95 == 1 & (shock_lag1 == 1 | shock_lag2 == 1)),
    # Isolated: shock in t, no shock in t-1 or t-2
    isolated_strict = as.integer(maxwind_95 == 1 & shock_lag1 == 0),
    isolated_broad  = as.integer(maxwind_95 == 1 & shock_lag1 == 0 & shock_lag2 == 0),
    # Adjacent: shock in t AND (shock in t-1 OR t+1)
    adjacent = as.integer(maxwind_95 == 1 & (shock_lag1 == 1 | shock_lead1 == 1))
  ) %>%
  ungroup()

cat("\n========== SHOCK CLASSIFICATION ==========\n")
cat(sprintf("Total shock-years: %d\n", sum(d_treated$maxwind_95)))
cat(sprintf("Compound (strict, t-1): %d (%.1f%%)\n",
            sum(d_treated$compound_strict),
            mean(d_treated$compound_strict[d_treated$maxwind_95 == 1]) * 100))
cat(sprintf("Compound (broad, t-1 or t-2): %d (%.1f%%)\n",
            sum(d_treated$compound_broad),
            mean(d_treated$compound_broad[d_treated$maxwind_95 == 1]) * 100))
cat(sprintf("Isolated (strict, no t-1): %d (%.1f%%)\n",
            sum(d_treated$isolated_strict),
            mean(d_treated$isolated_strict[d_treated$maxwind_95 == 1]) * 100))
cat(sprintf("Isolated (broad, no t-1 or t-2): %d (%.1f%%)\n",
            sum(d_treated$isolated_broad),
            mean(d_treated$isolated_broad[d_treated$maxwind_95 == 1]) * 100))

# ------------------------------------------------------------------
# 3. Country-level breakdown
# ------------------------------------------------------------------
country_shock <- d_treated %>%
  filter(maxwind_95 == 1) %>%
  group_by(countrycode, country_name, continent) %>%
  summarise(
    n_shocks       = n(),
    n_compound     = sum(compound_strict),
    n_isolated     = sum(isolated_strict),
    pct_compound   = mean(compound_strict) * 100,
    first_shock_yr = min(year),
    last_shock_yr  = max(year),
    .groups = "drop"
  ) %>%
  arrange(desc(n_shocks))

cat("\n========== COUNTRY-LEVEL SHOCK PATTERNS ==========\n")
cat(sprintf("%-5s %-20s %-10s %6s %6s %6s %6s\n",
            "ISO", "Country", "Continent", "Total", "Comp", "Isol", "%Comp"))
cat(paste(rep("-", 70), collapse = ""), "\n")
for (i in 1:nrow(country_shock)) {
  r <- country_shock[i, ]
  cat(sprintf("%-5s %-20s %-10s %6d %6d %6d %5.0f%%\n",
              r$countrycode, substr(r$country_name, 1, 20), r$continent,
              r$n_shocks, r$n_compound, r$n_isolated, r$pct_compound))
}

# ------------------------------------------------------------------
# 4. Countries with BOTH compound and isolated shocks (within-variation)
# ------------------------------------------------------------------
both_types <- country_shock %>%
  filter(n_compound > 0 & n_isolated > 0)

cat(sprintf("\n%d countries have BOTH compound and isolated shocks (within-country variation)\n",
            nrow(both_types)))
print(both_types %>% select(countrycode, country_name, n_shocks, n_compound, n_isolated))

# ------------------------------------------------------------------
# 5. Covariate comparison: compound-heavy vs isolated-heavy countries
# ------------------------------------------------------------------
# Classify countries by their dominant shock type
country_shock <- country_shock %>%
  mutate(
    shock_type = case_when(
      pct_compound >= 50 ~ "Compound-heavy",
      TRUE ~ "Isolated-heavy"
    )
  )

# Pre-treatment characteristics
pre_chars <- d_treated %>%
  filter(year <= 1990) %>%
  group_by(countrycode) %>%
  summarise(
    mean_loggdp = mean(loggdp, na.rm = TRUE),
    mean_growth = mean(gdp_diff, na.rm = TRUE),
    mean_hc     = mean(hc, na.rm = TRUE),
    mean_pop    = mean(log(pop + 1), na.rm = TRUE),
    .groups = "drop"
  )

country_compare <- country_shock %>%
  left_join(pre_chars, by = "countrycode") %>%
  filter(complete.cases(.))

cat("\n========== COMPOUND-HEAVY vs ISOLATED-HEAVY COUNTRIES ==========\n")
compare_tab <- country_compare %>%
  group_by(shock_type) %>%
  summarise(
    n         = n(),
    loggdp    = mean(mean_loggdp, na.rm = TRUE),
    growth    = mean(mean_growth, na.rm = TRUE),
    hc        = mean(mean_hc, na.rm = TRUE),
    pop       = mean(mean_pop, na.rm = TRUE),
    n_shocks  = mean(n_shocks),
    .groups   = "drop"
  )
print(compare_tab)

# ------------------------------------------------------------------
# 6. Temporal patterns
# ------------------------------------------------------------------
shock_by_year <- d_treated %>%
  filter(maxwind_95 == 1) %>%
  group_by(year) %>%
  summarise(
    n_shocks   = n(),
    n_compound = sum(compound_strict),
    n_isolated = sum(isolated_strict),
    pct_compound = mean(compound_strict) * 100
  )

cat("\n========== TEMPORAL PATTERNS ==========\n")
cat("Pre-1990:\n")
pre <- shock_by_year %>% filter(year <= 1990)
cat(sprintf("  Total shocks: %d, Compound: %d (%.1f%%), Isolated: %d\n",
            sum(pre$n_shocks), sum(pre$n_compound),
            sum(pre$n_compound) / sum(pre$n_shocks) * 100,
            sum(pre$n_isolated)))
cat("Post-1990:\n")
post <- shock_by_year %>% filter(year > 1990)
cat(sprintf("  Total shocks: %d, Compound: %d (%.1f%%), Isolated: %d\n",
            sum(post$n_shocks), sum(post$n_compound),
            sum(post$n_compound) / sum(post$n_shocks) * 100,
            sum(post$n_isolated)))

# ------------------------------------------------------------------
# 7. Figures
# ------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# 7a. Shock timeline for selected countries
selected <- c("CHN", "PHL", "JPN", "USA", "MEX", "NIC", "BHS", "DOM",
              "KOR", "VNM", "CUB", "AUS")
timeline_df <- d_treated %>%
  filter(countrycode %in% selected, maxwind_95 == 1) %>%
  mutate(
    type = case_when(
      compound_strict == 1 ~ "Compound",
      TRUE ~ "Isolated"
    ),
    country_name = countrycode(countrycode, "iso3c", "country.name"),
    country_name = fct_reorder(country_name, -as.numeric(factor(countrycode)))
  )

p_timeline <- ggplot(timeline_df, aes(x = year, y = country_name,
                                       color = type, shape = type)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Compound" = "firebrick", "Isolated" = "steelblue")) +
  scale_shape_manual(values = c("Compound" = 17, "Isolated" = 16)) +
  geom_vline(xintercept = 1990, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(x = "Year", y = NULL,
       title = "Shock Timeline: Compound vs Isolated",
       subtitle = "Compound = shock in year t with shock also in t-1",
       color = NULL, shape = NULL) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "bottom"
  )

ggsave("figures/shock_timeline.png", p_timeline,
       width = 12, height = 6, dpi = 300)
cat("\nSaved: figures/shock_timeline.png\n")

# 7b. Compound share over time
p_temporal <- ggplot(shock_by_year, aes(x = year)) +
  geom_col(aes(y = n_shocks), fill = "grey80", width = 0.8) +
  geom_col(aes(y = n_compound), fill = "firebrick", alpha = 0.7, width = 0.8) +
  geom_vline(xintercept = 1990, linetype = "dashed", color = "grey30") +
  labs(x = "Year", y = "Number of shock events",
       title = "Shock Events by Year",
       subtitle = "Red = compound (shock in t-1 also); Grey = all shocks") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40")
  )

ggsave("figures/shock_temporal.png", p_temporal,
       width = 10, height = 5, dpi = 300)
cat("Saved: figures/shock_temporal.png\n")

# 7c. Gap distribution
gap_df <- d_treated %>%
  filter(maxwind_95 == 1) %>%
  group_by(countrycode) %>%
  mutate(gap = year - lag(year)) %>%
  filter(!is.na(gap))

p_gaps <- ggplot(gap_df, aes(x = gap)) +
  geom_histogram(binwidth = 1, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "firebrick", linewidth = 0.5) +
  annotate("text", x = 2.5, y = max(table(gap_df$gap)) * 0.9,
           label = "gap = 1: compound", color = "firebrick", size = 3.5, hjust = 0) +
  labs(x = "Years between consecutive shocks (same country)",
       y = "Count",
       title = "Distribution of Inter-Shock Gaps",
       subtitle = "Among ever-treated countries with multiple shocks") +
  scale_x_continuous(breaks = seq(0, 25, 5)) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40")
  )

ggsave("figures/gap_distribution.png", p_gaps,
       width = 8, height = 5, dpi = 300)
cat("Saved: figures/gap_distribution.png\n")

# 7d. Country characteristics by shock type
char_plot_df <- country_compare %>%
  select(countrycode, shock_type, mean_loggdp, mean_growth, mean_hc, mean_pop) %>%
  pivot_longer(starts_with("mean_"), names_to = "variable", values_to = "value") %>%
  mutate(variable = case_when(
    variable == "mean_loggdp" ~ "Log GDP per capita",
    variable == "mean_growth" ~ "GDP Growth",
    variable == "mean_hc"     ~ "Human Capital",
    variable == "mean_pop"    ~ "Log Population"
  ))

p_chars <- ggplot(char_plot_df, aes(x = value, fill = shock_type, color = shock_type)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  facet_wrap(~variable, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("Compound-heavy" = "firebrick",
                                "Isolated-heavy" = "steelblue")) +
  scale_color_manual(values = c("Compound-heavy" = "firebrick",
                                 "Isolated-heavy" = "steelblue")) +
  labs(x = NULL, y = "Density",
       title = "Pre-1990 Characteristics by Shock Pattern",
       subtitle = "Compound-heavy (>=50% compound) vs Isolated-heavy countries",
       fill = NULL, color = NULL) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "bottom"
  )

ggsave("figures/chars_by_shock_type.png", p_chars,
       width = 9, height = 6, dpi = 300)
cat("Saved: figures/chars_by_shock_type.png\n")

# Save processed data
saveRDS(d_treated, "treated_panel.rds")
saveRDS(country_shock, "country_shock_summary.rds")
cat("\nSaved: treated_panel.rds, country_shock_summary.rds\n")
cat("Done.\n")
