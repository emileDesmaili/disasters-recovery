# investigation_counterfactual.R
# Decompose the TWFE counterfactual mechanism behind the post-1990 divergence.
#
# Key insight: In a TWFE LP, year fixed effects absorb the average GDP path
# across ALL countries. When 142/182 countries never experience a p95 cyclone
# shock, those countries dominate the year FE. If their GDP trends shifted
# around 1990 (end of Cold War, globalization, China's rise), this shifts the
# implicit counterfactual that treated countries are compared against.
#
# Analyses:
#   1. Year FE decomposition: full vs ever-treated samples
#   2. Counterfactual GDP paths: never-treated vs ever-treated
#   3. Placebo test on never-treated countries
#   4. Differential trends test
#   5. Synthetic counterfactual with continent x year FE

library(tidyverse)
library(fixest)
library(haven)
library(countrycode)

source("emileRegs.R")

setFixest_notes(FALSE)

cat("============================================================\n")
cat("  COUNTERFACTUAL DECOMPOSITION INVESTIGATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------------
# 1. Load and construct data (identically to other investigation scripts)
# ------------------------------------------------------------------
pwt  <- read_stata("raw_data/pwt_clean.dta")
tcs  <- read_stata("raw_data/ibtracs_clean.dta")
data <- pwt %>% left_join(tcs, by = c("year", "countrycode"))

knot_to_ms <- 0.514444

data <- data %>%
  mutate(
    maxwind  = replace_na(max_ann_wind_i_nob, 0) * knot_to_ms,
    energy   = replace_na(sum_ann_energy_i_nob, 0),
    nlands   = replace_na(sum_a_lands_nob, 0),
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

data <- data %>%
  mutate(maxwind_sd = maxwind / sd(maxwind[maxwind > 0], na.rm = TRUE))

# ------------------------------------------------------------------
# LP settings (matching lp_gdp_shocks.R)
# ------------------------------------------------------------------
outcome  <- "loggdp"
horizon  <- 10
fe       <- "countrycode[year] + countrycode[year2] + year"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

make_controls <- function(shock_var) {
  paste0("l(gdp_diff,1:2) + l(", shock_var, ",1:2)")
}

controls <- make_controls("maxwind_95")

# ------------------------------------------------------------------
# Identify treatment groups
# ------------------------------------------------------------------
ever_treated_countries <- data %>%
  filter(maxwind_95 == 1) %>%
  pull(countrycode) %>%
  unique() %>%
  sort()

all_countries <- unique(data$countrycode)
never_treated_countries <- sort(setdiff(all_countries, ever_treated_countries))

data <- data %>%
  mutate(
    ever_treated = as.integer(countrycode %in% ever_treated_countries)
  )

cat(sprintf("Total countries: %d\n", length(all_countries)))
cat(sprintf("Ever-treated:    %d\n", length(ever_treated_countries)))
cat(sprintf("Never-treated:   %d\n", length(never_treated_countries)))

# Add continent variable
data <- data %>%
  mutate(
    continent = countrycode(countrycode, "iso3c", "continent"),
    region    = countrycode(countrycode, "iso3c", "region")
  )

# Check continent coverage
cat("\n--- Continent distribution ---\n")
data %>%
  distinct(countrycode, continent, ever_treated) %>%
  group_by(continent) %>%
  summarise(
    n_countries   = n(),
    n_treated     = sum(ever_treated),
    n_never_treat = sum(1 - ever_treated),
    .groups = "drop"
  ) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# Subsets
data_full         <- data
data_ever_treated <- data %>% filter(countrycode %in% ever_treated_countries)

# Plotting theme
irf_theme <- theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40")
  )


# ==================================================================
# ANALYSIS 1: Year FE Decomposition
# ==================================================================
cat("\n============================================================\n")
cat("  ANALYSIS 1: Year Fixed Effect Decomposition\n")
cat("============================================================\n\n")

# Run LP at h=5 on full sample and ever-treated only, extracting year FE
h <- 5

fml_str <- paste0(
  "f(loggdp, ", h, ") - l(loggdp, 1) ~ ",
  "i(period, maxwind_95) + l(gdp_diff, 1:2) + l(maxwind_95, 1:2)",
  " | countrycode[year] + countrycode[year2] + year"
)

cat("Estimating LP at h=5 on full sample...\n")
mod_full <- feols(
  as.formula(fml_str),
  data = data_full,
  panel.id = panel_id,
  vcov = vcov_fm
)

cat("Estimating LP at h=5 on ever-treated sample...\n")
mod_ever <- feols(
  as.formula(fml_str),
  data = data_ever_treated,
  panel.id = panel_id,
  vcov = vcov_fm
)

# Extract year fixed effects
fe_full <- fixef(mod_full)$year
fe_ever <- fixef(mod_ever)$year

fe_df <- bind_rows(
  data.frame(
    year   = as.integer(names(fe_full)),
    year_fe = as.numeric(fe_full),
    sample  = "Full Sample (182 countries)"
  ),
  data.frame(
    year   = as.integer(names(fe_ever)),
    year_fe = as.numeric(fe_ever),
    sample  = "Ever-Treated Only (40 countries)"
  )
)

# Normalize: subtract the mean of each series so they're centered at 0
fe_df <- fe_df %>%
  group_by(sample) %>%
  mutate(year_fe_centered = year_fe - mean(year_fe, na.rm = TRUE)) %>%
  ungroup()

# Print comparison
cat("\n--- Year FE comparison at key years ---\n")
fe_wide <- fe_df %>%
  select(year, sample, year_fe_centered) %>%
  pivot_wider(names_from = sample, values_from = year_fe_centered)
print(as.data.frame(fe_wide), row.names = FALSE, digits = 4)

# Plot year FE
p_yfe <- ggplot(fe_df, aes(x = year, y = year_fe_centered, color = sample)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = 1990, linetype = "dashed", color = "grey40") +
  annotate("text", x = 1990.5, y = max(fe_df$year_fe_centered, na.rm = TRUE),
           label = "1990", hjust = 0, vjust = 1, size = 3.5, color = "grey40") +
  scale_color_manual(values = c(
    "Full Sample (182 countries)"        = "steelblue",
    "Ever-Treated Only (40 countries)"   = "firebrick"
  )) +
  labs(
    x     = "Year",
    y     = "Year FE (centered)",
    color = "Sample",
    title = "Year Fixed Effects from LP at h=5",
    subtitle = "Full sample vs ever-treated only | maxwind_95 with period interaction"
  ) +
  irf_theme +
  theme(legend.position = "top")

ggsave("investigation/figures/year_fe_decomposition.png", p_yfe,
       width = 10, height = 6, dpi = 300)
cat("Saved: investigation/figures/year_fe_decomposition.png\n")

# Compute the correlation and divergence
fe_merged <- fe_df %>%
  select(year, sample, year_fe_centered) %>%
  pivot_wider(names_from = sample, values_from = year_fe_centered) %>%
  drop_na()

cor_val <- cor(fe_merged[[2]], fe_merged[[3]])
cat(sprintf("\nCorrelation between full-sample and ever-treated year FE: %.4f\n", cor_val))

# Pre vs post 1990 divergence
fe_merged <- fe_merged %>%
  mutate(
    diff = .[[2]] - .[[3]],
    period = ifelse(year <= 1990, "Pre-1990", "Post-1990")
  )

cat("\n--- Mean year FE difference (Full - EverTreated) by period ---\n")
fe_merged %>%
  group_by(period) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    sd_diff   = sd(diff, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  as.data.frame() %>%
  print(row.names = FALSE, digits = 4)


# ==================================================================
# ANALYSIS 2: Counterfactual GDP Paths
# ==================================================================
cat("\n\n============================================================\n")
cat("  ANALYSIS 2: Counterfactual GDP Paths\n")
cat("============================================================\n\n")

# Compute annual average GDP growth for never-treated vs ever-treated
# For ever-treated, compute separately for shock vs non-shock years
gdp_paths <- data %>%
  filter(!is.na(gdp_diff)) %>%
  mutate(
    group = case_when(
      ever_treated == 0 ~ "Never-treated",
      ever_treated == 1 & maxwind_95 == 0 ~ "Ever-treated (non-shock years)",
      ever_treated == 1 & maxwind_95 == 1 ~ "Ever-treated (shock years)"
    )
  ) %>%
  group_by(year, group) %>%
  summarise(
    mean_growth  = mean(gdp_diff, na.rm = TRUE),
    median_growth = median(gdp_diff, na.rm = TRUE),
    n_obs        = n(),
    .groups      = "drop"
  )

# Also compute just never vs ever treated (all years)
gdp_paths_simple <- data %>%
  filter(!is.na(gdp_diff)) %>%
  mutate(group = ifelse(ever_treated == 1, "Ever-treated", "Never-treated")) %>%
  group_by(year, group) %>%
  summarise(
    mean_growth = mean(gdp_diff, na.rm = TRUE),
    n_obs       = n(),
    .groups     = "drop"
  )

# Summary statistics by period
cat("--- Average GDP growth by group and period ---\n")
gdp_summary <- data %>%
  filter(!is.na(gdp_diff)) %>%
  mutate(group = ifelse(ever_treated == 1, "Ever-treated", "Never-treated")) %>%
  group_by(group, period) %>%
  summarise(
    mean_growth   = mean(gdp_diff, na.rm = TRUE),
    median_growth = median(gdp_diff, na.rm = TRUE),
    sd_growth     = sd(gdp_diff, na.rm = TRUE),
    n_obs         = n(),
    .groups       = "drop"
  )
print(as.data.frame(gdp_summary), row.names = FALSE, digits = 4)

# Difference-in-difference of growth rates
growth_did <- gdp_summary %>%
  select(group, period, mean_growth) %>%
  pivot_wider(names_from = period, values_from = mean_growth) %>%
  mutate(change = `Post-1990` - `Pre-1990`)
cat("\n--- GDP growth rate change (Post - Pre 1990) ---\n")
print(as.data.frame(growth_did), row.names = FALSE, digits = 4)
cat(sprintf("\nDiD of growth rates: %.4f\n",
            growth_did$change[growth_did$group == "Ever-treated"] -
            growth_did$change[growth_did$group == "Never-treated"]))

# Plot: annual mean GDP growth by group with loess smoother
p_gdp <- ggplot(gdp_paths_simple, aes(x = year, y = mean_growth, color = group)) +
  geom_line(alpha = 0.4, linewidth = 0.5) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, linewidth = 1.3) +
  geom_vline(xintercept = 1990, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  annotate("text", x = 1990.5, y = max(gdp_paths_simple$mean_growth, na.rm = TRUE),
           label = "1990", hjust = 0, vjust = 1, size = 3.5, color = "grey40") +
  scale_color_manual(values = c(
    "Never-treated"  = "steelblue",
    "Ever-treated"   = "firebrick"
  )) +
  labs(
    x     = "Year",
    y     = "Mean GDP growth rate",
    color = "Group",
    title = "Average GDP Growth Paths: Never-Treated vs Ever-Treated",
    subtitle = paste0("Never-treated: ", length(never_treated_countries),
                      " countries | Ever-treated: ", length(ever_treated_countries),
                      " countries\nSmooth trend (loess) overlaid on annual means")
  ) +
  irf_theme +
  theme(legend.position = "top")

ggsave("investigation/figures/counterfactual_gdp_paths.png", p_gdp,
       width = 10, height = 6, dpi = 300)
cat("Saved: investigation/figures/counterfactual_gdp_paths.png\n")


# ==================================================================
# ANALYSIS 3: Placebo Test on Never-Treated Countries
# ==================================================================
cat("\n\n============================================================\n")
cat("  ANALYSIS 3: Placebo Test on Never-Treated Countries\n")
cat("============================================================\n\n")

# Assign random placebo shocks at the same frequency as real shocks
data_never <- data %>% filter(countrycode %in% never_treated_countries)

# Real shock frequency
real_shock_freq <- mean(data$maxwind_95, na.rm = TRUE)
cat(sprintf("Real shock frequency: %.4f (%.2f%%)\n", real_shock_freq, 100 * real_shock_freq))
cat(sprintf("N observations in never-treated sample: %d\n", nrow(data_never)))

n_sims <- 500
set.seed(42)

cat(sprintf("Running %d placebo simulations...\n", n_sims))

# For each simulation, assign placebo shocks, run period-split LP at h=5,
# and store the pre/post coefficient gap
placebo_results <- vector("list", n_sims)

for (sim in 1:n_sims) {
  if (sim %% 50 == 0) cat(sprintf("  Simulation %d/%d\n", sim, n_sims))

  # Assign placebo shocks
  data_never$placebo_shock <- rbinom(nrow(data_never), 1, real_shock_freq)

  # Run LP at h=5 with period interaction
  fml_placebo <- as.formula(paste0(
    "f(loggdp, 5) - l(loggdp, 1) ~ ",
    "i(period, placebo_shock) + l(gdp_diff, 1:2) + l(placebo_shock, 1:2)",
    " | countrycode[year] + countrycode[year2] + year"
  ))

  tryCatch({
    mod_placebo <- feols(
      fml_placebo,
      data = data_never,
      panel.id = panel_id,
      vcov = vcov_fm
    )

    coefs <- coef(mod_placebo)
    # Extract Pre-1990 and Post-1990 coefficients
    pre_name  <- grep("Pre.*placebo", names(coefs), value = TRUE)
    post_name <- grep("Post.*placebo", names(coefs), value = TRUE)

    if (length(pre_name) == 1 && length(post_name) == 1) {
      placebo_results[[sim]] <- data.frame(
        sim       = sim,
        pre_coef  = coefs[pre_name],
        post_coef = coefs[post_name],
        gap       = coefs[post_name] - coefs[pre_name]
      )
    }
  }, error = function(e) {
    # Skip failed simulations
  })
}

placebo_df <- bind_rows(placebo_results)
cat(sprintf("Successful simulations: %d/%d\n", nrow(placebo_df), n_sims))

# Now get the actual gap from the real data at h=5
cat("\nEstimating actual LP at h=5 on full sample for comparison...\n")
irf_actual <- lp_panel_inter(
  data = data, outcome = outcome, main_var = "maxwind_95",
  interact_var = "period", controls = controls,
  horizon = 5, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
)

# Clean category names
irf_actual <- irf_actual %>%
  mutate(
    category = str_replace_all(category, "_", "-"),
    category = case_when(
      str_detect(category, "Pre")  ~ "Pre-1990",
      str_detect(category, "Post") ~ "Post-1990",
      TRUE ~ category
    )
  )

actual_pre  <- irf_actual %>% filter(horizon == 5, category == "Pre-1990") %>% pull(irf_mean)
actual_post <- irf_actual %>% filter(horizon == 5, category == "Post-1990") %>% pull(irf_mean)
actual_gap  <- actual_post - actual_pre

cat(sprintf("\nActual coefficient gap at h=5 (Post - Pre): %.4f\n", actual_gap))
cat(sprintf("Placebo gap mean:   %.4f\n", mean(placebo_df$gap, na.rm = TRUE)))
cat(sprintf("Placebo gap sd:     %.4f\n", sd(placebo_df$gap, na.rm = TRUE)))
cat(sprintf("Placebo gap median: %.4f\n", median(placebo_df$gap, na.rm = TRUE)))

# P-value: fraction of placebo gaps more extreme than actual
p_val <- mean(abs(placebo_df$gap) >= abs(actual_gap), na.rm = TRUE)
cat(sprintf("P-value (two-sided): %.4f\n", p_val))

# Plot distribution of placebo gaps
p_placebo <- ggplot(placebo_df, aes(x = gap)) +
  geom_histogram(bins = 50, fill = "grey70", color = "grey50", alpha = 0.8) +
  geom_vline(xintercept = actual_gap, color = "firebrick", linewidth = 1.2, linetype = "solid") +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.5, linetype = "dashed") +
  annotate("text",
           x = actual_gap, y = Inf,
           label = sprintf("Actual gap = %.3f", actual_gap),
           hjust = -0.1, vjust = 2, color = "firebrick", size = 4, fontface = "bold") +
  annotate("text",
           x = max(placebo_df$gap, na.rm = TRUE) * 0.7, y = Inf,
           label = sprintf("p-value = %.3f\n%d simulations", p_val, nrow(placebo_df)),
           hjust = 0, vjust = 3.5, size = 3.5, color = "grey30") +
  labs(
    x     = "Placebo post-pre gap at h=5",
    y     = "Count",
    title = "Placebo Test: Random Shocks on Never-Treated Countries",
    subtitle = paste0("Distribution of (Post-1990 - Pre-1990) coefficient gap from ",
                      n_sims, " simulations\nRed line = actual gap from real data")
  ) +
  irf_theme

ggsave("investigation/figures/placebo_never_treated.png", p_placebo,
       width = 10, height = 6, dpi = 300)
cat("Saved: investigation/figures/placebo_never_treated.png\n")


# ==================================================================
# ANALYSIS 4: Differential Trends Test
# ==================================================================
cat("\n\n============================================================\n")
cat("  ANALYSIS 4: Differential Trends Test\n")
cat("============================================================\n\n")

cat("Testing whether ever-treated and never-treated countries have\n")
cat("differential GDP growth trend breaks around 1990.\n\n")

# Regression: gdp_diff ~ period * ever_treated | countrycode
# This directly tests whether the two groups have different trend breaks

data_trends <- data %>%
  filter(!is.na(gdp_diff)) %>%
  mutate(
    post1990     = as.integer(year > 1990),
    ever_treated = as.integer(countrycode %in% ever_treated_countries)
  )

# Model 1: Simple DiD with country FE
cat("--- Model 1: GDP growth ~ Post1990 * EverTreated | countrycode ---\n\n")
mod_trends1 <- feols(
  gdp_diff ~ post1990 * ever_treated | countrycode,
  data = data_trends,
  vcov = DK ~ year
)
print(summary(mod_trends1))

# Model 2: With year trends
cat("\n--- Model 2: GDP growth ~ Post1990 * EverTreated | countrycode[year] ---\n\n")
mod_trends2 <- feols(
  gdp_diff ~ post1990 * ever_treated | countrycode[year],
  data = data_trends,
  vcov = DK ~ year
)
print(summary(mod_trends2))

# Model 3: With year FE (absorbed, testing interaction)
cat("\n--- Model 3: GDP growth ~ Post1990:EverTreated | countrycode + year ---\n\n")
mod_trends3 <- feols(
  gdp_diff ~ post1990:ever_treated | countrycode + year,
  data = data_trends,
  vcov = DK ~ year
)
print(summary(mod_trends3))

# Model 4: Log GDP level regression (same as LP at h=0 basically)
cat("\n--- Model 4: Log GDP ~ Post1990:EverTreated | countrycode + year ---\n\n")
mod_trends4 <- feols(
  loggdp ~ post1990:ever_treated | countrycode + year,
  data = data_trends %>% filter(!is.na(loggdp)),
  vcov = DK ~ year
)
print(summary(mod_trends4))

cat("\n--- Interpretation ---\n")
cat("If the post1990:ever_treated interaction is significant, it means\n")
cat("the two groups experienced different GDP trajectory shifts around 1990.\n")
cat("This would mean the year FE in the TWFE LP are picking up different\n")
cat("counterfactual paths for treated vs untreated countries.\n\n")

# Extract the key coefficient
did_coef <- coef(mod_trends3)["post1990:ever_treated"]
did_se   <- sqrt(vcov(mod_trends3)["post1990:ever_treated", "post1990:ever_treated"])
did_t    <- did_coef / did_se
cat(sprintf("DiD coefficient (post1990 x ever_treated, with country + year FE):\n"))
cat(sprintf("  Estimate:  %.5f\n", did_coef))
cat(sprintf("  SE:        %.5f\n", did_se))
cat(sprintf("  t-stat:    %.3f\n", did_t))
cat(sprintf("  p-value:   %.4f\n", 2 * pt(abs(did_t), df = mod_trends3$nobs - length(coef(mod_trends3)), lower.tail = FALSE)))


# ==================================================================
# ANALYSIS 5: Synthetic Counterfactual (Continent x Year FE)
# ==================================================================
cat("\n\n============================================================\n")
cat("  ANALYSIS 5: Continent x Year FE (Flexible Counterfactual)\n")
cat("============================================================\n\n")

# Create continent x year interaction for FE
# Handle missing continent codes
data <- data %>%
  mutate(
    continent = ifelse(is.na(continent), "Other", continent),
    continent_year = paste0(continent, "_", year)
  )

# Check how many countries per continent
cat("--- Countries per continent (with treatment status) ---\n")
data %>%
  distinct(countrycode, continent, ever_treated) %>%
  group_by(continent) %>%
  summarise(
    n_total   = n(),
    n_treated = sum(ever_treated),
    .groups   = "drop"
  ) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# A) Standard year FE (baseline)
cat("\n\nA) Running period-split LP with standard year FE...\n")
irf_standard <- lp_panel_inter(
  data = data, outcome = outcome, main_var = "maxwind_95",
  interact_var = "period", controls = controls,
  horizon = horizon, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
) %>%
  mutate(
    category = str_replace_all(category, "_", "-"),
    category = case_when(
      str_detect(category, "Pre")  ~ "Pre-1990",
      str_detect(category, "Post") ~ "Post-1990",
      TRUE ~ category
    ),
    spec = "Year FE"
  )

# B) Continent x Year FE
cat("B) Running period-split LP with continent x year FE...\n")
fe_continent <- "countrycode[year] + countrycode[year2] + continent^year"

irf_continent <- tryCatch({
  lp_panel_inter(
    data = data, outcome = outcome, main_var = "maxwind_95",
    interact_var = "period", controls = controls,
    horizon = horizon,
    fe = fe_continent,
    panel_id = panel_id, vcov_formula = vcov_fm
  ) %>%
    mutate(
      category = str_replace_all(category, "_", "-"),
      category = case_when(
        str_detect(category, "Pre")  ~ "Pre-1990",
        str_detect(category, "Post") ~ "Post-1990",
        TRUE ~ category
      ),
      spec = "Continent x Year FE"
    )
}, error = function(e) {
  cat(sprintf("  Error with continent x year FE: %s\n", e$message))
  cat("  Trying manual approach...\n")

  # Manual LP loop
  irf_list <- list()
  for (hh in 0:horizon) {
    fml_str <- paste0(
      "f(loggdp, ", hh, ") - l(loggdp, 1) ~ ",
      "i(period, maxwind_95) + l(gdp_diff, 1:2) + l(maxwind_95, 1:2)",
      " | countrycode[year] + countrycode[year2] + continent^year"
    )
    mod <- feols(as.formula(fml_str), data = data,
                 panel.id = panel_id, vcov = vcov_fm)

    coef_df <- data.frame(
      horizon  = hh,
      term     = names(coef(mod)),
      irf_mean = coef(mod),
      se       = sqrt(diag(vcov(mod))),
      stringsAsFactors = FALSE
    ) %>%
      filter(str_detect(term, "^period::")) %>%
      mutate(
        category = str_extract(term, "(?<=::).*?(?=:)"),
        main_var = "maxwind_95",
        irf_down = irf_mean - 1.96 * se,
        irf_up   = irf_mean + 1.96 * se
      )

    irf_list[[hh + 1]] <- coef_df
  }
  result <- bind_rows(irf_list) %>%
    mutate(
      category = str_replace_all(category, "_", "-"),
      category = case_when(
        str_detect(category, "Pre")  ~ "Pre-1990",
        str_detect(category, "Post") ~ "Post-1990",
        TRUE ~ category
      ),
      spec = "Continent x Year FE"
    )
  return(result)
})

# C) Region x Year FE (more granular)
cat("C) Running period-split LP with region x year FE...\n")

# Check region coverage
data <- data %>%
  mutate(region = ifelse(is.na(region), "Other", region))

cat("--- Regions in data ---\n")
data %>%
  distinct(countrycode, region, ever_treated) %>%
  group_by(region) %>%
  summarise(
    n_total   = n(),
    n_treated = sum(ever_treated),
    .groups   = "drop"
  ) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

fe_region <- "countrycode[year] + countrycode[year2] + region^year"

irf_region <- tryCatch({
  lp_panel_inter(
    data = data, outcome = outcome, main_var = "maxwind_95",
    interact_var = "period", controls = controls,
    horizon = horizon,
    fe = fe_region,
    panel_id = panel_id, vcov_formula = vcov_fm
  ) %>%
    mutate(
      category = str_replace_all(category, "_", "-"),
      category = case_when(
        str_detect(category, "Pre")  ~ "Pre-1990",
        str_detect(category, "Post") ~ "Post-1990",
        TRUE ~ category
      ),
      spec = "Region x Year FE"
    )
}, error = function(e) {
  cat(sprintf("  Error with region x year FE: %s\n", e$message))
  cat("  Trying manual approach...\n")

  irf_list <- list()
  for (hh in 0:horizon) {
    fml_str <- paste0(
      "f(loggdp, ", hh, ") - l(loggdp, 1) ~ ",
      "i(period, maxwind_95) + l(gdp_diff, 1:2) + l(maxwind_95, 1:2)",
      " | countrycode[year] + countrycode[year2] + region^year"
    )

    tryCatch({
      mod <- feols(as.formula(fml_str), data = data,
                   panel.id = panel_id, vcov = vcov_fm)

      coef_df <- data.frame(
        horizon  = hh,
        term     = names(coef(mod)),
        irf_mean = coef(mod),
        se       = sqrt(diag(vcov(mod))),
        stringsAsFactors = FALSE
      ) %>%
        filter(str_detect(term, "^period::")) %>%
        mutate(
          category = str_extract(term, "(?<=::).*?(?=:)"),
          main_var = "maxwind_95",
          irf_down = irf_mean - 1.96 * se,
          irf_up   = irf_mean + 1.96 * se
        )

      irf_list[[hh + 1]] <- coef_df
    }, error = function(e2) {
      cat(sprintf("    Skipping h=%d: %s\n", hh, e2$message))
    })
  }

  if (length(irf_list) > 0) {
    result <- bind_rows(irf_list) %>%
      mutate(
        category = str_replace_all(category, "_", "-"),
        category = case_when(
          str_detect(category, "Pre")  ~ "Pre-1990",
          str_detect(category, "Post") ~ "Post-1990",
          TRUE ~ category
        ),
        spec = "Region x Year FE"
      )
    return(result)
  } else {
    return(NULL)
  }
})

# Combine and plot
all_specs <- bind_rows(
  irf_standard,
  irf_continent,
  irf_region
) %>%
  filter(!is.na(spec))

cat("\n--- Coefficient comparison at h=5 ---\n")
all_specs %>%
  filter(horizon == 5) %>%
  select(spec, category, irf_mean, se, irf_down, irf_up) %>%
  arrange(category, spec) %>%
  as.data.frame() %>%
  print(row.names = FALSE, digits = 4)

cat("\n--- Coefficient comparison at h=10 ---\n")
all_specs %>%
  filter(horizon == 10) %>%
  select(spec, category, irf_mean, se, irf_down, irf_up) %>%
  arrange(category, spec) %>%
  as.data.frame() %>%
  print(row.names = FALSE, digits = 4)

# Pre-post gap comparison
cat("\n--- Pre-Post gap by specification ---\n")
gap_table <- all_specs %>%
  filter(horizon %in% c(5, 10)) %>%
  select(spec, category, horizon, irf_mean) %>%
  pivot_wider(names_from = category, values_from = irf_mean) %>%
  mutate(gap = `Post-1990` - `Pre-1990`)

print(as.data.frame(gap_table), row.names = FALSE, digits = 4)

# Plot comparison across specifications
specs_present <- unique(all_specs$spec)
n_specs <- length(specs_present)

# Use facets by specification
p_specs <- ggplot(all_specs, aes(x = horizon, y = irf_mean, color = category, fill = category)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~spec, ncol = n_specs) +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_fill_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x     = "Horizon (years)",
    y     = "GDP response to maxwind_95 shock",
    color = "Period", fill = "Period",
    title = "Period-Split IRFs: Standard vs Flexible Year FE Specifications",
    subtitle = "Testing whether heterogeneous regional trends drive the post-1990 divergence"
  ) +
  irf_theme +
  theme(legend.position = "top")

ggsave("investigation/figures/irf_continent_year_fe.png", p_specs,
       width = 14, height = 5, dpi = 300)
cat("Saved: investigation/figures/irf_continent_year_fe.png\n")


# ==================================================================
# FINAL SUMMARY
# ==================================================================
cat("\n\n============================================================\n")
cat("  OVERALL SUMMARY\n")
cat("============================================================\n\n")

cat("1. YEAR FE DECOMPOSITION:\n")
cat(sprintf("   Correlation between full-sample and ever-treated year FE: %.4f\n", cor_val))
cat("   If the year FE diverge substantially post-1990, the implicit\n")
cat("   counterfactual differs between the two samples.\n\n")

cat("2. COUNTERFACTUAL GDP PATHS:\n")
cat(sprintf("   Never-treated GDP growth: Pre=%.4f, Post=%.4f, change=%.4f\n",
            gdp_summary$mean_growth[gdp_summary$group == "Never-treated" & gdp_summary$period == "Pre-1990"],
            gdp_summary$mean_growth[gdp_summary$group == "Never-treated" & gdp_summary$period == "Post-1990"],
            growth_did$change[growth_did$group == "Never-treated"]))
cat(sprintf("   Ever-treated GDP growth:  Pre=%.4f, Post=%.4f, change=%.4f\n",
            gdp_summary$mean_growth[gdp_summary$group == "Ever-treated" & gdp_summary$period == "Pre-1990"],
            gdp_summary$mean_growth[gdp_summary$group == "Ever-treated" & gdp_summary$period == "Post-1990"],
            growth_did$change[growth_did$group == "Ever-treated"]))
cat("   If these paths diverge, the year FE absorb different things\n")
cat("   for the two groups.\n\n")

cat("3. PLACEBO TEST:\n")
cat(sprintf("   Actual post-pre gap at h=5: %.4f\n", actual_gap))
cat(sprintf("   Placebo gap mean: %.4f, SD: %.4f\n",
            mean(placebo_df$gap), sd(placebo_df$gap)))
cat(sprintf("   P-value: %.4f\n", p_val))
if (p_val < 0.05) {
  cat("   => The actual gap is UNLIKELY to arise by chance in never-treated data.\n")
  cat("   => The post-1990 divergence is NOT a pure artifact of the year FE.\n\n")
} else {
  cat("   => The actual gap CAN arise by chance in never-treated data.\n")
  cat("   => The year FE + never-treated composition may partially drive the result.\n\n")
}

cat("4. DIFFERENTIAL TRENDS TEST:\n")
cat(sprintf("   post1990 x ever_treated coefficient: %.5f (SE=%.5f, t=%.3f)\n",
            did_coef, did_se, did_t))
if (abs(did_t) > 1.96) {
  cat("   => SIGNIFICANT differential trends between groups.\n")
  cat("   => The year FE cannot fully absorb the trend difference.\n\n")
} else {
  cat("   => NO significant differential trends (conditional on FE).\n")
  cat("   => Year FE appropriately absorb common trends.\n\n")
}

cat("5. FLEXIBLE COUNTERFACTUAL (Continent x Year FE):\n")
for (sp in unique(gap_table$spec)) {
  gaps <- gap_table %>% filter(spec == sp)
  for (i in 1:nrow(gaps)) {
    cat(sprintf("   %s, h=%d: gap=%.4f\n", sp, gaps$horizon[i], gaps$gap[i]))
  }
}
cat("   If continent x year FE substantially reduces the post-pre gap,\n")
cat("   heterogeneous regional trends (correlated with treatment) drive the result.\n")

cat("\n\nFigures saved:\n")
cat("  - investigation/figures/year_fe_decomposition.png\n")
cat("  - investigation/figures/counterfactual_gdp_paths.png\n")
cat("  - investigation/figures/placebo_never_treated.png\n")
cat("  - investigation/figures/irf_continent_year_fe.png\n")

cat("\nDone.\n")
