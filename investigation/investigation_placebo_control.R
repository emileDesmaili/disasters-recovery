# investigation_placebo_control.R
# Control group resampling placebo test
#
# The hypothesis: the large post-pre gap in the period-interacted LP is driven
# by control group composition. 142/182 countries never experience cyclone
# shocks, and they are geographically concentrated (Europe/Africa) with
# different GDP trends from treated countries (Americas/Asia).
#
# The test:
#   1. Keep the 40 treated countries (ever-treated by p95 maxwind shock) and
#      their REAL shocks fixed.
#   2. For each simulation, randomly sample a subset of never-treated countries
#      to serve as the control group.
#   3. Combine treated + sampled control -> run the period-interacted LP at h=5
#      -> record the post-pre gap.
#   4. Repeat 500 times for each control group size (30, 70, 100 of 142).
#   5. Plot the distribution of gaps.
#
# Prediction: if control group composition matters, different random subsets
# should produce substantially different gaps, and the spread should increase
# as the sample size shrinks.

library(tidyverse)
library(fixest)
library(haven)
source("../emileRegs.R")

setFixest_notes(FALSE)

set.seed(42)

cat("============================================================\n")
cat("  CONTROL GROUP RESAMPLING PLACEBO TEST\n")
cat("============================================================\n\n")

# ------------------------------------------------------------------
# 1. Load and construct data (identically to lp_gdp_shocks.R)
# ------------------------------------------------------------------
pwt <- read_stata("../raw_data/pwt_clean.dta")
tcs <- read_stata("../raw_data/ibtracs_clean.dta")
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
    maxwind_p95 = quantile(maxwind[year >= 1970 & year <= 2014], 0.95, na.rm = TRUE),
    energy_p95  = quantile(energy[year >= 1970 & year <= 2014],  0.95, na.rm = TRUE),
    nlands_p95  = quantile(nlands[year >= 1970 & year <= 2014],  0.95, na.rm = TRUE)
  ) %>%
  mutate(
    maxwind_95 = as.integer(maxwind_p95 > 0 & maxwind >= maxwind_p95),
    energy_95  = as.integer(energy_p95  > 0 & energy  >= energy_p95),
    nlands_95  = as.integer(nlands_p95  > 0 & nlands  >= nlands_p95)
  )

# ------------------------------------------------------------------
# 2. Identify treated vs never-treated countries
# ------------------------------------------------------------------
ever_treated <- data %>%
  filter(maxwind_95 == 1) %>%
  pull(countrycode) %>%
  unique() %>%
  sort()

all_countries <- unique(data$countrycode)
never_treated <- sort(setdiff(all_countries, ever_treated))

cat(sprintf("Total countries:        %d\n", length(all_countries)))
cat(sprintf("Ever-treated:           %d\n", length(ever_treated)))
cat(sprintf("Never-treated:          %d\n", length(never_treated)))

# Split data
data_treated     <- data %>% filter(countrycode %in% ever_treated)
data_never       <- data %>% filter(countrycode %in% never_treated)

# ------------------------------------------------------------------
# 3. LP settings
# ------------------------------------------------------------------
outcome  <- "loggdp"
h_target <- 5
controls <- "l(gdp_diff,1:2) + l(maxwind_95,1:2)"
fe       <- "countrycode[year] + countrycode[year2] + year"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

# ------------------------------------------------------------------
# 4. Compute the actual full-sample gap at h=5
# ------------------------------------------------------------------
cat("\nEstimating full-sample period-interacted LP at h=5...\n")

irf_full <- lp_panel_inter(
  data     = data,
  outcome  = outcome,
  main_var = "maxwind_95",
  interact_var = "period",
  controls = controls,
  horizon  = h_target,
  fe       = fe,
  panel_id = panel_id,
  vcov_formula = vcov_fm
)

# Extract post-pre gap at h=5
irf_h5 <- irf_full %>%
  filter(horizon == h_target) %>%
  mutate(
    category = case_when(
      str_detect(category, "Pre")  ~ "Pre-1990",
      str_detect(category, "Post") ~ "Post-1990",
      TRUE ~ category
    )
  )

post_coef <- irf_h5 %>% filter(category == "Post-1990") %>% pull(irf_mean)
pre_coef  <- irf_h5 %>% filter(category == "Pre-1990")  %>% pull(irf_mean)
actual_gap <- post_coef - pre_coef

cat(sprintf("\nFull-sample results at h=%d:\n", h_target))
cat(sprintf("  Pre-1990 coefficient:   %.4f (%.2f%%)\n", pre_coef, pre_coef * 100))
cat(sprintf("  Post-1990 coefficient:  %.4f (%.2f%%)\n", post_coef, post_coef * 100))
cat(sprintf("  Post - Pre gap:         %.4f (%.2f pp)\n", actual_gap, actual_gap * 100))

# ------------------------------------------------------------------
# 5. Helper function: run one simulation
# ------------------------------------------------------------------
run_one_sim <- function(n_control, data_treated, data_never, never_treated) {
  # Sample n_control never-treated countries (without replacement)
  sampled_controls <- sample(never_treated, size = n_control, replace = FALSE)

  # Combine treated + sampled control
  data_sim <- bind_rows(
    data_treated,
    data_never %>% filter(countrycode %in% sampled_controls)
  )

  # Run period-interacted LP at h=5 only
  irf_sim <- lp_panel_inter(
    data     = data_sim,
    outcome  = outcome,
    main_var = "maxwind_95",
    interact_var = "period",
    controls = controls,
    horizon  = h_target,
    fe       = fe,
    panel_id = panel_id,
    vcov_formula = vcov_fm
  )

  # Extract post-pre gap at h=5
  irf_h5_sim <- irf_sim %>%
    filter(horizon == h_target) %>%
    mutate(
      category = case_when(
        str_detect(category, "Pre")  ~ "Pre-1990",
        str_detect(category, "Post") ~ "Post-1990",
        TRUE ~ category
      )
    )

  post_val <- irf_h5_sim %>% filter(category == "Post-1990") %>% pull(irf_mean)
  pre_val  <- irf_h5_sim %>% filter(category == "Pre-1990")  %>% pull(irf_mean)

  if (length(post_val) == 1 & length(pre_val) == 1) {
    return(post_val - pre_val)
  } else {
    return(NA_real_)
  }
}

# ------------------------------------------------------------------
# 6. Run simulations for each control group size
# ------------------------------------------------------------------
n_sims <- 500
control_sizes <- c(30, 70, 100)

results_list <- list()

for (n_ctrl in control_sizes) {
  cat(sprintf("\nRunning %d simulations with %d/%d never-treated countries...\n",
              n_sims, n_ctrl, length(never_treated)))

  gaps <- numeric(n_sims)

  for (s in 1:n_sims) {
    if (s %% 50 == 0) cat(sprintf("  sim %d/%d\n", s, n_sims))

    gaps[s] <- tryCatch(
      run_one_sim(n_ctrl, data_treated, data_never, never_treated),
      error = function(e) {
        cat(sprintf("  [Warning] sim %d failed: %s\n", s, e$message))
        return(NA_real_)
      }
    )
  }

  results_list[[as.character(n_ctrl)]] <- tibble(
    sim           = 1:n_sims,
    gap           = gaps,
    control_size  = n_ctrl,
    control_label = sprintf("%d of %d", n_ctrl, length(never_treated))
  )
}

results <- bind_rows(results_list)

# ------------------------------------------------------------------
# 7. Summary statistics
# ------------------------------------------------------------------
cat("\n============================================================\n")
cat("  SUMMARY STATISTICS\n")
cat("============================================================\n\n")

cat(sprintf("Actual full-sample gap at h=%d: %.4f (%.2f pp)\n\n",
            h_target, actual_gap, actual_gap * 100))

summary_stats <- results %>%
  filter(!is.na(gap)) %>%
  group_by(control_label) %>%
  summarise(
    n_valid  = n(),
    mean_gap = mean(gap),
    sd_gap   = sd(gap),
    min_gap  = min(gap),
    p05_gap  = quantile(gap, 0.05),
    p25_gap  = quantile(gap, 0.25),
    median_gap = median(gap),
    p75_gap  = quantile(gap, 0.75),
    p95_gap  = quantile(gap, 0.95),
    max_gap  = max(gap),
    .groups  = "drop"
  )

cat("Distribution of post-pre gap by control group size:\n\n")
print(as.data.frame(summary_stats), row.names = FALSE)

cat(sprintf("\nInterpretation:\n"))
cat(sprintf("  If the SD of the gap is large relative to the actual gap (%.2f pp),\n",
            actual_gap * 100))
cat(sprintf("  then control group composition is an important driver.\n"))

for (i in seq_len(nrow(summary_stats))) {
  row <- summary_stats[i, ]
  cat(sprintf("\n  %s: SD = %.2f pp, range = [%.2f, %.2f] pp\n",
              row$control_label,
              row$sd_gap * 100,
              row$min_gap * 100,
              row$max_gap * 100))
}

# ------------------------------------------------------------------
# 8. Plot: distribution of gaps for each control group size
# ------------------------------------------------------------------
cat("\nGenerating figure...\n")

results_plot <- results %>%
  filter(!is.na(gap)) %>%
  mutate(
    control_label = factor(
      control_label,
      levels = sprintf("%d of %d", control_sizes, length(never_treated))
    )
  )

p <- ggplot(results_plot, aes(x = gap * 100)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "steelblue", alpha = 0.6, color = "white") +
  geom_density(linewidth = 0.8, color = "navy") +
  geom_vline(xintercept = actual_gap * 100, color = "firebrick",
             linewidth = 1.2, linetype = "dashed") +
  annotate("text",
           x = actual_gap * 100, y = Inf,
           label = sprintf("Full-sample gap\n%.2f pp", actual_gap * 100),
           color = "firebrick", hjust = -0.1, vjust = 1.5, size = 3.5,
           fontface = "bold") +
  facet_wrap(~control_label, ncol = 1, scales = "free_y",
             strip.position = "top") +
  labs(
    x = "Post-1990 minus Pre-1990 gap at h=5 (percentage points)",
    y = "Density",
    title = "Control Group Resampling Placebo Test",
    subtitle = paste0(
      "Distribution of post-pre gap across ", n_sims,
      " random subsets of never-treated countries\n",
      "Red dashed line = actual full-sample gap"
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey30")
  )

ggsave("figures/placebo_control_resampling.png", p,
       width = 8, height = 10, dpi = 300)
cat("Saved: figures/placebo_control_resampling.png\n")

cat("\nDone.\n")
