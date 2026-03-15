# investigation_placebo_control_regionfe.R
# Control group resampling placebo test with region x year FE
#
# Same design as investigation_placebo_control.R but using region^year FE
# instead of plain year FE. Also tests whether unit-specific polynomial
# controls (countrycode[year] + countrycode[year2]) matter with region^year FE.

library(tidyverse)
library(fixest)
library(haven)
library(countrycode)

source("../emileRegs.R")

setFixest_notes(FALSE)
set.seed(42)

cat("============================================================\n")
cat("  CONTROL GROUP RESAMPLING (REGION x YEAR FE)\n")
cat("============================================================\n\n")

# ------------------------------------------------------------------
# 1. Load and construct data
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
    maxwind_p95 = quantile(maxwind[year >= 1970 & year <= 2014], 0.95, na.rm = TRUE)
  ) %>%
  mutate(
    maxwind_95 = as.integer(maxwind_p95 > 0 & maxwind >= maxwind_p95)
  )

# Add region variable
data <- data %>%
  mutate(
    region = countrycode(countrycode, "iso3c", "region"),
    region = ifelse(is.na(region), "Other", region)
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

data_treated <- data %>% filter(countrycode %in% ever_treated)
data_never   <- data %>% filter(countrycode %in% never_treated)

# ------------------------------------------------------------------
# 3. LP settings — two FE specifications
# ------------------------------------------------------------------
outcome  <- "loggdp"
h_target <- 5
controls <- "l(gdp_diff,1:2) + l(maxwind_95,1:2)"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

fe_specs <- list(
  "Region x Year + Country Polynomials" = "countrycode[year] + countrycode[year2] + region^year",
  "Region x Year Only"                  = "region^year"
)

# ------------------------------------------------------------------
# 4. Helper: extract post-pre gap
# ------------------------------------------------------------------
extract_gap <- function(irf_df, h = 5) {
  irf_df <- irf_df %>%
    filter(horizon == h) %>%
    mutate(
      category = case_when(
        str_detect(category, "Pre")  ~ "Pre-1990",
        str_detect(category, "Post") ~ "Post-1990",
        TRUE ~ category
      )
    )
  post <- irf_df %>% filter(category == "Post-1990") %>% pull(irf_mean)
  pre  <- irf_df %>% filter(category == "Pre-1990")  %>% pull(irf_mean)
  if (length(post) == 1 & length(pre) == 1) {
    return(post - pre)
  } else {
    return(NA_real_)
  }
}

# ------------------------------------------------------------------
# 5. Run for each FE specification
# ------------------------------------------------------------------
n_sims <- 500
control_sizes <- c(30, 70, 100)
all_results <- list()

for (spec_name in names(fe_specs)) {
  fe <- fe_specs[[spec_name]]
  cat(sprintf("\n=== %s ===\n", spec_name))
  cat(sprintf("  FE: %s\n", fe))

  # Actual full-sample gap
  cat("  Computing actual gap...\n")
  irf_full <- lp_panel_inter(
    data = data, outcome = outcome, main_var = "maxwind_95",
    interact_var = "period", controls = controls,
    horizon = h_target, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
  )
  actual_gap <- extract_gap(irf_full)
  cat(sprintf("  Actual gap: %.4f (%.2f pp)\n", actual_gap, actual_gap * 100))

  for (n_ctrl in control_sizes) {
    cat(sprintf("\n  Running %d sims with %d/%d controls...\n",
                n_sims, n_ctrl, length(never_treated)))
    t_start <- Sys.time()
    gaps <- numeric(n_sims)

    for (s in 1:n_sims) {
      if (s %% 100 == 0) {
        elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
        cat(sprintf("    Sim %d/%d  (%.1f min)\n", s, n_sims, elapsed))
      }

      sampled_controls <- sample(never_treated, size = n_ctrl, replace = FALSE)
      data_sim <- bind_rows(
        data_treated,
        data_never %>% filter(countrycode %in% sampled_controls)
      )

      gaps[s] <- tryCatch({
        irf_sim <- lp_panel_inter(
          data = data_sim, outcome = outcome, main_var = "maxwind_95",
          interact_var = "period", controls = controls,
          horizon = h_target, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
        )
        extract_gap(irf_sim)
      }, error = function(e) NA_real_)
    }

    elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    valid_gaps <- gaps[!is.na(gaps)]
    cat(sprintf("  Done (%.1f min). Mean: %.2f pp, SD: %.2f pp\n",
                elapsed_total, mean(valid_gaps) * 100, sd(valid_gaps) * 100))

    all_results[[paste(spec_name, n_ctrl)]] <- tibble(
      sim           = seq_along(valid_gaps),
      gap           = valid_gaps,
      control_size  = n_ctrl,
      control_label = sprintf("%d of %d", n_ctrl, length(never_treated)),
      spec          = spec_name,
      actual_gap    = actual_gap
    )
  }
}

results <- bind_rows(all_results)

# ------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------
cat("\n============================================================\n")
cat("  SUMMARY\n")
cat("============================================================\n\n")

summary_stats <- results %>%
  group_by(spec, control_label) %>%
  summarise(
    actual   = first(actual_gap) * 100,
    n_valid  = n(),
    mean_gap = mean(gap) * 100,
    sd_gap   = sd(gap) * 100,
    min_gap  = min(gap) * 100,
    max_gap  = max(gap) * 100,
    .groups  = "drop"
  )

print(as.data.frame(summary_stats), row.names = FALSE)

# ------------------------------------------------------------------
# 7. Plot
# ------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

results_plot <- results %>%
  mutate(
    spec = factor(spec, levels = names(fe_specs)),
    control_label = factor(
      control_label,
      levels = sprintf("%d of %d", control_sizes, length(never_treated))
    )
  )

actual_df <- results_plot %>%
  distinct(spec, actual_gap)

p <- ggplot(results_plot, aes(x = gap * 100)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "steelblue", alpha = 0.6, color = "white") +
  geom_density(linewidth = 0.8, color = "navy") +
  geom_vline(data = actual_df, aes(xintercept = actual_gap * 100),
             color = "firebrick", linewidth = 1.2, linetype = "dashed") +
  facet_grid(control_label ~ spec, scales = "free") +
  labs(
    x = "Post-1990 minus Pre-1990 gap at h=5 (percentage points)",
    y = "Density",
    title = "Control Group Resampling with Region x Year FE",
    subtitle = sprintf(
      "%d simulations per cell | Red dashed line = actual full-sample gap",
      n_sims
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey30")
  )

ggsave("figures/placebo_control_regionfe.png", p,
       width = 12, height = 10, dpi = 300)
cat("\nSaved: figures/placebo_control_regionfe.png\n")
cat("Done.\n")
