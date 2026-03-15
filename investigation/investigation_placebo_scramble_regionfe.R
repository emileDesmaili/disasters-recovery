# investigation_placebo_scramble_regionfe.R
# Shock scrambling placebo test with region x year FE
#
# Same design as investigation_placebo_scramble.R but using region^year FE
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
cat("  PLACEBO TEST: SHOCK SCRAMBLING (REGION x YEAR FE)\n")
cat("============================================================\n\n")

# ------------------------------------------------------------------
# 1. Load and construct data (identically to lp_gdp_shocks.R)
# ------------------------------------------------------------------
pwt  <- read_stata("../raw_data/pwt_clean.dta")
tcs  <- read_stata("../raw_data/ibtracs_clean.dta")
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
# 2. LP settings — two FE specifications
# ------------------------------------------------------------------
outcome  <- "loggdp"
horizon  <- 5
controls <- "l(gdp_diff,1:2) + l(maxwind_95,1:2)"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

fe_specs <- list(
  "Region x Year + Country Polynomials" = "countrycode[year] + countrycode[year2] + region^year",
  "Region x Year Only"                  = "region^year"
)

# ------------------------------------------------------------------
# 3. Identify treated countries
# ------------------------------------------------------------------
ever_treated <- data %>%
  filter(maxwind_95 == 1) %>%
  pull(countrycode) %>%
  unique()

treated_idx <- which(data$countrycode %in% ever_treated)

cat(sprintf("Total countries:       %d\n", n_distinct(data$countrycode)))
cat(sprintf("Ever-treated:          %d\n", length(ever_treated)))
cat(sprintf("Treated country-years: %d\n", length(treated_idx)))
cat(sprintf("Total shocks:          %d\n\n", sum(data$maxwind_95)))

# ------------------------------------------------------------------
# 4. Helper: extract post-pre gap at h=5 from lp_panel_inter output
# ------------------------------------------------------------------
extract_gap <- function(irf_df) {
  irf_df <- irf_df %>%
    mutate(
      category = str_replace_all(category, "_", "-"),
      category = case_when(
        str_detect(category, "Pre")  ~ "Pre-1990",
        str_detect(category, "Post") ~ "Post-1990",
        TRUE ~ category
      )
    )
  h5 <- irf_df %>% filter(horizon == 5)
  post <- h5 %>% filter(category == "Post-1990") %>% pull(irf_mean)
  pre  <- h5 %>% filter(category == "Pre-1990")  %>% pull(irf_mean)
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
all_results <- list()

for (spec_name in names(fe_specs)) {
  fe <- fe_specs[[spec_name]]
  cat(sprintf("=== %s ===\n", spec_name))
  cat(sprintf("  FE: %s\n", fe))

  # Actual gap
  cat("  Computing actual gap...\n")
  irf_actual <- lp_panel_inter(
    data = data, outcome = outcome, main_var = "maxwind_95",
    interact_var = "period", controls = controls,
    horizon = horizon, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
  )
  actual_gap <- extract_gap(irf_actual)
  cat(sprintf("  Actual gap: %.4f (%.2f pp)\n", actual_gap, actual_gap * 100))

  # Scrambling
  cat(sprintf("  Running %d scrambling simulations...\n", n_sims))
  scrambled_gaps <- numeric(n_sims)
  t_start <- Sys.time()

  for (s in 1:n_sims) {
    if (s %% 100 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      cat(sprintf("    Sim %d/%d  (%.1f min)\n", s, n_sims, elapsed))
    }

    d_sim <- data
    d_sim$maxwind_95[treated_idx] <- sample(d_sim$maxwind_95[treated_idx])

    scrambled_gaps[s] <- tryCatch({
      irf_sim <- lp_panel_inter(
        data = d_sim, outcome = outcome, main_var = "maxwind_95",
        interact_var = "period", controls = controls,
        horizon = horizon, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
      )
      extract_gap(irf_sim)
    }, error = function(e) NA_real_)
  }

  elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  cat(sprintf("  Done. %.1f minutes\n", elapsed_total))

  valid_gaps <- scrambled_gaps[!is.na(scrambled_gaps)]
  p_value <- mean(abs(valid_gaps) >= abs(actual_gap))

  cat(sprintf("  Scrambled mean: %.4f (%.2f pp)\n", mean(valid_gaps), mean(valid_gaps) * 100))
  cat(sprintf("  Scrambled SD:   %.4f (%.2f pp)\n", sd(valid_gaps), sd(valid_gaps) * 100))
  cat(sprintf("  Two-sided p:    %.4f\n\n", p_value))

  all_results[[spec_name]] <- tibble(
    gap        = valid_gaps,
    spec       = spec_name,
    actual_gap = actual_gap,
    p_value    = p_value
  )
}

# ------------------------------------------------------------------
# 6. Plot: panel of histograms, one per FE specification
# ------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

plot_df <- bind_rows(all_results) %>%
  mutate(spec = factor(spec, levels = names(fe_specs)))

actual_df <- plot_df %>%
  distinct(spec, actual_gap, p_value)

p <- ggplot(plot_df, aes(x = gap)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "grey70", color = "grey40", linewidth = 0.3) +
  geom_vline(data = actual_df, aes(xintercept = actual_gap),
             color = "firebrick", linewidth = 1.2) +
  geom_text(data = actual_df,
            aes(x = actual_gap, y = Inf,
                label = sprintf("Actual = %.2f pp\np = %.3f",
                                actual_gap * 100, p_value)),
            color = "firebrick", fontface = "bold", size = 4,
            hjust = -0.1, vjust = 2) +
  facet_wrap(~spec, ncol = 1, scales = "free_y") +
  scale_x_continuous(labels = function(x) paste0(round(x * 100, 1), " pp")) +
  labs(
    x = "Post-1990 minus Pre-1990 gap at h = 5",
    y = "Density",
    title = "Shock Scrambling Placebo Test with Region x Year FE",
    subtitle = sprintf("%d simulations | Shocks permuted within %d treated countries",
                       n_sims, length(ever_treated))
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40")
  )

ggsave("figures/placebo_scramble_regionfe.png", p,
       width = 9, height = 8, dpi = 300)
cat("\nFigure saved: figures/placebo_scramble_regionfe.png\n")
cat("Done.\n")
