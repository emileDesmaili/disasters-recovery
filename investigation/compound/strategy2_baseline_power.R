# strategy2_baseline_power.R
# ------------------------------------------------------------------
# Robustness grid for the Strategy-2 baseline average-IRF (spec E0).
# Goal: identify specifications where beta(h) is more cleanly
# estimated, so that downstream interaction analyses have more power.
#
# Variants grouped by lever:
#   Treatment:  Cat 2+, Cat 3+, Cat 1+/p95, continuous wind, scaled
#   Sample/FE:  same-continent, treated-only, all+region-year FE
#   Functional: scalar beta vs quadratic beta(h)
#
# Each variant reports beta at h = 5 with clustered SE. Output is a
# coefplot + IRF overlay.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(haven)
  library(countrycode)
})

setFixest_notes(FALSE)

out_fig <- "figures/strategy2_power"
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

# ==================================================================
# 1. DATA
# ==================================================================
knot_to_ms  <- 0.514444
H           <- 10

pwt <- read_stata("../../raw_data/pwt_clean.dta")
tcs <- read_stata("../../raw_data/ibtracs_clean.dta")

d_full <- pwt %>%
  left_join(tcs, by = c("year", "countrycode")) %>%
  mutate(
    maxwind   = replace_na(max_ann_wind_i_nob, 0) * knot_to_ms,
    loggdp    = ln_real_gdp_usd_pc,
    gdp_diff  = growth_real_gdp_usd_pc,
    year2     = year^2
  ) %>%
  filter(year >= 1969) %>%
  mutate(
    continent = countrycode(countrycode, "iso3c", "continent"),
    continent = ifelse(is.na(continent), "Other", continent),
    region    = countrycode(countrycode, "iso3c", "region"),
    region    = ifelse(is.na(region), "Other", region),
    shock_cat1 = as.integer(maxwind >= 64 * knot_to_ms),  # true Cat 1+
    shock_cat2 = as.integer(maxwind >= 83 * knot_to_ms),
    shock_cat3 = as.integer(maxwind >= 96 * knot_to_ms),
    maxwind_10 = maxwind / 10  # per 10 m/s
  )

treated_cc   <- d_full %>% filter(shock_cat2 == 1) %>% pull(countrycode) %>% unique()
treated_cont <- d_full %>% filter(countrycode %in% treated_cc) %>%
  pull(continent) %>% unique()

d_samecont <- d_full %>% filter(continent %in% treated_cont)
d_treated  <- d_full %>% filter(countrycode %in% treated_cc)
d_all      <- d_full

cat(sprintf("Samples: treated=%d, same-continent=%d, all=%d (obs)\n",
            nrow(d_treated), nrow(d_samecont), nrow(d_all)))
cat(sprintf("Shock counts (same-continent): Cat1=%d Cat2=%d Cat3=%d\n",
            sum(d_samecont$shock_cat1), sum(d_samecont$shock_cat2),
            sum(d_samecont$shock_cat3)))

# ==================================================================
# 2. STACKED-PANEL BUILDER
# ==================================================================
build_stacked <- function(dat, treat_var) {
  out <- vector("list", H + 1)
  for (h in 0:H) {
    out[[h + 1]] <- dat %>%
      group_by(countrycode) %>%
      arrange(year) %>%
      mutate(
        dy          = lead(loggdp, h) - lag(loggdp, 1),
        gdp_diff_l1 = lag(gdp_diff, 1),
        gdp_diff_l2 = lag(gdp_diff, 2),
        shock_l1    = lag(.data[[treat_var]], 1),
        shock_l2    = lag(.data[[treat_var]], 2)
      ) %>%
      ungroup() %>%
      filter(!is.na(dy), !is.na(gdp_diff_l1), !is.na(gdp_diff_l2)) %>%
      mutate(h = h,
             T_ = .data[[treat_var]])
  }
  bind_rows(out) %>%
    mutate(
      fh = factor(h),
      h1 = h, h2 = h^2,
      T_h0 = T_,
      T_h1 = T_ * h1,
      T_h2 = T_ * h2,
      country_h = interaction(countrycode, fh, drop = TRUE),
      region_h  = interaction(region,      fh, drop = TRUE)
    )
}

# ==================================================================
# 3. ESTIMATORS
# ==================================================================
fit_quadratic <- function(sd, fe_formula, label) {
  fml <- as.formula(paste0(
    "dy ~ T_h0 + T_h1 + T_h2",
    " + i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2)",
    " + i(fh, shock_l1) + i(fh, shock_l2) | ", fe_formula
  ))
  mod <- feols(fml, data = sd, vcov = ~countrycode)
  g <- coef(mod)[c("T_h0","T_h1","T_h2")]
  V <- vcov(mod)[c("T_h0","T_h1","T_h2"), c("T_h0","T_h1","T_h2"), drop = FALSE]
  # beta at h = 5
  x <- c(1, 5, 25)
  b <- sum(x * g); se <- as.numeric(sqrt(t(x) %*% V %*% x))
  # smooth IRF curve
  h_seq <- seq(0, H, by = 0.25)
  curve <- data.frame(
    horizon = h_seq,
    irf = sapply(h_seq, function(hh) sum(c(1, hh, hh^2) * g)),
    se  = sapply(h_seq, function(hh) {
      xv <- c(1, hh, hh^2)
      as.numeric(sqrt(t(xv) %*% V %*% xv))
    })
  )
  list(label = label, beta5 = b, se5 = se, curve = curve,
       poly = "quadratic", nobs = nobs(mod))
}

fit_scalar <- function(sd, fe_formula, label) {
  fml <- as.formula(paste0(
    "dy ~ T_h0",
    " + i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2)",
    " + i(fh, shock_l1) + i(fh, shock_l2) | ", fe_formula
  ))
  mod <- feols(fml, data = sd, vcov = ~countrycode)
  b <- coef(mod)["T_h0"]; se <- sqrt(vcov(mod)["T_h0","T_h0"])
  curve <- data.frame(horizon = seq(0, H, by = 0.25),
                      irf = b, se = se)
  list(label = label, beta5 = b, se5 = se, curve = curve,
       poly = "scalar", nobs = nobs(mod))
}

fe_base      <- "country_h[year] + country_h[year2] + year^fh"
fe_region    <- "country_h[year] + country_h[year2] + region_h^year"
# Note: same-continent already includes ~all 182 countries (every
# continent contains at least one Cat 2+ treated country), so a
# pure "all-countries" variant is identical to same-continent.
# Use continent-year FE as a coarser-than-region alternative.
fe_contin_yr <- "country_h[year] + country_h[year2] + continent^year^fh"

# ==================================================================
# 4. RUN VARIANT GRID
# ==================================================================
runs <- list()

cat("\n--- A. Treatment variable (same-continent, quadratic) ---\n")
runs[["A1 Cat 2+ (reference)"]]   <- fit_quadratic(build_stacked(d_samecont, "shock_cat2"), fe_base,    "A1 Cat 2+ (reference)")
runs[["A2 Cat 3+"]]                <- fit_quadratic(build_stacked(d_samecont, "shock_cat3"), fe_base,    "A2 Cat 3+")
runs[["A3 Cat 1+ (64 kts)"]]       <- fit_quadratic(build_stacked(d_samecont, "shock_cat1"), fe_base,    "A3 Cat 1+ (64 kts)")
runs[["A4 Continuous wind (m/s)"]] <- fit_quadratic(build_stacked(d_samecont, "maxwind"),    fe_base,    "A4 Continuous wind (m/s)")
runs[["A5 Continuous wind / 10"]]  <- fit_quadratic(build_stacked(d_samecont, "maxwind_10"), fe_base,    "A5 Continuous wind / 10")

cat("\n--- B. Sample / FE (Cat 2+, quadratic) ---\n")
runs[["B1 Treated only"]]              <- fit_quadratic(build_stacked(d_treated,  "shock_cat2"), fe_base,      "B1 Treated only")
runs[["B2 + region-year FE"]]          <- fit_quadratic(build_stacked(d_samecont, "shock_cat2"), fe_region,    "B2 + region-year FE")
runs[["B3 + continent-year-h FE"]]     <- fit_quadratic(build_stacked(d_samecont, "shock_cat2"), fe_contin_yr, "B3 + continent-year-h FE")

cat("\n--- C. Functional form ---\n")
runs[["C1 Cat 2+ scalar beta"]]             <- fit_scalar(build_stacked(d_samecont, "shock_cat2"), fe_base, "C1 Cat 2+ scalar beta")
runs[["C2 Continuous wind scalar beta"]]    <- fit_scalar(build_stacked(d_samecont, "maxwind"),    fe_base, "C2 Continuous wind scalar beta")
runs[["C3 Continuous wind / 10 scalar"]]    <- fit_scalar(build_stacked(d_samecont, "maxwind_10"), fe_base, "C3 Continuous wind / 10 scalar")

cat("\n--- D. Combined: most-powerful candidates ---\n")
runs[["D1 Cat 1+ treated-only"]]            <- fit_quadratic(build_stacked(d_treated, "shock_cat1"), fe_base, "D1 Cat 1+ treated-only")
runs[["D2 Cat 1+ scalar beta"]]             <- fit_scalar(build_stacked(d_samecont, "shock_cat1"),   fe_base, "D2 Cat 1+ scalar beta")
runs[["D3 Cat 1+ treated + scalar"]]        <- fit_scalar(build_stacked(d_treated,  "shock_cat1"),   fe_base, "D3 Cat 1+ treated + scalar")

# Summary table
summary_df <- bind_rows(lapply(runs, function(r) {
  data.frame(label = r$label, beta5 = r$beta5, se5 = r$se5,
             t_stat = r$beta5 / r$se5, poly = r$poly, nobs = r$nobs)
})) %>% arrange(se5)

print(summary_df, row.names = FALSE)

# ==================================================================
# 5. FIGURES
# ==================================================================
# --- Coefplot of beta(h=5), sorted by SE ---
coef_df <- summary_df %>%
  mutate(label = factor(label, levels = rev(label)),
         lo = beta5 - 1.96 * se5, hi = beta5 + 1.96 * se5,
         group = substr(label, 1, 1))

p_coef <- ggplot(coef_df, aes(label, beta5, color = group)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.6) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_color_manual(values = c("A" = "#d73027", "B" = "#4575b4",
                                "C" = "#7b3294", "D" = "#1a9850"),
                     labels = c("A" = "Treatment variant",
                                "B" = "Sample / FE",
                                "C" = "Functional form",
                                "D" = "Combined (best-power)"),
                     name = NULL) +
  labs(x = NULL, y = expression(beta(h == 5)),
       title = "Baseline average IRF at h = 5: variants sorted by SE",
       subtitle = "Quadratic and scalar specs on the same axis. Narrower CI = more power.") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(out_fig, "beta5_coefplot.png"), p_coef,
       width = 11, height = 6, dpi = 200)

# --- IRF overlay: quadratic specs only (continuous-wind specs on a
#      separate panel because units differ) ---
binary_labels <- c("A1 Cat 2+ (reference)", "A2 Cat 3+", "A3 Cat 1+ (64 kts)",
                   "B1 Treated only", "B2 + region-year FE",
                   "B3 + continent-year-h FE", "D1 Cat 1+ treated-only")
cont_labels   <- c("A4 Continuous wind (m/s)", "A5 Continuous wind / 10")

bin_curves <- bind_rows(lapply(runs[binary_labels], function(r) {
  r$curve %>% mutate(label = r$label)
})) %>%
  mutate(lo = irf - 1.96 * se, hi = irf + 1.96 * se)

p_bin <- ggplot(bin_curves, aes(horizon, irf)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  facet_wrap(~label, ncol = 3) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Horizon (years)", y = "Cumulative effect on log GDP pc",
       title = "Baseline IRFs: binary treatments across sample/FE variants") +
  theme_classic(base_size = 11)
ggsave(file.path(out_fig, "irf_binary_variants.png"), p_bin,
       width = 13, height = 7, dpi = 200)

cont_curves <- bind_rows(lapply(runs[cont_labels], function(r) {
  r$curve %>% mutate(label = r$label)
})) %>%
  mutate(lo = irf - 1.96 * se, hi = irf + 1.96 * se)

p_cont <- ggplot(cont_curves, aes(horizon, irf)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "darkorange", alpha = 0.25) +
  geom_line(color = "darkorange", linewidth = 0.8) +
  facet_wrap(~label, ncol = 2, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "Horizon (years)", y = "Cumulative effect on log GDP pc per unit wind",
       title = "Baseline IRFs: continuous wind treatment")
ggsave(file.path(out_fig, "irf_continuous_variants.png"), p_cont,
       width = 11, height = 4.5, dpi = 200)

# ==================================================================
# 6. DIAGNOSTICS
# ==================================================================
diag_lines <- c(
  "Strategy 2 baseline power grid",
  "=================================",
  sprintf("Sample sizes: treated=%d obs, same-continent=%d obs, all=%d obs",
          nrow(d_treated), nrow(d_samecont), nrow(d_all)),
  sprintf("Same-continent shock counts: Cat1=%d Cat2=%d Cat3=%d",
          sum(d_samecont$shock_cat1), sum(d_samecont$shock_cat2),
          sum(d_samecont$shock_cat3)),
  "",
  "beta(h=5) across variants (sorted by SE, smallest first):",
  paste0("  ", capture.output(print(summary_df, row.names = FALSE)))
)
writeLines(diag_lines, "strategy2_baseline_power.txt")

saveRDS(list(runs = runs, summary = summary_df),
        "strategy2_baseline_power.rds")

cat("\nSaved: strategy2_baseline_power.txt\n")
cat("Saved: strategy2_baseline_power.rds\nDone.\n")
