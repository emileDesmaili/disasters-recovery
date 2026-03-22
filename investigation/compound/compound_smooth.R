# compound_smooth.R
# Smooth parametric extensions:
#   Part A: Standard LP with continuous gap interaction (horizon-by-horizon)
#   Part B: Parametric (polynomial-in-h) IRF with gap interaction (stacked LP)
#   Part C: Country-specific quadratic parametric IRFs
#
# Each part is run for two samples:
#   (i)  Treated only (ever-hit countries)
#   (ii) Same-continent controls (all countries on a treated continent)

library(tidyverse)
library(fixest)
library(haven)
library(countrycode)

source("../../emileRegs.R")
setFixest_notes(FALSE)

# ==================================================================
# 1. BUILD BOTH SAMPLES
# ==================================================================

# --- Threshold: Category 2+ (83 knots) ---
knot_to_ms <- 0.514444
cat2_thresh <- 83 * knot_to_ms  # 42.7 m/s

# --- Full panel ---
pwt <- read_stata("../../raw_data/pwt_clean.dta")
tcs <- read_stata("../../raw_data/ibtracs_clean.dta")

d_full <- pwt %>%
  left_join(tcs, by = c("year", "countrycode")) %>%
  mutate(
    maxwind  = replace_na(max_ann_wind_i_nob, 0) * knot_to_ms,
    loggdp   = ln_real_gdp_usd_pc,
    gdp_diff = growth_real_gdp_usd_pc,
    year2    = year^2,
    period   = ifelse(year <= 1990, "Pre-1990", "Post-1990")
  ) %>%
  filter(year >= 1969) %>%
  mutate(
    maxwind_95  = as.integer(maxwind >= cat2_thresh),
    continent   = countrycode(countrycode, "iso3c", "continent"),
    continent   = ifelse(is.na(continent), "Other", continent),
    region      = countrycode(countrycode, "iso3c", "region"),
    region      = ifelse(is.na(region), "Other", region)
  )

cat(sprintf("Threshold: Cat 2+ (83 knots / %.1f m/s)\n", cat2_thresh))
cat(sprintf("Total shock-years: %d\n", sum(d_full$maxwind_95)))

# Identify treated countries and continents
treated_countries <- d_full %>%
  filter(maxwind_95 == 1) %>%
  pull(countrycode) %>%
  unique()

treated_continents <- d_full %>%
  filter(countrycode %in% treated_countries) %>%
  pull(continent) %>%
  unique()

# --- Treated sample: ever-hit countries only ---
d_treated <- d_full %>%
  filter(countrycode %in% treated_countries)

# --- Same-continent sample ---
d_continent <- d_full %>%
  filter(continent %in% treated_continents)

cat(sprintf("Treated countries: %d\n", length(treated_countries)))
cat(sprintf("Treated continents: %s\n", paste(treated_continents, collapse = ", ")))
cat(sprintf("Same-continent sample: %d countries, %d obs\n",
            n_distinct(d_continent$countrycode), nrow(d_continent)))

# --- Add shock lags and gap variable to both samples ---
add_gap_vars <- function(dat) {
  dat <- dat %>%
    group_by(countrycode) %>%
    arrange(year) %>%
    mutate(
      shock_lag1 = lag(maxwind_95, 1, default = 0L),
      shock_lag2 = lag(maxwind_95, 2, default = 0L)
    ) %>%
    ungroup()

  shock_gaps <- dat %>%
    filter(maxwind_95 == 1) %>%
    group_by(countrycode) %>%
    arrange(year) %>%
    mutate(gap = year - lag(year)) %>%
    ungroup()

  first_years <- dat %>%
    group_by(countrycode) %>%
    summarise(first_year = min(year), .groups = "drop")

  shock_gaps <- shock_gaps %>%
    left_join(first_years, by = "countrycode") %>%
    mutate(gap = ifelse(is.na(gap), year - first_year, gap)) %>%
    select(countrycode, year, gap)

  dat %>%
    left_join(shock_gaps, by = c("countrycode", "year")) %>%
    mutate(
      gap = replace_na(gap, 0),
      gap = pmax(gap, 1),
      gap_capped = pmin(gap, 20)
    )
}

d_treated   <- add_gap_vars(d_treated)
d_continent <- add_gap_vars(d_continent)

# Ensure d_treated has continent
if (!"continent" %in% names(d_treated)) {
  d_treated <- d_treated %>%
    mutate(
      continent = countrycode(countrycode, "iso3c", "continent"),
      continent = ifelse(is.na(continent), "Other", continent)
    )
}

samples <- list(
  "Treated only"   = d_treated,
  "Same continent" = d_continent
)

cat(sprintf("\nSample sizes: Treated = %d obs (%d countries), Continent = %d obs (%d countries)\n",
            nrow(d_treated), n_distinct(d_treated$countrycode),
            nrow(d_continent), n_distinct(d_continent$countrycode)))

# ==================================================================
# HELPER FUNCTIONS
# ==================================================================
H <- 10
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

# --- Part A helper: standard LP with gap interaction ---
run_standard_lp <- function(dat, label) {
  cat(sprintf("\n--- Standard LP: %s ---\n", label))
  irf_list <- list()
  coef_list <- list()

  for (h in 0:H) {
    fml <- as.formula(paste0(
      "f(loggdp, ", h, ") - l(loggdp, 1) ~ ",
      "maxwind_95 + maxwind_95:gap_capped + ",
      "l(gdp_diff, 1:2) + l(maxwind_95, 1:2) | ",
      "countrycode[year] + countrycode[year2] + year"
    ))

    mod <- feols(fml, data = dat, panel.id = panel_id, vcov = vcov_fm)

    cn <- names(coef(mod))
    shock_nm <- cn[cn == "maxwind_95"]
    inter_nm <- cn[grepl("maxwind_95.*gap_capped|gap_capped.*maxwind_95", cn)]

    b_shock <- coef(mod)[shock_nm]
    b_gap   <- coef(mod)[inter_nm]
    V <- vcov(mod)

    coef_list[[h + 1]] <- data.frame(
      horizon = h, beta_shock = b_shock, beta_gap = b_gap,
      se_shock = sqrt(V[shock_nm, shock_nm]),
      se_gap   = sqrt(V[inter_nm, inter_nm])
    )

    for (g in c(1, 3, 5, 10, 15)) {
      irf_val <- b_shock + b_gap * g
      irf_se  <- sqrt(V[shock_nm, shock_nm] + g^2 * V[inter_nm, inter_nm] +
                       2 * g * V[shock_nm, inter_nm])
      irf_list[[length(irf_list) + 1]] <- data.frame(
        horizon = h, gap = g, irf_mean = irf_val, se = irf_se,
        irf_down = irf_val - 1.96 * irf_se, irf_up = irf_val + 1.96 * irf_se
      )
    }
  }

  irf_df  <- bind_rows(irf_list) %>% mutate(sample = label)
  coef_df <- bind_rows(coef_list) %>% mutate(sample = label)

  h5 <- coef_df %>% filter(horizon == 5)
  cat(sprintf("  h=5: beta_gap = %.4f (SE = %.4f), gap=1: %.3f, gap=10: %.3f\n",
              h5$beta_gap, h5$se_gap,
              h5$beta_shock + h5$beta_gap * 1,
              h5$beta_shock + h5$beta_gap * 10))

  list(irf = irf_df, coefs = coef_df)
}

# --- Part B helper: stacked data + parametric estimation ---
build_stacked <- function(dat) {
  stacked_list <- vector("list", H + 1)
  for (h in 0:H) {
    stacked_list[[h + 1]] <- dat %>%
      group_by(countrycode) %>%
      arrange(year) %>%
      mutate(
        dy          = lead(loggdp, h) - lag(loggdp, 1),
        gdp_diff_l1 = lag(gdp_diff, 1),
        gdp_diff_l2 = lag(gdp_diff, 2),
        shock_l1    = lag(maxwind_95, 1),
        shock_l2    = lag(maxwind_95, 2)
      ) %>%
      ungroup() %>%
      filter(!is.na(dy), !is.na(gdp_diff_l1), !is.na(gdp_diff_l2)) %>%
      mutate(h = h)
  }
  bind_rows(stacked_list) %>%
    mutate(
      fh = factor(h),
      h1 = h, h2 = h^2, h3 = h^3,
      shock_h0 = maxwind_95,
      shock_h1 = maxwind_95 * h1,
      shock_h2 = maxwind_95 * h2,
      shock_h3 = maxwind_95 * h3,
      shock_gap_h0 = maxwind_95 * gap_capped,
      shock_gap_h1 = maxwind_95 * gap_capped * h1,
      shock_gap_h2 = maxwind_95 * gap_capped * h2,
      shock_gap_h3 = maxwind_95 * gap_capped * h3,
      country_h = interaction(countrycode, fh, drop = TRUE)
    )
}

estimate_parametric <- function(degree, stacked_df, label = "") {
  deg_label <- if (degree == 2) "Quadratic" else "Cubic"
  cat(sprintf("\n--- %s %s ---\n", deg_label, label))

  shock_base <- paste0("shock_h", 0:degree)
  shock_gap  <- paste0("shock_gap_h", 0:degree)
  sv <- c(shock_base, shock_gap)

  fml <- as.formula(paste0(
    "dy ~ ", paste(sv, collapse = " + "),
    " + i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2)",
    " + i(fh, shock_l1) + i(fh, shock_l2)",
    " | country_h[year] + country_h[year2] + year^fh"
  ))

  mod <- tryCatch({
    feols(fml, data = stacked_df, vcov = ~countrycode)
  }, error = function(e) {
    cat("  Note: Country-trend FE failed, trying without trends...\n")
    fml_s <- as.formula(paste0(
      "dy ~ ", paste(sv, collapse = " + "),
      " + i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2)",
      " + i(fh, shock_l1) + i(fh, shock_l2)",
      " | countrycode^fh + year^fh"
    ))
    feols(fml_s, data = stacked_df, vcov = ~countrycode)
  })

  theta <- coef(mod)[sv]
  V     <- vcov(mod)[sv, sv]

  # Predict smooth IRFs
  h_seq <- seq(0, H, by = 0.25)
  pred <- list()
  for (g in c(1, 3, 5, 10, 15)) {
    for (hh in h_seq) {
      x_base <- hh^(0:degree)
      x <- c(x_base, g * x_base)
      irf_val <- sum(x * theta)
      irf_se  <- as.numeric(sqrt(t(x) %*% V %*% x))
      pred[[length(pred) + 1]] <- data.frame(
        horizon = hh, gap = g, irf_mean = irf_val, se = irf_se,
        irf_down = irf_val - 1.96 * irf_se, irf_up = irf_val + 1.96 * irf_se,
        degree = degree
      )
    }
  }

  pred_df <- bind_rows(pred)

  cat(sprintf("  IRF at h = 5:\n"))
  for (g in c(1, 5, 10, 15)) {
    x_base <- 5^(0:degree)
    x <- c(x_base, g * x_base)
    val <- sum(x * theta)
    se  <- as.numeric(sqrt(t(x) %*% V %*% x))
    cat(sprintf("    gap = %2d: %.3f (SE = %.3f)\n", g, val, se))
  }

  list(pred = pred_df, theta = theta, V = V, model = mod)
}

# ==================================================================
# PART A: Standard LP with gap interaction — both samples
# ==================================================================
cat("\n=== Part A: Standard LP with gap interaction ===\n")

res_A <- lapply(names(samples), function(nm) run_standard_lp(samples[[nm]], nm))
names(res_A) <- names(samples)

irf_gap_df  <- bind_rows(lapply(res_A, `[[`, "irf"))
gap_coef_df <- bind_rows(lapply(res_A, `[[`, "coefs"))

# ==================================================================
# PART B: Parametric IRF — both samples, quadratic & cubic
# ==================================================================
cat("\n=== Part B: Parametric IRF (stacked LP) ===\n")

stacked_samples <- lapply(samples, build_stacked)

res_B <- list()
for (nm in names(samples)) {
  cat(sprintf("\nSample: %s (%d rows)\n", nm,  nrow(stacked_samples[[nm]])))
  for (deg in c(2, 3)) {
    key <- paste(nm, deg)
    res_B[[key]] <- estimate_parametric(deg, stacked_samples[[nm]], label = nm)
    res_B[[key]]$pred$sample <- nm
  }
}

pred_irf_df <- bind_rows(lapply(res_B, `[[`, "pred"))

# ==================================================================
# PART C: Country-specific quadratic IRFs (pooled, treated sample)
# ==================================================================
cat("\n=== Part C: Country-specific quadratic IRFs (pooled) ===\n")

# Use same-continent stacked data (more countries for year FE identification)
stacked_cont <- stacked_samples[["Same continent"]]

country_shocks <- d_continent %>%
  filter(maxwind_95 == 1) %>%
  count(countrycode, name = "n_shocks") %>%
  filter(n_shocks >= 5) %>%
  arrange(desc(n_shocks))

cat(sprintf("Countries with >= 5 shocks: %d\n", nrow(country_shocks)))

stacked_sub <- stacked_cont %>%
  filter(countrycode %in% country_shocks$countrycode) %>%
  mutate(cc = factor(countrycode))

mod_country <- tryCatch({
  feols(
    dy ~ i(cc, shock_h0) + i(cc, shock_h1) + i(cc, shock_h2) +
      i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2) +
      i(fh, shock_l1) + i(fh, shock_l2) |
      country_h[year] + country_h[year2] + year^fh,
    data = stacked_sub,
    vcov = ~countrycode
  )
}, error = function(e) {
  cat("  Note: Trend FE failed, trying simpler...\n")
  feols(
    dy ~ i(cc, shock_h0) + i(cc, shock_h1) + i(cc, shock_h2) +
      i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2) +
      i(fh, shock_l1) + i(fh, shock_l2) |
      countrycode^fh + year^fh,
    data = stacked_sub,
    vcov = ~countrycode
  )
})

all_coefs <- coef(mod_country)
all_vcov  <- vcov(mod_country)
h_seq <- seq(0, H, by = 0.25)
country_irf_list <- list()

for (cc in country_shocks$countrycode) {
  n_shk <- country_shocks$n_shocks[country_shocks$countrycode == cc]
  sv <- paste0("cc::", cc, ":shock_h", 0:2)
  available <- sv[sv %in% names(all_coefs)]

  if (length(available) < 3) {
    cat(sprintf("  %s: only %d of 3 coefficients, skipping\n", cc, length(available)))
    next
  }

  theta_cc <- all_coefs[sv]
  V_cc     <- all_vcov[sv, sv]
  if (any(is.na(theta_cc))) next

  for (hh in h_seq) {
    x <- c(1, hh, hh^2)
    irf_val <- sum(x * theta_cc)
    irf_se  <- as.numeric(sqrt(t(x) %*% V_cc %*% x))
    country_irf_list[[length(country_irf_list) + 1]] <- data.frame(
      countrycode = cc, horizon = hh, irf_mean = irf_val, se = irf_se,
      irf_down = irf_val - 1.96 * irf_se, irf_up = irf_val + 1.96 * irf_se,
      n_shocks = n_shk
    )
  }
  cat(sprintf("  %s (%d shocks): IRF at h=5 = %.3f\n",
              cc, n_shk, sum(c(1, 5, 25) * theta_cc)))
}

country_irf_df <- bind_rows(country_irf_list) %>%
  mutate(
    country_name = countrycode(countrycode, "iso3c", "country.name"),
    label = paste0(country_name, " (n=", n_shocks, ")"),
    label = fct_reorder(label, -n_shocks)
  )

cat(sprintf("\nSuccessfully estimated for %d countries\n",
            n_distinct(country_irf_df$countrycode)))

# ==================================================================
# FIGURES
# ==================================================================
dir.create("figures", showWarnings = FALSE)

irf_theme <- theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "bottom"
  )

gap_colors <- c("1"  = "#d73027", "3"  = "#fc8d59", "5"  = "#fee08b",
                "10" = "#91bfdb", "15" = "#4575b4")
gap_labels <- c("1" = "Gap = 1 yr (compound)", "3" = "Gap = 3 yrs",
                "5" = "Gap = 5 yrs", "10" = "Gap = 10 yrs",
                "15" = "Gap = 15 yrs (isolated)")

# ------------------------------------------------------------------
# Figure A: Standard LP — treated vs same-continent, faceted
# ------------------------------------------------------------------
p_gap_std <- ggplot(irf_gap_df %>% mutate(gap = factor(gap)),
                     aes(x = horizon, y = irf_mean, color = gap, fill = gap)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.10, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~sample, ncol = 2) +
  scale_color_manual(values = gap_colors, labels = gap_labels) +
  scale_fill_manual(values = gap_colors, labels = gap_labels) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 0:10) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Shock Effect by Years Since Previous Shock",
    subtitle = "Standard LP: shock + shock x gap | Year FE, country trends, DK SEs",
    color = NULL, fill = NULL
  ) +
  irf_theme

ggsave("figures/irf_gap_standard.png", p_gap_std, width = 13, height = 5.5, dpi = 300)
cat("\nSaved: figures/irf_gap_standard.png\n")

# ------------------------------------------------------------------
# Figure B: Gap coefficient — treated vs same-continent
# ------------------------------------------------------------------
gap_coef_df <- gap_coef_df %>%
  mutate(ci_lo = beta_gap - 1.96 * se_gap, ci_hi = beta_gap + 1.96 * se_gap)

p_gap_coef <- ggplot(gap_coef_df, aes(x = horizon, y = beta_gap)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, fill = "steelblue") +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  facet_wrap(~sample, ncol = 2) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 0:10) +
  labs(
    x = "Horizon (years)",
    y = expression(beta[gap] ~ "(per additional year)"),
    title = "Gap Interaction Coefficient Across Horizons",
    subtitle = "Coefficient on shock x gap | Negative = longer gaps more damaging"
  ) +
  irf_theme

ggsave("figures/irf_gap_coefficient.png", p_gap_coef, width = 12, height = 5, dpi = 300)
cat("Saved: figures/irf_gap_coefficient.png\n")

# ------------------------------------------------------------------
# Figure C: Parametric IRF — quadratic, treated vs same-continent
# ------------------------------------------------------------------
pred_quad <- pred_irf_df %>% filter(degree == 2)

p_gap_param <- ggplot(pred_quad %>% mutate(gap = factor(gap)),
                       aes(x = horizon, y = irf_mean, color = gap, fill = gap)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.10, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~sample, ncol = 2) +
  scale_color_manual(values = gap_colors, labels = gap_labels) +
  scale_fill_manual(values = gap_colors, labels = gap_labels) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 0:10) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Parametric IRF (Quadratic in h): Treated vs Same-Continent Controls",
    subtitle = "Stacked LP with gap interaction | Country x horizon trends, year x horizon FE",
    color = NULL, fill = NULL
  ) +
  irf_theme

ggsave("figures/irf_gap_parametric.png", p_gap_param, width = 13, height = 5.5, dpi = 300)
cat("Saved: figures/irf_gap_parametric.png\n")

# ------------------------------------------------------------------
# Figure D: Quadratic vs Cubic — same-continent sample only
# ------------------------------------------------------------------
degree_labels <- c("2" = "Quadratic in h", "3" = "Cubic in h")

pred_cont <- pred_irf_df %>% filter(sample == "Same continent")

p_deg_compare <- ggplot(pred_cont %>% mutate(gap = factor(gap), degree = factor(degree)),
                         aes(x = horizon, y = irf_mean, color = gap, fill = gap)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.10, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~degree, ncol = 2, labeller = labeller(degree = degree_labels)) +
  scale_color_manual(values = gap_colors, labels = gap_labels) +
  scale_fill_manual(values = gap_colors, labels = gap_labels) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 0:10) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Quadratic vs Cubic Polynomial (Same-Continent Sample)",
    subtitle = "Stacked LP with gap interaction | Country trends, clustered by country",
    color = NULL, fill = NULL
  ) +
  irf_theme

ggsave("figures/irf_gap_degree_compare.png", p_deg_compare, width = 13, height = 5.5, dpi = 300)
cat("Saved: figures/irf_gap_degree_compare.png\n")

# ------------------------------------------------------------------
# Figure E: Country-specific quadratic IRFs, faceted
# ------------------------------------------------------------------
n_countries <- n_distinct(country_irf_df$countrycode)
ncol_facet  <- min(5, n_countries)
nrow_facet  <- ceiling(n_countries / ncol_facet)

p_country <- ggplot(country_irf_df, aes(x = horizon, y = irf_mean)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.2, fill = "steelblue") +
  geom_line(linewidth = 0.8, color = "steelblue") +
  facet_wrap(~label, scales = "free_y", ncol = ncol_facet) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = seq(0, 10, 5)) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Country-Specific Parametric IRFs (Quadratic in h)",
    subtitle = "Pooled regression with country-interacted shock polynomial | Same-continent sample"
  ) +
  irf_theme +
  theme(strip.text = element_text(size = 8))

fig_height <- max(4, nrow_facet * 2.5)
ggsave("figures/irf_country_parametric.png", p_country,
       width = 14, height = fig_height, dpi = 300)
cat("Saved: figures/irf_country_parametric.png\n")

# ==================================================================
# PART D: Parametric IRF across treatment thresholds
# ==================================================================
cat("\n=== Part D: Parametric IRF across thresholds ===\n")

thresholds <- c(
  "p95 (60 kts)"    = 60 * knot_to_ms,
  "Cat 1+ (64 kts)" = 64 * knot_to_ms,
  "Cat 2+ (83 kts)" = 83 * knot_to_ms,
  "Cat 3+ (96 kts)" = 96 * knot_to_ms
)

# We already have d_full (the raw panel before threshold is applied).
# For each threshold: redefine the shock, rebuild same-continent sample,
# build stacked data, estimate quadratic parametric model.

thresh_pred_list <- list()

for (tname in names(thresholds)) {
  tval <- thresholds[tname]
  cat(sprintf("\n--- %s (%.1f m/s) ---\n", tname, tval))

  # Redefine shock at this threshold
  d_t <- d_full %>%
    mutate(maxwind_95 = as.integer(maxwind >= tval))

  n_shocks <- sum(d_t$maxwind_95)
  treated_cc <- d_t %>% filter(maxwind_95 == 1) %>% pull(countrycode) %>% unique()
  treated_cont <- d_t %>% filter(countrycode %in% treated_cc) %>%
    pull(continent) %>% unique()

  cat(sprintf("  Shocks: %d, Treated countries: %d\n", n_shocks, length(treated_cc)))

  # Same-continent sample
  d_samp <- d_t %>% filter(continent %in% treated_cont)

  # Add gap variables
  d_samp <- add_gap_vars(d_samp)

  # Build stacked data
  stacked_t <- build_stacked(d_samp)

  # Estimate quadratic parametric model
  res_t <- tryCatch(
    estimate_parametric(2, stacked_t, label = tname),
    error = function(e) {
      cat("  FAILED: ", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(res_t)) {
    thresh_pred_list[[tname]] <- res_t$pred %>%
      mutate(threshold = tname, n_shocks = n_shocks,
             n_treated = length(treated_cc))
  }
}

thresh_pred_df <- bind_rows(thresh_pred_list) %>%
  mutate(
    threshold = factor(threshold, levels = names(thresholds)),
    thresh_label = paste0(threshold, "\n(", n_shocks, " shocks, ", n_treated, " countries)")
  )

# Reorder labels to match threshold order
thresh_pred_df <- thresh_pred_df %>%
  mutate(thresh_label = fct_reorder(thresh_label, as.numeric(threshold)))

# Figure F: Parametric IRF faceted by threshold
p_thresh <- ggplot(thresh_pred_df %>% mutate(gap = factor(gap)),
                    aes(x = horizon, y = irf_mean, color = gap, fill = gap)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.10, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~thresh_label, ncol = 4) +
  scale_color_manual(values = gap_colors, labels = gap_labels) +
  scale_fill_manual(values = gap_colors, labels = gap_labels) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = seq(0, 10, 5)) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Parametric IRF Across Treatment Thresholds",
    subtitle = "Quadratic in h, same-continent sample | Country x horizon trends, year x horizon FE",
    color = NULL, fill = NULL
  ) +
  irf_theme

ggsave("figures/irf_gap_thresholds.png", p_thresh, width = 15, height = 5, dpi = 300)
cat("Saved: figures/irf_gap_thresholds.png\n")

# ------------------------------------------------------------------
# Save all results
# ------------------------------------------------------------------
saveRDS(list(
  irf_gap_standard   = irf_gap_df,
  irf_gap_parametric = pred_irf_df,
  gap_coefs          = gap_coef_df,
  country_irfs       = country_irf_df,
  thresh_parametric  = thresh_pred_df,
  treated_continents = treated_continents,
  n_continent_countries = n_distinct(d_continent$countrycode)
), "smooth_results.rds")

cat("\nSaved: smooth_results.rds\nDone.\n")
