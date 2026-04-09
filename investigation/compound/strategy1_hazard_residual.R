# strategy1_hazard_residual.R
# ------------------------------------------------------------------
# Strategy 1: residualize the shock x gap interaction against the
# country-specific cyclone hazard rate, to isolate the within-country
# (Poisson arrival noise) component of realized inter-arrival times.
#
# Runs four LP specs on the same-continent sample:
#   (0) Baseline:        shock only (no gap interaction)      -- context
#   (1) Baseline w/ gap: shock + shock x gap                  -- replication
#   (2) Residualized:    shock + shock x gap_resid            -- main Strategy 1
#   (3) Lazy FWL check:  shock + shock x gap + shock x lambda -- equivalent to (2)
#   (4) Region x year:   as (2) but with region x year x h FE -- clean controls
#
# Plus EDA diagnostics for the Poisson-arrival assumption.
#
# Outputs:
#   figures/strategy1/*.png
#   strategy1_results.rds
#   strategy1_diagnostics.txt
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(haven)
  library(countrycode)
  library(MASS, include.only = "glm.nb")
})

setFixest_notes(FALSE)

out_fig <- "figures/strategy1"
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

# ==================================================================
# 1. DATA (mirrors compound_smooth.R so results are comparable)
# ==================================================================
knot_to_ms  <- 0.514444
# Baseline threshold: Cat 1+ (64 knots). Lower threshold chosen after
# the baseline-power grid showed Cat 1+ has the highest t-stat on
# beta(h=5) because the 374-event sample dominates Cat 2+'s 191.
cat2_thresh <- 64 * knot_to_ms

pwt <- read_stata("../../raw_data/pwt_clean.dta")
tcs <- read_stata("../../raw_data/ibtracs_clean.dta")

d_full <- pwt %>%
  left_join(tcs, by = c("year", "countrycode")) %>%
  mutate(
    maxwind   = replace_na(max_ann_wind_i_nob, 0) * knot_to_ms,
    n_landf   = replace_na(sum_a_lands_nob, 0),
    loggdp    = ln_real_gdp_usd_pc,
    gdp_diff  = growth_real_gdp_usd_pc,
    year2     = year^2
  ) %>%
  filter(year >= 1969) %>%
  mutate(
    shock     = as.integer(maxwind >= cat2_thresh),
    continent = countrycode(countrycode, "iso3c", "continent"),
    continent = ifelse(is.na(continent), "Other", continent),
    region    = countrycode(countrycode, "iso3c", "region"),
    region    = ifelse(is.na(region), "Other", region)
  )

treated_cc   <- d_full %>% filter(shock == 1) %>% pull(countrycode) %>% unique()
treated_cont <- d_full %>% filter(countrycode %in% treated_cc) %>%
  pull(continent) %>% unique()

d_samp <- d_full %>% filter(continent %in% treated_cont)

cat(sprintf("Sample: %d countries, %d obs, %d Cat 1+ shocks\n",
            n_distinct(d_samp$countrycode), nrow(d_samp), sum(d_samp$shock)))

# --- Add gap (years since last shock) ---
add_gap <- function(dat) {
  gaps <- dat %>%
    filter(shock == 1) %>%
    group_by(countrycode) %>%
    arrange(year) %>%
    mutate(gap = year - lag(year)) %>%
    ungroup() %>%
    dplyr::select(countrycode, year, gap)
  firsty <- dat %>% group_by(countrycode) %>%
    summarise(first_year = min(year), .groups = "drop")
  gaps <- gaps %>% left_join(firsty, by = "countrycode") %>%
    mutate(gap = ifelse(is.na(gap), year - first_year, gap)) %>%
    dplyr::select(countrycode, year, gap)
  dat %>% left_join(gaps, by = c("countrycode", "year")) %>%
    mutate(
      gap        = replace_na(gap, 0),
      gap        = pmax(gap, 1),
      gap_capped = pmin(gap, 20)
    )
}
d_samp <- add_gap(d_samp)

# ==================================================================
# 2. HAZARD ESTIMATION
# ==================================================================
# Two hazards per country:
#   lambda_cat2 : mean annual count of Cat 2+ shock-years
#   lambda_all  : mean annual count of all landfalling storms (sub-threshold too)
hazard_df <- d_samp %>%
  group_by(countrycode) %>%
  summarise(
    nyrs         = n(),
    n_shocks     = sum(shock),
    lambda_cat2  = mean(shock),
    lambda_all   = mean(n_landf > 0),       # annual prob of any landfall
    mean_landf   = mean(n_landf),
    .groups      = "drop"
  ) %>%
  mutate(
    egap_cat2 = ifelse(lambda_cat2 > 0, 1 / lambda_cat2, NA_real_),
    egap_all  = ifelse(lambda_all  > 0, 1 / lambda_all,  NA_real_)
  )

d_samp <- d_samp %>%
  left_join(hazard_df %>%
              dplyr::select(countrycode, lambda_cat2, lambda_all,
                            egap_cat2, egap_all, n_shocks),
            by = "countrycode") %>%
  mutate(
    # headline residual: realized gap - expected gap (under Poisson(lambda_cat2))
    gap_resid     = ifelse(shock == 1, gap_capped - egap_cat2, 0),
    gap_resid     = replace_na(gap_resid, 0),
    gap_resid_cap = pmin(pmax(gap_resid, -10), 10),
    # FWL control for the "lazy" spec: interacting shock with 1/lambda.
    # Adding shock x (1/lambda) alongside shock x gap makes the coefficient
    # on shock x gap equal to the spec (2) residualized estimate.
    egap_ctrl     = ifelse(is.na(egap_cat2), 0, egap_cat2)
  )

# Sample mean of 1/lambda across shock-years. Used for IRF prediction
# in specs that use gap_resid (so the plot shows "raw gap = g", not
# "gap_resid = g") and for specs that include shock x (1/lambda) as a
# control (so the control is evaluated at its sample mean, not zero).
egap_mean <- mean(d_samp$egap_cat2[d_samp$shock == 1], na.rm = TRUE)
cat(sprintf("egap_mean (over shock-years) = %.3f years\n", egap_mean))

# ==================================================================
# 3. EDA / DIAGNOSTICS
# ==================================================================
theme_set(theme_classic(base_size = 12))

## 3a. Hazard histogram (Cat 1+ vs full landfall)
haz_long <- hazard_df %>%
  filter(n_shocks > 0) %>%
  dplyr::select(countrycode, lambda_cat2, lambda_all) %>%
  pivot_longer(-countrycode, names_to = "rate", values_to = "lambda") %>%
  mutate(rate = recode(rate,
                       lambda_cat2 = "Cat 1+ shock rate",
                       lambda_all  = "Any-landfall rate"))
p1 <- ggplot(haz_long, aes(lambda)) +
  geom_histogram(bins = 25, fill = "steelblue", color = "white") +
  facet_wrap(~rate, scales = "free") +
  labs(x = expression(hat(lambda)[i]), y = "Countries",
       title = "Country-specific cyclone hazard estimates")
ggsave(file.path(out_fig, "hazard_hist.png"), p1, width = 9, height = 4, dpi = 200)

## 3b. Realized gap vs 1/lambda (between-country confound)
shocks_only <- d_samp %>% filter(shock == 1)
p2 <- ggplot(shocks_only, aes(egap_cat2, gap_capped)) +
  geom_point(alpha = 0.4) +
  geom_abline(linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
  labs(x = expression("Expected gap " * 1/hat(lambda)[i]),
       y = "Realized gap (years, capped at 20)",
       title = "Cross-country confound: realized gap vs expected gap",
       subtitle = "Each point is a shock-year. Slope ~ share of variation that is between-country.")
ggsave(file.path(out_fig, "gap_vs_hazard.png"), p2, width = 7, height = 5, dpi = 200)

## 3c. Poisson assumption: QQ of inter-arrival times for countries with >=5 shocks
ia_df <- d_samp %>%
  filter(shock == 1, n_shocks >= 5) %>%
  group_by(countrycode) %>%
  arrange(year) %>%
  mutate(ia = year - lag(year)) %>%
  ungroup() %>%
  filter(!is.na(ia))
p3 <- ia_df %>%
  group_by(countrycode) %>%
  arrange(ia) %>%
  mutate(
    u_emp  = (seq_along(ia) - 0.5) / n(),
    u_theo = pexp(ia, rate = first(lambda_cat2))
  ) %>%
  ggplot(aes(u_theo, u_emp)) +
  geom_abline(linetype = "dashed", color = "red") +
  geom_point(alpha = 0.6) +
  facet_wrap(~countrycode) +
  coord_equal() +
  labs(x = "Exp(lambda) theoretical CDF",
       y = "Empirical CDF",
       title = "Poisson-arrival check: QQ of inter-arrival times",
       subtitle = "Countries with >= 5 Cat 1+ shocks")
ggsave(file.path(out_fig, "poisson_qq.png"), p3,
       width = 10, height = 7, dpi = 200)

## 3d. Overdispersion test: Poisson vs NegBin on country-year shock counts
od_test <- tryCatch({
  m_pois <- glm(shock ~ factor(countrycode), family = "poisson", data = d_samp)
  m_nb   <- glm.nb(shock ~ factor(countrycode), data = d_samp)
  list(
    dispersion = sum(residuals(m_pois, "pearson")^2) / df.residual(m_pois),
    theta      = m_nb$theta,
    lr_stat    = 2 * (logLik(m_nb) - logLik(m_pois))[1]
  )
}, error = function(e) list(dispersion = NA, theta = NA, lr_stat = NA))

## 3e. Variance decomposition: between vs within share of Var(gap|shock)
var_dec <- shocks_only %>%
  group_by(countrycode) %>%
  mutate(gap_dm = gap_capped - mean(gap_capped)) %>%
  ungroup() %>%
  summarise(
    var_total   = var(gap_capped),
    var_within  = var(gap_dm),
    var_between = var_total - var_within
  )

## 3f. Residualized gap vs pre-period state (channels 2-3 check)
p4 <- shocks_only %>%
  mutate(gap_resid_chk = gap_capped - egap_cat2) %>%
  ggplot(aes(gap_resid_chk, gdp_diff)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", color = "steelblue") +
  labs(x = "gap_resid (realized - expected)",
       y = "Pre-shock GDP growth",
       title = "Channel check: does gap_resid correlate with pre-state?",
       subtitle = "Flat line = Strategy 1 isolates clean noise; sloped = residual confound")
ggsave(file.path(out_fig, "gap_resid_vs_prestate.png"), p4,
       width = 7, height = 5, dpi = 200)

# ==================================================================
# 4. STACKED LP
# ==================================================================
H <- 10

build_stacked <- function(dat) {
  out <- vector("list", H + 1)
  for (h in 0:H) {
    out[[h + 1]] <- dat %>%
      group_by(countrycode) %>%
      arrange(year) %>%
      mutate(
        dy          = lead(loggdp, h) - lag(loggdp, 1),
        gdp_diff_l1 = lag(gdp_diff, 1),
        gdp_diff_l2 = lag(gdp_diff, 2),
        shock_l1    = lag(shock, 1),
        shock_l2    = lag(shock, 2)
      ) %>%
      ungroup() %>%
      filter(!is.na(dy), !is.na(gdp_diff_l1), !is.na(gdp_diff_l2)) %>%
      mutate(h = h)
  }
  bind_rows(out) %>%
    mutate(
      fh = factor(h),
      h1 = h, h2 = h^2,
      shock_h0        = shock,
      shock_h1        = shock * h1,
      shock_h2        = shock * h2,
      shock_gap_h0    = shock * gap_capped,
      shock_gap_h1    = shock * gap_capped * h1,
      shock_gap_h2    = shock * gap_capped * h2,
      shock_gapres_h0 = shock * gap_resid_cap,
      shock_gapres_h1 = shock * gap_resid_cap * h1,
      shock_gapres_h2 = shock * gap_resid_cap * h2,
      shock_lam_h0    = shock * egap_ctrl,
      shock_lam_h1    = shock * egap_ctrl * h1,
      shock_lam_h2    = shock * egap_ctrl * h2,
      country_h       = interaction(countrycode, fh, drop = TRUE),
      region_h        = interaction(region, fh, drop = TRUE)
    )
}

sd_stack <- build_stacked(d_samp)
cat(sprintf("Stacked rows: %d\n", nrow(sd_stack)))

# --- non-parametric estimator: horizon dummies for beta(h) and phi(h) ---
# Uses i(fh, x) so coefficients come out as fh::0:x, fh::1:x, ..., fh::H:x.
# gap_offset: subtract from g before multiplying (for residualized specs,
#             pass egap_mean so the plot reads as "raw gap = g")
# extra_mult: multiplier for the extra block (for specs that include
#             shock x (1/lambda) as a control, pass egap_mean so the
#             control is evaluated at its sample mean, not zero)
estimate_spec <- function(stacked_df, main_var, gap_var = NULL, extra_var = NULL,
                          fe_formula, gap_values = c(1, 3, 5, 10, 15),
                          gap_offset = 0, extra_mult = 0,
                          label = "") {
  int_vars <- c(main_var)
  if (!is.null(gap_var))   int_vars <- c(int_vars, gap_var)
  if (!is.null(extra_var)) int_vars <- c(int_vars, extra_var)

  rhs <- paste(sprintf("i(fh, %s)", int_vars), collapse = " + ")
  fml <- as.formula(paste0(
    "dy ~ ", rhs,
    " + i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2)",
    " + i(fh, shock_l1) + i(fh, shock_l2) | ", fe_formula
  ))
  mod <- feols(fml, data = stacked_df, vcov = ~countrycode)

  cf <- coef(mod)
  Vm <- vcov(mod)

  pred_list <- list()
  for (h in 0:H) {
    idx <- sprintf("fh::%d:%s", h, int_vars)
    # Skip horizons where a coefficient failed to identify
    if (!all(idx %in% names(cf))) next
    theta_h <- cf[idx]
    V_h     <- Vm[idx, idx, drop = FALSE]

    if (is.null(gap_var)) {
      # Spec 0: single beta(h)
      val <- theta_h[1]; s <- sqrt(V_h[1, 1])
      pred_list[[length(pred_list) + 1]] <- data.frame(
        horizon = h, gap = NA_real_, irf_mean = val, se = s,
        irf_down = val - 1.96 * s, irf_up = val + 1.96 * s
      )
    } else {
      for (g in gap_values) {
        # Linear combination weights
        x <- c(1, g - gap_offset)
        if (!is.null(extra_var)) x <- c(x, extra_mult)
        val <- sum(x * theta_h)
        s   <- as.numeric(sqrt(t(x) %*% V_h %*% x))
        pred_list[[length(pred_list) + 1]] <- data.frame(
          horizon = h, gap = g, irf_mean = val, se = s,
          irf_down = val - 1.96 * s, irf_up = val + 1.96 * s
        )
      }
    }
  }

  # Anchor: dy_{-1} = 0 by construction
  if (is.null(gap_var)) {
    pred_list <- c(list(data.frame(horizon = -1, gap = NA_real_,
                                   irf_mean = 0, se = 0,
                                   irf_down = 0, irf_up = 0)), pred_list)
  } else {
    for (g in gap_values) {
      pred_list <- c(list(data.frame(horizon = -1, gap = g,
                                     irf_mean = 0, se = 0,
                                     irf_down = 0, irf_up = 0)), pred_list)
    }
  }

  list(pred  = bind_rows(pred_list) %>% mutate(spec = label),
       model = mod, int_vars = int_vars)
}

# FE formula components
fe_base   <- "country_h[year] + country_h[year2] + year^fh"
fe_region <- "country_h[year] + country_h[year2] + region_h^year"

# Variable names (single interaction variable per block; horizon dummies
# are added inside estimate_spec via i(fh, var)).
main_var <- "shock_h0"        # shock
gap_var  <- "shock_gap_h0"    # shock * gap_capped
res_var  <- "shock_gapres_h0" # shock * gap_resid_cap
lam_var  <- "shock_lam_h0"    # shock * egap_ctrl (= 1/lambda)

cat("\n--- Spec 0: baseline (no gap interaction) ---\n")
r0 <- estimate_spec(sd_stack, main_var,
                    fe_formula = fe_base, label = "(0) Baseline: no gap")

cat("--- Spec 1: baseline + shock x gap ---\n")
r1 <- estimate_spec(sd_stack, main_var, gap_var = gap_var,
                    fe_formula = fe_base,
                    label = "(1) Baseline + gap")

cat("--- Spec 2: residualized gap (Strategy 1 main) ---\n")
r2 <- estimate_spec(sd_stack, main_var, gap_var = res_var,
                    fe_formula = fe_base,
                    gap_offset = egap_mean,
                    label = "(2) Residualized gap")

cat("--- Spec 3: FWL equivalent (gap + shock x 1/lambda) ---\n")
r3 <- estimate_spec(sd_stack, main_var, gap_var = gap_var, extra_var = lam_var,
                    fe_formula = fe_base,
                    extra_mult = egap_mean,
                    label = "(3) gap + shock x 1/lambda")

# ------------------------------------------------------------------
# Scalar-phi variants: gap interaction is a single scalar (constant
# across horizons). Main shock effect stays non-parametric via fh dummies.
# ------------------------------------------------------------------
fit_scalar <- function(stacked_df, gap_var, extra_vars = NULL,
                       fe_formula = fe_base, label = "") {
  # beta(h) via i(fh, shock); scalar phi on gap_var; optional scalar
  # control on extra_vars.
  rhs <- c("i(fh, shock_h0)", gap_var)
  if (!is.null(extra_vars)) rhs <- c(rhs, extra_vars)
  fml <- as.formula(paste0(
    "dy ~ ", paste(rhs, collapse = " + "),
    " + i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2)",
    " + i(fh, shock_l1) + i(fh, shock_l2) | ", fe_formula
  ))
  mod <- feols(fml, data = stacked_df, vcov = ~countrycode)
  b <- coef(mod)[gap_var]
  s <- sqrt(vcov(mod)[gap_var, gap_var])
  data.frame(spec = label, phi = b, se = s,
             lo = b - 1.96*s, hi = b + 1.96*s,
             row.names = NULL)
}

cat("--- Scalar spec 1s: replication, scalar phi ---\n")
s1 <- fit_scalar(sd_stack, "shock_gap_h0",    label = "(1s) Replication")
cat("--- Scalar spec 2s: residualized, scalar phi ---\n")
s2 <- fit_scalar(sd_stack, "shock_gapres_h0", label = "(2s) Residualized")
cat("--- Scalar spec 3s: control for 1/lambda, scalar phi ---\n")
s3 <- fit_scalar(sd_stack, "shock_gap_h0",
                 extra_vars = "shock_lam_h0",
                 label = "(3s) + 1/lambda control")
cat("--- Scalar spec 4s: residualized + region x year ---\n")
s4 <- fit_scalar(sd_stack, "shock_gapres_h0",
                 fe_formula = fe_region,
                 label = "(4s) + region x year")

phi_scalar_df <- bind_rows(s1, s2, s3, s4)
print(phi_scalar_df, row.names = FALSE)

p_phi_scalar <- ggplot(phi_scalar_df, aes(spec, phi)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_pointrange(aes(ymin = lo, ymax = hi), color = "steelblue",
                  size = 0.7) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  labs(x = NULL, y = expression("Scalar " * phi),
       title = expression("Scalar " * phi * ": gap effect constant across horizons"),
       subtitle = "Replaces phi(h) polynomial with a single parameter. ~sqrt(3) tighter SE.")
ggsave(file.path(out_fig, "phi_scalar_compare.png"), p_phi_scalar,
       width = 9, height = 4, dpi = 200)

cat("--- Spec 4: region x year FE, residualized gap ---\n")
r4 <- tryCatch(
  estimate_spec(sd_stack, main_var, gap_var = res_var,
                fe_formula = fe_region,
                gap_offset = egap_mean,
                label = "(4) Residualized + region x year"),
  error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); NULL }
)

all_pred <- bind_rows(
  r0$pred, r1$pred, r2$pred, r3$pred,
  if (!is.null(r4)) r4$pred else NULL
)

# ==================================================================
# 5. OUTPUT FIGURES
# ==================================================================
gap_colors <- c("1" = "#d73027", "3" = "#fc8d59", "5" = "#fee08b",
                "10" = "#91bfdb", "15" = "#4575b4")
gap_labels <- c("1" = "gap=1 (compound)", "3" = "gap=3", "5" = "gap=5",
                "10" = "gap=10", "15" = "gap=15 (isolated)")

## 5a. Baseline IRF on its own (context)
p_base <- ggplot(r0$pred, aes(horizon, irf_mean)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_vline(xintercept = -1, linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up),
              alpha = 0.2, fill = "steelblue") +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 2) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = c(-1, 0, 2, 4, 6, 8, 10)) +
  labs(x = "Horizon (years; h=-1 anchored at 0 by construction)",
       y = "Cumulative effect on log GDP per capita",
       title = "Baseline IRF: average Cat 1+ shock (no gap interaction)",
       subtitle = "Same-continent sample, non-parametric horizon dummies")
ggsave(file.path(out_fig, "irf_baseline.png"), p_base,
       width = 7, height = 4.5, dpi = 200)

## 5b. 2x2 facet across gap-interaction specs
gap_pred <- all_pred %>% filter(!is.na(gap))
p_cmp <- ggplot(gap_pred %>% mutate(gap = factor(gap)),
                aes(horizon, irf_mean, color = gap, fill = gap)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_vline(xintercept = -1, linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up),
              alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  facet_wrap(~spec, ncol = 4) +
  scale_color_manual(values = gap_colors, labels = gap_labels) +
  scale_fill_manual(values = gap_colors, labels = gap_labels) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = c(-1, 0, 2, 4, 6, 8, 10)) +
  labs(x = "Horizon (years; h=-1 anchored at 0 by construction)",
       y = "Cumulative effect on log GDP pc",
       color = NULL, fill = NULL,
       title = "Strategy 1: IRF under four specifications (non-parametric horizon dummies)",
       subtitle = sprintf("Gap values shown are RAW years. Residualized specs plot at gap_resid = g - E[1/lambda] (mean = %.2f yrs).",
                          egap_mean)) +
  theme(legend.position = "bottom")
ggsave(file.path(out_fig, "irf_compare_specs.png"), p_cmp,
       width = 18, height = 5.5, dpi = 200)

## 5c. phi(h) summary: extract the gap-interaction dummy coefficients
# directly from the non-parametric fit. The "gap var" is the 2nd block
# of int_vars (main is [1], gap is [2], optional extra is [3]).
phi_tab <- function(res, target_h = 5) {
  if (is.null(res) || length(res$int_vars) < 2) return(NULL)
  gvar  <- res$int_vars[2]
  nm    <- sprintf("fh::%d:%s", target_h, gvar)
  if (!(nm %in% names(coef(res$model)))) return(NULL)
  val   <- coef(res$model)[nm]
  s     <- sqrt(vcov(res$model)[nm, nm])
  data.frame(spec = unique(res$pred$spec),
             phi5 = val, se = s,
             lo = val - 1.96*s, hi = val + 1.96*s,
             row.names = NULL)
}
phi_df <- bind_rows(phi_tab(r1), phi_tab(r2), phi_tab(r3), phi_tab(r4)) %>%
  filter(!is.na(spec))

# --- phi(h) curve across horizons ---
phi_curve <- function(res) {
  if (is.null(res) || length(res$int_vars) < 2) return(NULL)
  gvar <- res$int_vars[2]
  cf <- coef(res$model); V <- vcov(res$model)
  rows <- list()
  for (h in 0:H) {
    nm <- sprintf("fh::%d:%s", h, gvar)
    if (!(nm %in% names(cf))) next
    rows[[length(rows)+1]] <- data.frame(
      horizon = h,
      phi     = cf[nm],
      se      = sqrt(V[nm, nm]),
      spec    = unique(res$pred$spec)
    )
  }
  bind_rows(rows)
}

phi_curve_df <- bind_rows(phi_curve(r1), phi_curve(r2),
                          phi_curve(r3), phi_curve(r4)) %>%
  mutate(lo = phi - 1.96 * se, hi = phi + 1.96 * se)

p_phi_curve <- ggplot(phi_curve_df, aes(horizon, phi)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, color = "steelblue") +
  geom_line(color = "steelblue", linewidth = 0.6) +
  geom_point(color = "steelblue", size = 2) +
  facet_wrap(~spec, ncol = 4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = 0:H) +
  labs(x = "Horizon h (years)",
       y = expression(phi(h) ~ "per additional year of gap"),
       title = expression("Gap interaction slope " * phi(h) * " across horizons"),
       subtitle = "Non-parametric horizon dummies. Positive = longer gap implies less damage.") +
  theme(strip.text = element_text(face = "bold"))
ggsave(file.path(out_fig, "phi_curve_compare.png"), p_phi_curve,
       width = 18, height = 5.5, dpi = 200)

p_phi <- ggplot(phi_df, aes(spec, phi5)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_pointrange(aes(ymin = lo, ymax = hi), color = "steelblue") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = expression(phi(h == 5)),
       title = "Gap-interaction slope at h=5 across specifications",
       subtitle = "Negative = longer gap -> more damage (compound-less-damaging gradient)")
ggsave(file.path(out_fig, "phi_h5_compare.png"), p_phi,
       width = 8, height = 3.5, dpi = 200)

# ==================================================================
# 6. DIAGNOSTICS FILE
# ==================================================================
diag_lines <- c(
  sprintf("Strategy 1 diagnostics"),
  sprintf("======================"),
  sprintf("Sample: %d countries, %d obs, %d Cat 1+ shocks",
          n_distinct(d_samp$countrycode), nrow(d_samp), sum(d_samp$shock)),
  "",
  sprintf("Hazard estimates (Cat 1+): mean=%.3f  median=%.3f  max=%.3f",
          mean(hazard_df$lambda_cat2), median(hazard_df$lambda_cat2),
          max(hazard_df$lambda_cat2)),
  sprintf("Countries with >=5 Cat 1+ shocks: %d",
          sum(hazard_df$n_shocks >= 5)),
  "",
  sprintf("Poisson overdispersion test (shock | country FE):"),
  sprintf("  dispersion ratio = %.3f (1 = Poisson)", od_test$dispersion),
  sprintf("  NB theta         = %.3f",               od_test$theta),
  sprintf("  LR NB vs Poisson = %.2f",               od_test$lr_stat),
  "",
  sprintf("Gap variance decomposition (shock-years only):"),
  sprintf("  Var(gap)         = %.2f", var_dec$var_total),
  sprintf("  within-country   = %.2f (%.1f%%)",
          var_dec$var_within, 100 * var_dec$var_within / var_dec$var_total),
  sprintf("  between-country  = %.2f (%.1f%%)",
          var_dec$var_between, 100 * var_dec$var_between / var_dec$var_total),
  "",
  sprintf("phi(h=5) across specs:"),
  paste0("  ", capture.output(print(phi_df, row.names = FALSE)))
)
writeLines(diag_lines, "strategy1_diagnostics.txt")
cat("\nSaved: strategy1_diagnostics.txt\n")

# ==================================================================
# 7. SAVE RESULTS
# ==================================================================
saveRDS(list(
  hazard_df = hazard_df,
  var_dec   = var_dec,
  od_test   = od_test,
  pred      = all_pred,
  phi_df    = phi_df,
  specs     = list(r0 = r0, r1 = r1, r2 = r2, r3 = r3, r4 = r4)
), "strategy1_results.rds")

cat("Saved: strategy1_results.rds\nDone.\n")
