# strategy2_event_stack.R
# ------------------------------------------------------------------
# Strategy 2: event-stacked design using FORWARD time-to-next-shock
# as the treatment heterogeneity variable.
#
# Identification: conditional on country (which absorbs the Poisson
# hazard lambda_i), the time until the next Cat 1+ shock is an
# exponential draw, exogenous to the economic state at the index
# shock. This side-steps the selection / state-dependence problems
# with backward-looking `gap since last shock'.
#
# Specs (all on same-continent sample, stacked LP with non-parametric
# horizon dummies for both beta(h) and the interaction):
#   (E0) Baseline: average IRF, no next-shock interaction
#   (E1) Continuous: shock x time_to_next, non-parametric phi(h)
#   (E2) Binary compound (time_to_next <= 3), non-parametric psi(h)
#   (E3) Binary compound, scalar psi (horizon-invariant)
#
# Outputs:
#   figures/strategy2/*.png
#   strategy2_results.rds
#   strategy2_diagnostics.txt
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(haven)
  library(countrycode)
})

setFixest_notes(FALSE)

out_fig <- "figures/strategy2"
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

# ==================================================================
# 1. DATA (mirrors compound_smooth.R / strategy1 for comparability)
# ==================================================================
knot_to_ms  <- 0.514444
# Baseline threshold: Cat 1+ (64 knots), selected from the power grid.
cat2_thresh <- 64 * knot_to_ms
H           <- 10
COMPOUND_K  <- 3   # binary compound cutoff: time_to_next <= K
CENSOR_CAP  <- 15  # cap time_to_next at this

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

# ==================================================================
# 2. FORWARD TIME-TO-NEXT-SHOCK
# ==================================================================
# For each country-year (t, i), time_to_next_it = years until the
# next Cat 1+ shock strictly after year t in country i.
# NA if no further shock in the sample (right-censored).
d_samp <- d_samp %>%
  group_by(countrycode) %>%
  arrange(year) %>%
  mutate(
    time_to_next = {
      sy <- year[shock == 1]
      vapply(year, function(y) {
        fut <- sy[sy > y]
        if (length(fut) == 0) NA_real_ else min(fut) - y
      }, numeric(1))
    }
  ) %>%
  ungroup() %>%
  mutate(
    # For E1 (continuous) we cap time_to_next at CENSOR_CAP and treat
    # censored events as "isolated" (diagnostics report the censored
    # count). For E2/E3 (binary) we need the K-year forward window to
    # be observed -- otherwise compound_next is genuinely unknown.
    ttn_censored       = as.integer(is.na(time_to_next)),
    time_to_next_cap   = pmin(replace_na(time_to_next, CENSOR_CAP), CENSOR_CAP),
    compound_next      = as.integer(!is.na(time_to_next) & time_to_next <= COMPOUND_K),
    # window_obs: event year t_e is far enough from sample end that
    # the K-year compound window is observable.
    window_obs         = (year + COMPOUND_K) <= max(year)
  )

n_events         <- sum(d_samp$shock == 1)
n_events_compnd  <- sum(d_samp$shock == 1 & d_samp$compound_next == 1)
n_events_censrd  <- sum(d_samp$shock == 1 & d_samp$ttn_censored == 1)
n_events_iso     <- n_events - n_events_compnd
cat(sprintf("Events (shock-years): %d\n", n_events))
cat(sprintf("  compound (next shock within %d yrs): %d (%.1f%%)\n",
            COMPOUND_K, n_events_compnd, 100 * n_events_compnd / n_events))
cat(sprintf("  isolated / censored:                 %d (%.1f%%)\n",
            n_events_iso,  100 * n_events_iso  / n_events))
cat(sprintf("  right-censored (no future shock in sample): %d\n", n_events_censrd))

# ==================================================================
# 3. EDA / DIAGNOSTICS
# ==================================================================
theme_set(theme_classic(base_size = 12))

shock_rows <- d_samp %>% filter(shock == 1)

## 3a. Distribution of time_to_next across events (pre-censoring)
p_ttn <- ggplot(shock_rows %>% filter(!is.na(time_to_next)),
                aes(time_to_next)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  geom_vline(xintercept = COMPOUND_K + 0.5, linetype = "dashed",
             color = "red") +
  annotate("text", x = COMPOUND_K + 1, y = Inf, vjust = 1.5,
           label = sprintf("compound cutoff\n(<= %d yrs)", COMPOUND_K),
           hjust = 0, size = 3, color = "red") +
  labs(x = "Years until next Cat 1+ shock (same country)",
       y = "Events",
       title = "Distribution of forward time-to-next-shock",
       subtitle = sprintf("%d events, %d right-censored and dropped from histogram",
                          n_events, n_events_censrd))
ggsave(file.path(out_fig, "ttn_hist.png"), p_ttn,
       width = 8, height = 4.5, dpi = 200)

## 3b. Within- vs between-country variance of time_to_next_cap
var_dec <- shock_rows %>%
  group_by(countrycode) %>%
  mutate(ttn_dm = time_to_next_cap - mean(time_to_next_cap)) %>%
  ungroup() %>%
  summarise(
    var_total   = var(time_to_next_cap),
    var_within  = var(ttn_dm),
    var_between = var_total - var_within
  )
cat(sprintf("\nVar(time_to_next_cap) = %.2f\n", var_dec$var_total))
cat(sprintf("  within-country  = %.2f (%.1f%%)\n",
            var_dec$var_within, 100*var_dec$var_within/var_dec$var_total))
cat(sprintf("  between-country = %.2f (%.1f%%)\n",
            var_dec$var_between, 100*var_dec$var_between/var_dec$var_total))

## 3c. Compound-next rate by country (sanity)
by_cc <- shock_rows %>%
  group_by(countrycode) %>%
  summarise(
    n_events       = n(),
    compound_share = mean(compound_next),
    .groups        = "drop"
  ) %>%
  filter(n_events >= 3) %>%
  arrange(desc(n_events))

p_cc <- ggplot(by_cc, aes(reorder(countrycode, n_events), compound_share,
                          fill = n_events)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "D", name = "# events") +
  labs(x = NULL, y = "Share of events that are 'compound' (next shock within 3 yrs)",
       title = "Compound-next share by country",
       subtitle = "Countries with >= 3 Cat 1+ events")
ggsave(file.path(out_fig, "compound_by_country.png"), p_cc,
       width = 8, height = 6, dpi = 200)

# ==================================================================
# 4. STACKED LP
# ==================================================================
# Same structure as strategy1 / compound_smooth.R. The `shock` main
# effect is polynomial-in-h; the next-shock interaction layers on top.
build_stacked <- function(dat) {
  out <- vector("list", H + 1)
  for (h in 0:H) {
    out[[h + 1]] <- dat %>%
      group_by(countrycode) %>%
      arrange(year) %>%
      mutate(
        # Same LP normalization as compound_smooth.R. The identification
        # difference vs Strategy 1 is not the dy formula but the regressor:
        # time_to_next is forward-looking so y_{t-1} at an index shock is
        # the pre-shock level, not a level already depressed by a
        # sample-selected prior shock. No mechanical baseline bias.
        dy          = lead(loggdp, h) - lag(loggdp, 1),  # y_{t+h} - y_{t-1}
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
      shock_h0        = shock,
      # Continuous forward-time interaction
      shock_ttn_h0    = shock * time_to_next_cap,
      # Binary compound-next interaction
      shock_cmp_h0    = shock * compound_next,
      country_h       = interaction(countrycode, fh, drop = TRUE),
      region_h        = interaction(region,      fh, drop = TRUE)
    )
}

sd_stack <- build_stacked(d_samp)
cat(sprintf("\nStacked rows: %d\n", nrow(sd_stack)))

fe_base <- "country_h[year] + country_h[year2] + year^fh"

# --- Non-parametric estimator: horizon dummies for beta(h) and phi(h) ---
# int_var optional: if supplied, adds i(fh, int_var) and predicts IRFs on
# a grid of interaction values.
estimate_np <- function(stacked_df, main_var = "shock_h0", int_var = NULL,
                        fe_formula, grid = c(1, 3, 5, 10, 15),
                        label = "") {
  vars <- c(main_var)
  if (!is.null(int_var)) vars <- c(vars, int_var)
  rhs <- paste(sprintf("i(fh, %s)", vars), collapse = " + ")
  fml <- as.formula(paste0(
    "dy ~ ", rhs,
    " + i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2)",
    " + i(fh, shock_l1) + i(fh, shock_l2) | ", fe_formula
  ))
  mod <- feols(fml, data = stacked_df, vcov = ~countrycode)
  cf  <- coef(mod); Vm <- vcov(mod)

  pred <- list()
  # Anchor at h=-1: by construction dy = 0
  if (is.null(int_var)) {
    pred[[length(pred)+1]] <- data.frame(horizon = -1, grp = NA_character_,
      irf_mean = 0, se = 0, irf_down = 0, irf_up = 0)
  } else {
    for (g in grid) {
      pred[[length(pred)+1]] <- data.frame(horizon = -1, grp = as.character(g),
        irf_mean = 0, se = 0, irf_down = 0, irf_up = 0)
    }
  }

  for (h in 0:H) {
    idx <- sprintf("fh::%d:%s", h, vars)
    if (!all(idx %in% names(cf))) next
    th <- cf[idx]; V <- Vm[idx, idx, drop = FALSE]

    if (is.null(int_var)) {
      val <- th[1]; se <- sqrt(V[1, 1])
      pred[[length(pred)+1]] <- data.frame(horizon = h, grp = NA_character_,
        irf_mean = val, se = se,
        irf_down = val - 1.96*se, irf_up = val + 1.96*se)
    } else {
      for (g in grid) {
        x <- c(1, g)
        val <- sum(x * th); se <- as.numeric(sqrt(t(x) %*% V %*% x))
        pred[[length(pred)+1]] <- data.frame(horizon = h, grp = as.character(g),
          irf_mean = val, se = se,
          irf_down = val - 1.96*se, irf_up = val + 1.96*se)
      }
    }
  }
  list(pred = bind_rows(pred) %>% mutate(spec = label),
       model = mod, int_var = int_var)
}

cat("\n--- Spec E0: baseline (no next-shock interaction) ---\n")
e0 <- estimate_np(sd_stack, int_var = NULL,
                  fe_formula = fe_base,
                  label = "(E0) Baseline")

cat("--- Spec E1: shock x time_to_next (continuous) ---\n")
e1 <- estimate_np(sd_stack, int_var = "shock_ttn_h0",
                  fe_formula = fe_base,
                  grid = c(1, 3, 5, 10, 15),
                  label = "(E1) shock x time_to_next")

cat("--- Spec E2: shock x compound_next (binary) ---\n")
# Drop rows where the current event's K-year compound window is not
# observed (i.e. shock_it=1 with window_obs=FALSE). Non-shock rows
# are retained to identify FE.
sd_stack_e2 <- sd_stack %>% filter(shock == 0 | window_obs)
cat(sprintf("  Dropped %d stacked rows (events with unobservable K-year window)\n",
            nrow(sd_stack) - nrow(sd_stack_e2)))
e2 <- estimate_np(sd_stack_e2, int_var = "shock_cmp_h0",
                  fe_formula = fe_base,
                  grid = c(0, 1),
                  label = "(E2) shock x compound_next")

# --- Spec E3: scalar binary (no h dummies on the interaction) ---
cat("--- Spec E3: shock x compound_next, scalar psi ---\n")
fml_e3 <- as.formula(paste0(
  "dy ~ i(fh, shock_h0) + shock_cmp_h0",
  " + i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2)",
  " + i(fh, shock_l1) + i(fh, shock_l2) | ", fe_base
))
mod_e3 <- feols(fml_e3, data = sd_stack_e2, vcov = ~countrycode)
psi_hat <- coef(mod_e3)["shock_cmp_h0"]
psi_se  <- sqrt(vcov(mod_e3)["shock_cmp_h0", "shock_cmp_h0"])
cat(sprintf("  psi = %.4f  (SE %.4f) t = %.2f\n",
            psi_hat, psi_se, psi_hat / psi_se))

# ==================================================================
# 5. FIGURES
# ==================================================================

## 5a. Baseline IRF (spec E0), with pre-index normalization
p_base <- ggplot(e0$pred, aes(horizon, irf_mean)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up),
              alpha = 0.2, fill = "steelblue") +
  geom_line(color = "steelblue", linewidth = 1) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Horizon (years)",
       y = "Cumulative effect on log GDP per capita",
       title = "Spec E0: baseline IRF (no next-shock interaction)",
       subtitle = "Same-continent sample, pre-index normalization")
ggsave(file.path(out_fig, "irf_E0_baseline.png"), p_base,
       width = 7, height = 4.5, dpi = 200)

## 5b. Continuous: IRFs by time_to_next value (spec E1)
ttn_colors <- c("1" = "#d73027", "3" = "#fc8d59", "5" = "#fee08b",
                "10" = "#91bfdb", "15" = "#4575b4")
ttn_labels <- c("1" = "next shock in 1 yr (compound)",
                "3" = "next shock in 3 yrs",
                "5" = "next shock in 5 yrs",
                "10" = "next shock in 10 yrs",
                "15" = "no near shock (isolated)")

p_e1 <- ggplot(e1$pred %>% mutate(grp = factor(grp,
                                               levels = c("1","3","5","10","15"))),
               aes(horizon, irf_mean, color = grp, fill = grp)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up),
              alpha = 0.10, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = ttn_colors, labels = ttn_labels) +
  scale_fill_manual(values = ttn_colors, labels = ttn_labels) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Horizon (years)", y = "Cumulative effect on log GDP pc",
       color = NULL, fill = NULL,
       title = "Spec E1: IRFs by forward time-to-next-shock",
       subtitle = "Continuous shock x time_to_next interaction, non-parametric in h") +
  theme(legend.position = "bottom")
ggsave(file.path(out_fig, "irf_E1_continuous.png"), p_e1,
       width = 9, height = 5, dpi = 200)

## 5c. Binary: compound vs isolated trajectories (spec E2)
cmp_labels <- c("0" = "isolated (no next shock in 3 yrs)",
                "1" = "compound (next shock in 3 yrs)")
p_e2 <- ggplot(e2$pred %>% mutate(grp = factor(grp, levels = c("0","1"))),
               aes(horizon, irf_mean, color = grp, fill = grp)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("0" = "#4575b4", "1" = "#d73027"),
                     labels = cmp_labels) +
  scale_fill_manual(values = c("0" = "#4575b4", "1" = "#d73027"),
                    labels = cmp_labels) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Horizon (years)", y = "Cumulative effect on log GDP pc",
       color = NULL, fill = NULL,
       title = "Spec E2: IRFs by compound-next status",
       subtitle = sprintf("Binary shock x compound_next (cutoff: next shock within %d yrs)",
                          COMPOUND_K)) +
  theme(legend.position = "bottom")
ggsave(file.path(out_fig, "irf_E2_binary.png"), p_e2,
       width = 9, height = 5, dpi = 200)

## 5d. phi(h) curve: pull gap-interaction horizon dummies directly
phi_curve <- function(res) {
  if (is.null(res$int_var)) return(NULL)
  cf <- coef(res$model); V <- vcov(res$model)
  rows <- list()
  for (h in 0:H) {
    nm <- sprintf("fh::%d:%s", h, res$int_var)
    if (!(nm %in% names(cf))) next
    rows[[length(rows)+1]] <- data.frame(
      horizon = h, phi = cf[nm], se = sqrt(V[nm, nm])
    )
  }
  bind_rows(rows) %>% mutate(lo = phi - 1.96*se, hi = phi + 1.96*se)
}
phi_e1 <- phi_curve(e1) %>% mutate(spec = "(E1) time_to_next")
psi_e2 <- phi_curve(e2) %>% mutate(spec = "(E2) compound_next")

p_phi <- bind_rows(phi_e1, psi_e2) %>%
  ggplot(aes(horizon, phi)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, color = "steelblue") +
  geom_line(color = "steelblue", linewidth = 0.6) +
  geom_point(color = "steelblue", size = 2) +
  facet_wrap(~spec, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = 0:H) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "Horizon h (years)", y = expression(phi(h) ~ "or " * psi(h)),
       title = "Interaction slope across horizons (non-parametric)",
       subtitle = "Spec E1: per year of time_to_next. Spec E2: compound - isolated.")
ggsave(file.path(out_fig, "phi_curve_E1_E2.png"), p_phi,
       width = 11, height = 4.5, dpi = 200)

# ==================================================================
# 6. DIAGNOSTICS FILE
# ==================================================================
phi5 <- function(res, label) {
  if (is.null(res$int_var)) return(NULL)
  nm <- sprintf("fh::5:%s", res$int_var)
  cf <- coef(res$model); V <- vcov(res$model)
  if (!(nm %in% names(cf))) return(NULL)
  val <- cf[nm]; se <- sqrt(V[nm, nm])
  data.frame(spec = label, val = val, se = se,
             lo = val - 1.96*se, hi = val + 1.96*se,
             row.names = NULL)
}

phi5_df <- bind_rows(
  phi5(e1, "(E1) phi(h=5): time_to_next"),
  phi5(e2, "(E2) psi(h=5): compound_next"),
  data.frame(spec = "(E3) scalar psi", val = psi_hat, se = psi_se,
             lo = psi_hat - 1.96*psi_se, hi = psi_hat + 1.96*psi_se,
             row.names = NULL)
)

diag <- c(
  "Strategy 2 diagnostics (event-stacked, forward time-to-next)",
  "=============================================================",
  sprintf("Sample: %d countries, %d obs, %d Cat 1+ shocks (events)",
          n_distinct(d_samp$countrycode), nrow(d_samp), n_events),
  sprintf("Compound events (next shock within %d yrs): %d (%.1f%%)",
          COMPOUND_K, n_events_compnd, 100*n_events_compnd/n_events),
  sprintf("Right-censored events (no future shock in sample): %d (%.1f%%)",
          n_events_censrd, 100*n_events_censrd/n_events),
  "",
  "Variance decomposition of time_to_next_cap (over events):",
  sprintf("  Var(ttn) = %.2f", var_dec$var_total),
  sprintf("  within   = %.2f (%.1f%%)",
          var_dec$var_within, 100*var_dec$var_within/var_dec$var_total),
  sprintf("  between  = %.2f (%.1f%%)",
          var_dec$var_between, 100*var_dec$var_between/var_dec$var_total),
  "",
  "Interaction slope at h=5 (non-parametric specs) and scalar psi (E3):",
  paste0("  ", capture.output(print(phi5_df, row.names = FALSE)))
)
writeLines(diag, "strategy2_diagnostics.txt")
cat("\nSaved: strategy2_diagnostics.txt\n")

# ==================================================================
# 7. SAVE RESULTS
# ==================================================================
saveRDS(list(
  e0_pred  = e0$pred,
  e1_pred  = e1$pred,
  e2_pred  = e2$pred,
  e3_psi   = psi_hat,
  e3_se    = psi_se,
  var_dec  = var_dec,
  phi5_df  = phi5_df,
  n_events = n_events,
  n_compnd = n_events_compnd,
  n_censrd = n_events_censrd
), "strategy2_results.rds")

cat("Saved: strategy2_results.rds\nDone.\n")
