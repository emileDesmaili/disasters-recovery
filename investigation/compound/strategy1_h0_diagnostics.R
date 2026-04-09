# strategy1_h0_diagnostics.R
# ------------------------------------------------------------------
# Diagnostics for the anomalous positive h=0 IRF for large gap in the
# residualized spec (strategy1_hazard_residual.R, spec 2).
#
# D1: raw binscatter of (y_t - y_{t-1}) vs gap, with and without
#     country demeaning. Asks: is the positive association in the
#     data, or a model artefact?
# D2: non-parametric horizon dummies (fh:shock, fh:shock:gap_resid)
#     replacing the quadratic polynomial. If h=0 stays weird, it's
#     data; if it flattens, it's a polynomial artefact.
# D3: control for pre-event state (pre-event 3-year growth rate,
#     interacted with shock). Tests channel 3 (state-dependence).
# D4: anchored parameterization with h' = h+1 so beta(-1) = 0 is
#     imposed mechanically. Tests endpoint pathology.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(haven)
  library(countrycode)
  library(patchwork)
})

setFixest_notes(FALSE)
out_fig <- "figures/strategy1_h0"
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

# ==================================================================
# 1. DATA (copied from strategy1_hazard_residual.R, Cat 1+)
# ==================================================================
knot_to_ms  <- 0.514444
cat1_thresh <- 64 * knot_to_ms
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
    shock     = as.integer(maxwind >= cat1_thresh),
    continent = countrycode(countrycode, "iso3c", "continent"),
    continent = ifelse(is.na(continent), "Other", continent),
    region    = countrycode(countrycode, "iso3c", "region"),
    region    = ifelse(is.na(region), "Other", region)
  )

treated_cc   <- d_full %>% filter(shock == 1) %>% pull(countrycode) %>% unique()
treated_cont <- d_full %>% filter(countrycode %in% treated_cc) %>%
  pull(continent) %>% unique()
d_samp <- d_full %>% filter(continent %in% treated_cont)

add_gap <- function(dat) {
  gaps <- dat %>%
    filter(shock == 1) %>%
    group_by(countrycode) %>%
    arrange(year) %>%
    mutate(gap = year - lag(year)) %>%
    ungroup() %>% dplyr::select(countrycode, year, gap)
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

hazard_df <- d_samp %>%
  group_by(countrycode) %>%
  summarise(
    n_shocks    = sum(shock),
    lambda_cat2 = mean(shock),
    .groups     = "drop"
  ) %>%
  mutate(egap_cat2 = ifelse(lambda_cat2 > 0, 1 / lambda_cat2, NA_real_))

d_samp <- d_samp %>%
  left_join(hazard_df, by = "countrycode") %>%
  mutate(
    gap_resid     = ifelse(shock == 1, gap_capped - egap_cat2, 0),
    gap_resid     = replace_na(gap_resid, 0),
    gap_resid_cap = pmin(pmax(gap_resid, -10), 10)
  )

egap_mean <- mean(d_samp$egap_cat2[d_samp$shock == 1], na.rm = TRUE)
cat(sprintf("Sample: %d obs, %d shocks. egap_mean = %.3f\n",
            nrow(d_samp), sum(d_samp$shock), egap_mean))

theme_set(theme_classic(base_size = 12))

# ==================================================================
# D1. RAW BINSCATTER
# ==================================================================
cat("\n=== D1: raw same-year GDP change vs gap ===\n")

dy0_df <- d_samp %>%
  group_by(countrycode) %>%
  arrange(year) %>%
  mutate(dy0 = loggdp - lag(loggdp, 1)) %>%
  ungroup() %>%
  filter(shock == 1, !is.na(dy0)) %>%
  group_by(countrycode) %>%
  mutate(
    dy0_dm = dy0 - mean(dy0),
    gap_dm = gap_capped - mean(gap_capped)
  ) %>%
  ungroup()

cat(sprintf("D1 sample: %d shock-years\n", nrow(dy0_df)))

# Slopes: raw, within-country demeaned, country-FE with raw gap,
# country-FE with residualized gap
d1_raw <- lm(dy0 ~ gap_capped, data = dy0_df)
d1_dm  <- lm(dy0_dm ~ gap_dm, data = dy0_df)
d1_fe  <- feols(dy0 ~ gap_capped | countrycode, data = dy0_df, vcov = ~countrycode)
d1_res <- feols(dy0 ~ gap_resid_cap | countrycode, data = dy0_df, vcov = ~countrycode)

d1_summary <- data.frame(
  what  = c("Raw OLS", "Within-country demeaned",
            "Country FE (raw gap)", "Country FE (gap_resid)"),
  slope = c(coef(d1_raw)[2], coef(d1_dm)[2],
            coef(d1_fe)["gap_capped"], coef(d1_res)["gap_resid_cap"]),
  se    = c(summary(d1_raw)$coefficients[2, 2],
            summary(d1_dm)$coefficients[2, 2],
            se(d1_fe)["gap_capped"], se(d1_res)["gap_resid_cap"])
) %>% mutate(t = slope / se)
print(d1_summary, row.names = FALSE)

p_d1_raw <- ggplot(dy0_df, aes(gap_capped, dy0)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_point(alpha = 0.35, color = "grey40") +
  geom_smooth(method = "lm", color = "steelblue", se = TRUE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Realized gap (capped at 20)",
       y = expression(y[t] - y[t - 1] ~ "(raw)"),
       title = "D1a: raw cross-section",
       subtitle = "One point per shock-year")

p_d1_dm <- ggplot(dy0_df, aes(gap_dm, dy0_dm)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_vline(xintercept = 0, color = "grey60") +
  geom_point(alpha = 0.35, color = "grey40") +
  geom_smooth(method = "lm", color = "steelblue", se = TRUE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Gap (demeaned within country)",
       y = expression(paste(y[t] - y[t - 1], ", demeaned")),
       title = "D1b: within-country",
       subtitle = "Does the positive slope survive country FE?")

p_d1 <- p_d1_raw + p_d1_dm
ggsave(file.path(out_fig, "d1_binscatter.png"), p_d1,
       width = 12, height = 5, dpi = 200)

# ==================================================================
# 2. STACKED PANEL FOR D2, D3, D4
# ==================================================================
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
        shock_l2    = lag(shock, 2),
        # Pre-event 3-year average growth, measured as of t-1
        pre_growth  = (lag(loggdp, 1) - lag(loggdp, 4)) / 3
      ) %>%
      ungroup() %>%
      filter(!is.na(dy), !is.na(gdp_diff_l1), !is.na(gdp_diff_l2)) %>%
      mutate(h = h)
  }
  sd <- bind_rows(out) %>%
    mutate(
      fh        = factor(h),
      h1        = h, h2 = h^2,
      hp1       = h + 1, hp1sq = (h + 1)^2,
      country_h = interaction(countrycode, fh, drop = TRUE)
    )
  # Horizon-specific shock and shock*gap_resid columns for D2
  for (hh in 0:H) {
    sd[[paste0("s_h",  hh)]] <- as.integer(sd$h == hh) * sd$shock
    sd[[paste0("sg_h", hh)]] <- as.integer(sd$h == hh) * sd$shock * sd$gap_resid_cap
  }
  # Polynomial columns for D3, D4
  sd <- sd %>%
    mutate(
      shock_h0       = shock,
      shock_h1       = shock * h1,
      shock_h2       = shock * h2,
      shock_gr_h0    = shock * gap_resid_cap,
      shock_gr_h1    = shock * gap_resid_cap * h1,
      shock_gr_h2    = shock * gap_resid_cap * h2,
      # Anchored parameterization: h' = h+1, no intercept
      shock_hp1      = shock * hp1,
      shock_hp1sq    = shock * hp1sq,
      shock_gr_hp1   = shock * gap_resid_cap * hp1,
      shock_gr_hp1sq = shock * gap_resid_cap * hp1sq,
      # Pre-growth interaction for D3
      shock_pg_h0    = shock * pre_growth,
      shock_pg_h1    = shock * pre_growth * h1,
      shock_pg_h2    = shock * pre_growth * h2
    )
  sd
}

sd_stack <- build_stacked(d_samp) %>% filter(!is.na(pre_growth))
cat(sprintf("Stacked rows (with pre_growth): %d\n", nrow(sd_stack)))

fe_base <- "country_h[year] + country_h[year2] + year^fh"

# ==================================================================
# D2. NON-PARAMETRIC HORIZON DUMMIES
# ==================================================================
cat("\n=== D2: non-parametric horizon dummies ===\n")

s_vars  <- paste0("s_h",  0:H)
sg_vars <- paste0("sg_h", 0:H)

fml_d2 <- as.formula(paste0(
  "dy ~ ", paste(c(s_vars, sg_vars), collapse = " + "),
  " + i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2)",
  " + i(fh, shock_l1) + i(fh, shock_l2) | ", fe_base
))
mod_d2 <- feols(fml_d2, data = sd_stack, vcov = ~countrycode)

beta_hat <- coef(mod_d2)[s_vars]
beta_se  <- se(mod_d2)[s_vars]
phi_hat  <- coef(mod_d2)[sg_vars]
phi_se   <- se(mod_d2)[sg_vars]

d2_df <- bind_rows(
  data.frame(horizon = 0:H, est = beta_hat, se = beta_se, what = "beta(h)"),
  data.frame(horizon = 0:H, est = phi_hat,  se = phi_se,  what = "phi(h)")
) %>%
  mutate(lo = est - 1.96*se, hi = est + 1.96*se,
         what = factor(what, levels = c("beta(h)", "phi(h)")))

cat("D2 beta(h) and phi(h):\n")
print(d2_df, row.names = FALSE)

p_d2 <- ggplot(d2_df, aes(horizon, est)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_pointrange(aes(ymin = lo, ymax = hi), color = "steelblue") +
  geom_line(color = "steelblue") +
  facet_wrap(~what, scales = "free_y", ncol = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = 0:H) +
  labs(x = "Horizon h (years)", y = NULL,
       title = "D2: non-parametric horizon dummies",
       subtitle = "If h=0 weirdness is a polynomial artefact, the dummies will flatten it")
ggsave(file.path(out_fig, "d2_nonparametric.png"), p_d2,
       width = 11, height = 4.5, dpi = 200)

# ==================================================================
# D3. CONTROL FOR PRE-EVENT STATE
# ==================================================================
cat("\n=== D3: control for pre-event 3-year growth ===\n")

# Reference: spec 2 polynomial with gap_resid (no pre_growth control)
fml_base_poly <- as.formula(paste0(
  "dy ~ shock_h0 + shock_h1 + shock_h2",
  " + shock_gr_h0 + shock_gr_h1 + shock_gr_h2",
  " + i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2)",
  " + i(fh, shock_l1) + i(fh, shock_l2) | ", fe_base
))
mod_d3_ref <- feols(fml_base_poly, data = sd_stack, vcov = ~countrycode)

# With pre_growth * shock control
fml_d3 <- as.formula(paste0(
  "dy ~ shock_h0 + shock_h1 + shock_h2",
  " + shock_gr_h0 + shock_gr_h1 + shock_gr_h2",
  " + shock_pg_h0 + shock_pg_h1 + shock_pg_h2",
  " + i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2)",
  " + i(fh, shock_l1) + i(fh, shock_l2) | ", fe_base
))
mod_d3 <- feols(fml_d3, data = sd_stack, vcov = ~countrycode)

# Predict IRFs at h=0 for g in {1, 5, 15} under both specs
predict_poly <- function(mod, g_values, gap_offset, has_pg = FALSE) {
  sv <- c("shock_h0", "shock_h1", "shock_h2",
          "shock_gr_h0", "shock_gr_h1", "shock_gr_h2")
  if (has_pg) sv <- c(sv, "shock_pg_h0", "shock_pg_h1", "shock_pg_h2")
  theta <- coef(mod)[sv]; V <- vcov(mod)[sv, sv, drop = FALSE]
  h_seq <- seq(0, H, by = 0.25)
  out <- list()
  for (g in g_values) {
    for (hh in h_seq) {
      xb <- hh^(0:2)
      g_eff <- g - gap_offset
      x <- c(xb, g_eff * xb)
      if (has_pg) x <- c(x, rep(0, 3))  # pre_growth at its mean (0)
      val <- sum(x * theta)
      s   <- as.numeric(sqrt(t(x) %*% V %*% x))
      out[[length(out)+1]] <- data.frame(horizon = hh, gap = g,
                                         irf = val, se = s,
                                         lo = val - 1.96*s,
                                         hi = val + 1.96*s)
    }
  }
  bind_rows(out)
}

# Center pre_growth at zero so the "no-pg" interpretation is unchanged
sd_stack$pre_growth <- sd_stack$pre_growth - mean(sd_stack$pre_growth, na.rm = TRUE)
# Recompute interactions
sd_stack <- sd_stack %>%
  mutate(
    shock_pg_h0 = shock * pre_growth,
    shock_pg_h1 = shock * pre_growth * h1,
    shock_pg_h2 = shock * pre_growth * h2
  )
mod_d3 <- feols(fml_d3, data = sd_stack, vcov = ~countrycode)

d3_ref <- predict_poly(mod_d3_ref, c(1, 5, 15), egap_mean, has_pg = FALSE) %>%
  mutate(spec = "(ref) no pre_growth control")
d3_ctrl <- predict_poly(mod_d3,   c(1, 5, 15), egap_mean, has_pg = TRUE) %>%
  mutate(spec = "(D3) + pre_growth x shock")

d3_df <- bind_rows(d3_ref, d3_ctrl) %>%
  mutate(gap = factor(gap))

p_d3 <- ggplot(d3_df, aes(horizon, irf, color = gap, fill = gap)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~spec, ncol = 2) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_color_manual(values = c("1" = "#d73027", "5" = "#fee08b", "15" = "#4575b4")) +
  scale_fill_manual(values  = c("1" = "#d73027", "5" = "#fee08b", "15" = "#4575b4")) +
  labs(x = "Horizon (years)", y = "Cumulative effect on log GDP pc",
       color = "Raw gap", fill = "Raw gap",
       title = "D3: residualized IRF, with and without pre-event growth control",
       subtitle = "If channel 3 (state dependence) is driving the h=0 anomaly, adding the control flattens it") +
  theme(legend.position = "bottom")
ggsave(file.path(out_fig, "d3_pregrowth.png"), p_d3,
       width = 12, height = 5, dpi = 200)

# Print h=0 IRFs at gap=15 under both specs
cat("D3 IRF at h=0, gap=15:\n")
cat(sprintf("  no pre_growth:   %.4f%% (ref)\n",
            100 * d3_ref$irf[d3_ref$horizon == 0 & d3_ref$gap == 15]))
cat(sprintf("  + pre_growth:    %.4f%%\n",
            100 * d3_ctrl$irf[d3_ctrl$horizon == 0 & d3_ctrl$gap == 15]))

# ==================================================================
# D4. ANCHORED PARAMETERIZATION (beta(-1) = 0)
# ==================================================================
cat("\n=== D4: anchored parameterization (h' = h+1) ===\n")

fml_d4 <- as.formula(paste0(
  "dy ~ shock_hp1 + shock_hp1sq",
  " + shock_gr_hp1 + shock_gr_hp1sq",
  " + i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2)",
  " + i(fh, shock_l1) + i(fh, shock_l2) | ", fe_base, " - 1"
))
# Note: the -1 is inside fe_base formula parsing; instead drop intercept via formula update
# feols typically doesn't have intercept because of FE absorption; but to be safe, use noconst
mod_d4 <- feols(
  dy ~ shock_hp1 + shock_hp1sq + shock_gr_hp1 + shock_gr_hp1sq +
    i(fh, gdp_diff_l1) + i(fh, gdp_diff_l2) +
    i(fh, shock_l1)   + i(fh, shock_l2) |
    country_h[year] + country_h[year2] + year^fh,
  data = sd_stack, vcov = ~countrycode
)

sv4 <- c("shock_hp1", "shock_hp1sq", "shock_gr_hp1", "shock_gr_hp1sq")
theta4 <- coef(mod_d4)[sv4]
V4     <- vcov(mod_d4)[sv4, sv4, drop = FALSE]

# Predict IRF at (h, g): h' = h+1
#   val = theta4[1]*(h+1) + theta4[2]*(h+1)^2 + (g - egap_mean)*(theta4[3]*(h+1) + theta4[4]*(h+1)^2)
predict_d4 <- function(g_values) {
  h_seq <- seq(-1, H, by = 0.25)
  out <- list()
  for (g in g_values) {
    for (hh in h_seq) {
      hp <- hh + 1
      g_eff <- g - egap_mean
      x <- c(hp, hp^2, g_eff * hp, g_eff * hp^2)
      val <- sum(x * theta4)
      s   <- as.numeric(sqrt(t(x) %*% V4 %*% x))
      out[[length(out)+1]] <- data.frame(horizon = hh, gap = g,
                                         irf = val, se = s,
                                         lo = val - 1.96*s,
                                         hi = val + 1.96*s)
    }
  }
  bind_rows(out)
}

d4_df <- predict_d4(c(1, 5, 15)) %>% mutate(gap = factor(gap))

p_d4 <- ggplot(d4_df, aes(horizon, irf, color = gap, fill = gap)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_vline(xintercept = -1, linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = c(-1, 0, 2, 4, 6, 8, 10)) +
  scale_color_manual(values = c("1" = "#d73027", "5" = "#fee08b", "15" = "#4575b4")) +
  scale_fill_manual(values  = c("1" = "#d73027", "5" = "#fee08b", "15" = "#4575b4")) +
  labs(x = expression("Horizon h (years); anchored at h = -1 via " * h + 1),
       y = "Cumulative effect on log GDP pc",
       color = "Raw gap", fill = "Raw gap",
       title = "D4: anchored parameterization (no intercept in (h+1))",
       subtitle = "Forces beta(-1) = 0 and phi(-1) = 0 mechanically") +
  theme(legend.position = "bottom")
ggsave(file.path(out_fig, "d4_anchored.png"), p_d4,
       width = 9, height = 5, dpi = 200)

cat("D4 IRF at h=0, gap=15 (anchored):\n")
cat(sprintf("  %.4f%%\n",
            100 * d4_df$irf[d4_df$horizon == 0 & d4_df$gap == 15]))

# ==================================================================
# DIAGNOSTICS FILE
# ==================================================================
diag_lines <- c(
  "Strategy 1: h=0 anomaly diagnostics",
  "====================================",
  sprintf("egap_mean = %.3f", egap_mean),
  "",
  "D1: raw same-year GDP change vs gap",
  paste0("  ", capture.output(print(d1_summary, row.names = FALSE))),
  "",
  "D2: non-parametric horizon dummies",
  sprintf("  beta(h=0) = %.4f (SE %.4f)  phi(h=0) = %.4f (SE %.4f)",
          beta_hat[1], beta_se[1], phi_hat[1], phi_se[1]),
  sprintf("  beta(h=5) = %.4f (SE %.4f)  phi(h=5) = %.4f (SE %.4f)",
          beta_hat[6], beta_se[6], phi_hat[6], phi_se[6]),
  "",
  "D3: effect of pre_growth control on IRF at h=0, gap=15",
  sprintf("  no control: %.4f%%",
          100 * d3_ref$irf[d3_ref$horizon == 0 & d3_ref$gap == 15]),
  sprintf("  + control:  %.4f%%",
          100 * d3_ctrl$irf[d3_ctrl$horizon == 0 & d3_ctrl$gap == 15]),
  "",
  "D4: anchored parameterization, IRF at h=0, gap=15",
  sprintf("  %.4f%%",
          100 * d4_df$irf[d4_df$horizon == 0 & d4_df$gap == 15])
)
writeLines(diag_lines, "strategy1_h0_diagnostics.txt")
cat("\nSaved: strategy1_h0_diagnostics.txt\nDone.\n")
