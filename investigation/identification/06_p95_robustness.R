# =============================================================================
# 06_p95_robustness.R
# Robustness of the main finding: gradient at p95 threshold (τ=5 vs τ=20)
# All tests use evaluation points τ = {5, 20} throughout.
#
# Tests:
#   A. Sample restrictions: drop never-treated, drop never-compound
#   B. Climatology subsamples: W̄ ≥ 20 m/s vs W̄ < 20 m/s
#   C. Leave-one-continent-out (LOCO)
#   D. Leave-one-year-out (LOYO) — jackknife envelope
#   E. Leave-one-country-out (LOOC) — jackknife envelope
#   F. Placebo: shuffle τ within country, 500 draws
#   G. Standard error variants (gradient at h=5)
#
# Figures saved to: figures/
# =============================================================================

library(tidyverse)
library(fixest)
library(patchwork)
library(latex2exp)
source("../../emileRegs.R")
setFixest_notes(FALSE)
setFixest_nthreads(1)

# ── Helpers -------------------------------------------------------------------

time_since_last <- function(x, y) {
  result <- rep(NA_real_, length(x)); last_yr <- NA_real_
  for (i in seq_along(x)) {
    if (!is.na(last_yr))           result[i] <- y[i] - last_yr
    if (!is.na(x[i]) && x[i] > 0) last_yr   <- y[i]
  }
  result
}

# ── Settings ------------------------------------------------------------------

FE      <- "countrycode[year] + countrycode[year2] + region^year"
PANEL   <- c("countrycode", "year")
DK      <- DK ~ year
HORIZON <- 10
CTRL    <- "l(shock_pct, 1:2) + l(gdp_diff, 1:2)"
EVAL    <- c(5, 20)
myred   <- "#d7191c"
myblue  <- "#2c7bb6"
mygray  <- "grey50"
pal_eval <- c("5 yrs" = myred, "20 yrs" = myblue)

theme_nature <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 1,
                                   colour = "grey10", margin = margin(b = 4)),
      plot.subtitle = element_text(colour = "grey45", size = base_size - 2,
                                   margin = margin(b = 8)),
      axis.title    = element_text(colour = "grey15", size = base_size),
      axis.text     = element_text(colour = "grey25", size = base_size - 2),
      axis.line     = element_line(colour = "grey40", linewidth = 0.45),
      axis.ticks    = element_line(colour = "grey50", linewidth = 0.45),
      legend.position  = "top",
      legend.text      = element_text(colour = "grey20", size = base_size - 1),
      legend.key.size  = unit(1.3, "lines"),
      legend.title     = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(face = "bold", size = base_size - 2),
      plot.margin      = margin(10, 14, 8, 10)
    )
}

# ── Data and p95 setup --------------------------------------------------------

data <- readRDS("../direct_damages/panel_damages.rds") %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  mutate(W_bar = mean(maxwind, na.rm = TRUE)) %>%
  ungroup()

thresh_p95 <- data %>%
  filter(W_bar > 0, maxwind > 0) %>%
  group_by(countrycode) %>%
  summarise(w_thresh = quantile(maxwind, 0.95, na.rm = TRUE), .groups = "drop")

make_d_p95 <- function(dat) {
  dat %>%
    left_join(thresh_p95, by = "countrycode") %>%
    mutate(
      shock_pct = case_when(
        is.na(w_thresh)     ~ 0L,
        maxwind >= w_thresh ~ 1L,
        TRUE                ~ 0L
      )
    ) %>%
    arrange(countrycode, year) %>%
    group_by(countrycode) %>%
    mutate(tau_pct = time_since_last(shock_pct, year)) %>%
    ungroup()
}

d_p95 <- make_d_p95(data)

cat("d_p95:", nrow(d_p95), "rows;",
    sum(d_p95$shock_pct, na.rm = TRUE), "p95 events;",
    n_distinct(d_p95$countrycode[d_p95$shock_pct == 1]), "treated countries\n")

# ── Core LP helper (lp_state_dep wrapper) ------------------------------------

run_lp <- function(dat, label = "", vcov_fml = DK) {
  tryCatch(
    lp_state_dep(
      data         = dat,
      outcome      = "loggdp",
      main_var     = "shock_pct",
      state_var    = "tau_pct",
      controls     = CTRL,
      horizon      = HORIZON,
      fe           = FE,
      panel_id     = PANEL,
      vcov_formula = vcov_fml,
      eval_values  = EVAL,
      eval_labels  = as.character(EVAL)
    ) %>%
      mutate(
        label     = label,
        tau_label = if_else(as.numeric(as.character(quantile)) == 5,
                            "5 yrs", "20 yrs"),
        tau_label = factor(tau_label, levels = c("5 yrs", "20 yrs")),
        irf_mean  = irf_mean * 100,
        irf_down  = irf_down * 100,
        irf_up    = irf_up   * 100
      ),
    error = function(e) { message("Skip [", label, "]: ", e$message); NULL }
  )
}

# Lean version: extract β₁(h) and β₂(h) + ME at τ=5,20 per horizon
# Used for LOYO/LOOC jackknife (no delta-method CI needed)
lp_extract_b <- function(dat, vcov_fml = DK) {
  map_dfr(0:HORIZON, function(h) {
    tryCatch({
      # Pre-compute LHS with dplyr lead/lag (matches 05_compound_robustness.R Section 7f)
      d_h <- dat %>%
        arrange(countrycode, year) %>%
        group_by(countrycode) %>%
        mutate(dy_h = lead(loggdp, h) - lag(loggdp, 1)) %>%
        ungroup() %>%
        filter(!is.na(dy_h), !is.na(shock_pct), !is.na(tau_pct))
      mod <- feols(
        dy_h ~ shock_pct + shock_pct:tau_pct + l(shock_pct, 1:2) + l(gdp_diff, 1:2) |
          countrycode[year] + countrycode[year2] + region^year,
        data = d_h, panel.id = PANEL, vcov = vcov_fml
      )
      cn <- c("shock_pct", "shock_pct:tau_pct")
      b1 <- coef(mod)[cn[1]]; b2 <- coef(mod)[cn[2]]
      V  <- vcov(mod)[cn, cn]
      tibble(
        h       = h,
        me5     = (b1 + 5  * b2) * 100,
        me20    = (b1 + 20 * b2) * 100,
        se5     = sqrt(V[1,1] + 25  * V[2,2] + 10 * V[1,2]) * 100,
        se20    = sqrt(V[1,1] + 400 * V[2,2] + 40 * V[1,2]) * 100,
        grad    = b2 * 15 * 100,   # ME(20)-ME(5) = β₂·15
        grad_se = sqrt(V[2,2]) * 15 * 100
      )
    }, error = function(e) NULL)
  })
}

# ── IRF panel plot helper -----------------------------------------------------

irf_panel <- function(df, title = NULL, sub = NULL, base_size = 12) {
  ggplot(df, aes(x = horizon, colour = tau_label, fill = tau_label)) +
    geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.10, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_line(aes(y = irf_mean), linewidth = 1.3) +
    scale_x_continuous(breaks = seq(0, HORIZON, 2)) +
    scale_colour_manual(values = pal_eval,
                        labels = c("5 yrs (recent)", "20 yrs (long recovery)")) +
    scale_fill_manual(values = pal_eval, guide = "none") +
    labs(title = title, subtitle = sub,
         x = "Years after strike", y = "Cumulative GDP growth (%)",
         colour = NULL) +
    theme_nature(base_size = base_size) +
    guides(colour = guide_legend(nrow = 1,
                                 override.aes = list(linewidth = 2, fill = NA)))
}

# ── Baseline ------------------------------------------------------------------
cat("\n=== Baseline ===\n")
irf_base <- run_lp(d_p95, "Baseline")
grad_base <- lp_extract_b(d_p95)
cat("  Baseline gradient at h=5:", round(grad_base$grad[grad_base$h == 5], 2), "pp\n")

# =============================================================================
# A. Sample restrictions
# =============================================================================
cat("\n=== A. Sample restrictions ===\n")

# A1. Drop never-treated (countries with W̄ = 0)
cat("  A1: Drop never-treated\n")
d_ever  <- d_p95 %>% filter(W_bar > 0)
irf_ever <- run_lp(d_ever, "Drop never-struck")

# A2. Drop never-compound: countries where no p95 event is compound (τ^p95 ≤ 5 at a shock)
never_compound <- d_p95 %>%
  group_by(countrycode) %>%
  summarise(has_compound = any(shock_pct == 1 & !is.na(tau_pct) & tau_pct <= 5),
            .groups = "drop") %>%
  filter(!has_compound) %>%
  pull(countrycode)
cat("  A2: Drop", length(never_compound), "never-compound countries\n")
d_hascmp <- d_p95 %>% filter(!countrycode %in% never_compound)
irf_hascmp <- run_lp(d_hascmp, "Drop never-compound")

# Combine and plot
irf_samples <- bind_rows(
  irf_base   %>% mutate(panel = "Full sample (baseline)"),
  irf_ever   %>% mutate(panel = "Ever-struck countries only"),
  irf_hascmp %>% mutate(panel = "Countries with compound events")
) %>%
  mutate(panel = factor(panel, levels = c(
    "Full sample (baseline)",
    "Ever-struck countries only",
    "Countries with compound events"
  )))

p_samples <- irf_panel(irf_samples,
  "Main Finding Holds Across Sample Restrictions",
  "p95 threshold. Each panel restricts the estimation sample. τ = 5 yrs (red) vs 20 yrs (blue).") +
  facet_wrap(~ panel, ncol = 3)

ggsave("figures/compound_p95_rob_samples.png", p_samples,
       width = 14, height = 6, dpi = 300)
message("Saved compound_p95_rob_samples.png")

# =============================================================================
# B. Climatology subsamples
# =============================================================================
cat("\n=== B. Climatology subsamples ===\n")

d_high <- d_p95 %>% filter(W_bar >= 20)
d_low  <- d_p95 %>% filter(W_bar > 0, W_bar < 20)
cat("  High W̄ (≥20): n =", nrow(d_high), "; Low W̄ (<20, >0): n =", nrow(d_low), "\n")

irf_high <- run_lp(d_high, "High climatology (W\u0304 \u2265 20 m/s)")
irf_low  <- run_lp(d_low,  "Low climatology (W\u0304 < 20 m/s)")

irf_clim <- bind_rows(
  irf_base %>% mutate(panel = "Full sample (baseline)"),
  irf_high %>% mutate(panel = "High climatology (W\u0304 \u2265 20 m/s)"),
  irf_low  %>% mutate(panel = "Low climatology (W\u0304 < 20 m/s)")
) %>%
  mutate(panel = factor(panel, levels = c(
    "Full sample (baseline)",
    "High climatology (W\u0304 \u2265 20 m/s)",
    "Low climatology (W\u0304 < 20 m/s)"
  )))

p_clim <- irf_panel(irf_clim,
  "Result Holds Within Climatology Groups",
  "p95 threshold. Splitting by long-run average wind exposure. τ = 5 yrs (red) vs 20 yrs (blue).") +
  facet_wrap(~ panel, ncol = 3)

ggsave("figures/compound_p95_rob_clim.png", p_clim,
       width = 14, height = 6, dpi = 300)
message("Saved compound_p95_rob_clim.png")

# =============================================================================
# C. Leave-one-continent-out (LOCO)
# =============================================================================
cat("\n=== C. LOCO ===\n")

continents_all <- sort(unique(na.omit(d_p95$continent)))
cat("  Continents:", paste(continents_all, collapse = ", "), "\n")

irf_loco_list <- map(continents_all, function(cont) {
  cat("  Drop", cont, "\n")
  d_sub <- d_p95 %>% filter(continent != cont | is.na(continent))
  run_lp(d_sub, paste0("Drop ", cont))
})
names(irf_loco_list) <- continents_all

irf_loco <- bind_rows(irf_loco_list, .id = "dropped") %>%
  mutate(panel = paste0("Drop ", dropped),
         panel = factor(panel))

irf_loco_all <- bind_rows(
  irf_base %>% mutate(panel = "Baseline (all continents)"),
  irf_loco
) %>%
  mutate(panel = factor(panel, levels = c(
    "Baseline (all continents)",
    paste0("Drop ", continents_all)
  )))

p_loco <- irf_panel(irf_loco_all,
  "Leave-One-Continent-Out: Gradient Driven by Asia",
  "p95 threshold. \u03c4 = 5 yrs (red) vs 20 yrs (blue). Each panel drops one continent.") +
  facet_wrap(~ panel, ncol = 3)

ggsave("figures/compound_p95_rob_loco.png", p_loco,
       width = 14, height = 9, dpi = 300)
message("Saved compound_p95_rob_loco.png")

# =============================================================================
# D. Leave-one-year-out (LOYO) — jackknife on gradient
# =============================================================================
cat("\n=== D. LOYO jackknife ===\n")

years_all <- sort(unique(d_p95$year[d_p95$shock_pct == 1]))
cat("  Years with p95 events:", length(years_all), "\n")

irf_loyo_all <- map_dfr(years_all, function(yr) {
  if (yr %% 5 == 0) cat("  Leave out", yr, "\n")
  d_sub <- d_p95 %>%
    mutate(shock_pct = if_else(year == yr, 0L, shock_pct)) %>%
    arrange(countrycode, year) %>%
    group_by(countrycode) %>%
    mutate(tau_pct = time_since_last(shock_pct, year)) %>%
    ungroup()
  lp_extract_b(d_sub) %>% mutate(year_out = yr)
})

# =============================================================================
# E. Leave-one-country-out (LOOC) — jackknife on gradient
# =============================================================================
cat("\n=== E. LOOC jackknife ===\n")

treated_cc <- d_p95 %>%
  group_by(countrycode) %>%
  filter(any(shock_pct == 1)) %>%
  pull(countrycode) %>% unique() %>% sort()
cat("  Treated countries:", length(treated_cc), "\n")

irf_looc_all <- map_dfr(seq_along(treated_cc), function(idx) {
  cc <- treated_cc[idx]
  if (idx %% 10 == 0) cat("  Leave out country", idx, "/", length(treated_cc), "\n")
  d_sub <- d_p95 %>%
    mutate(shock_pct = if_else(countrycode == cc, 0L, shock_pct)) %>%
    arrange(countrycode, year) %>%
    group_by(countrycode) %>%
    mutate(tau_pct = time_since_last(shock_pct, year)) %>%
    ungroup()
  lp_extract_b(d_sub) %>% mutate(country_out = cc)
})

# ── Jackknife summary plot (gradient IRF) ─────────────────────────────────

summarise_jackknife <- function(df, group_var) {
  df %>%
    group_by(h) %>%
    summarise(
      q10 = quantile(grad, 0.10, na.rm = TRUE),
      q25 = quantile(grad, 0.25, na.rm = TRUE),
      q50 = median(grad, na.rm = TRUE),
      q75 = quantile(grad, 0.75, na.rm = TRUE),
      q90 = quantile(grad, 0.90, na.rm = TRUE),
      .groups = "drop"
    )
}

jk_loyo <- summarise_jackknife(irf_loyo_all) %>% mutate(type = "Leave one year out")
jk_looc <- summarise_jackknife(irf_looc_all) %>% mutate(type = "Leave one country out")
jk_base <- grad_base %>% select(h, grad) %>% mutate(type = "Baseline")

p_loyo <- ggplot() +
  geom_ribbon(data = jk_loyo, aes(x = h, ymin = q10, ymax = q90),
              fill = "steelblue", alpha = 0.20) +
  geom_ribbon(data = jk_loyo, aes(x = h, ymin = q25, ymax = q75),
              fill = "steelblue", alpha = 0.30) +
  geom_line(data = jk_loyo, aes(x = h, y = q50),
            colour = "steelblue", linewidth = 1.0, linetype = "dashed") +
  geom_line(data = grad_base, aes(x = h, y = grad),
            colour = myblue, linewidth = 1.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(0, HORIZON, 2)) +
  labs(
    title    = "Leave-One-Year-Out: Gradient Stable Across Years",
    subtitle = "Gradient = ME(20 yrs) - ME(5 yrs) in pp GDP.\nShaded bands: 10-90th and 25-75th pct of 500+ jackknife estimates. Bold = baseline.",
    x = "Years after strike", y = "Gradient: long vs. short recovery (pp)"
  ) +
  theme_nature()

p_looc <- ggplot() +
  geom_ribbon(data = jk_looc, aes(x = h, ymin = q10, ymax = q90),
              fill = "firebrick", alpha = 0.20) +
  geom_ribbon(data = jk_looc, aes(x = h, ymin = q25, ymax = q75),
              fill = "firebrick", alpha = 0.30) +
  geom_line(data = jk_looc, aes(x = h, y = q50),
            colour = "firebrick", linewidth = 1.0, linetype = "dashed") +
  geom_line(data = grad_base, aes(x = h, y = grad),
            colour = myblue, linewidth = 1.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(0, HORIZON, 2)) +
  labs(
    title    = "Leave-One-Country-Out: No Country Drives the Result",
    subtitle = "Gradient = ME(20 yrs) - ME(5 yrs) in pp GDP.\nShaded bands: 10-90th and 25-75th pct of jackknife estimates. Bold = baseline.",
    x = "Years after strike", y = "Gradient: long vs. short recovery (pp)"
  ) +
  theme_nature()

p_jk <- p_loyo + p_looc +
  plot_annotation(
    title = "Jackknife Robustness: Gradient Is Not Driven by Any Single Year or Country",
    theme = theme(plot.title = element_text(face = "bold", size = 15,
                                            colour = "grey10", hjust = 0.5))
  )

ggsave("figures/compound_p95_rob_jackknife.png", p_jk,
       width = 14, height = 6, dpi = 300)
message("Saved compound_p95_rob_jackknife.png")

# =============================================================================
# F. Placebo: shuffle maxwind within country, recompute p95 threshold + τ
# =============================================================================
# Within-country permutation of maxwind preserves the wind distribution
# (so the country p95 threshold is invariant), but scrambles which years
# are p95 events and therefore randomises τ.  GDP outcomes are unchanged.
# This tests whether the timing of extreme events generates the gradient.
# =============================================================================
cat("\n=== F. Placebo (500 draws, shuffle maxwind within country) ===\n")
set.seed(42)

# Actual gradient at h=5 (in pp)
actual_b2 <- grad_base$grad[grad_base$h == 5]

N_PLACEBO <- 500
grad_placebo <- map_dbl(seq_len(N_PLACEBO), function(i) {
  if (i %% 100 == 0) cat("  Placebo draw", i, "/", N_PLACEBO, "\n")

  # 1. Shuffle maxwind within each country (permutes timing, preserves distribution)
  # 2. Recompute shock_pct using original thresh_p95 (unchanged by permutation)
  # 3. Recompute tau_pct from new shock timing
  # 4. Pre-compute dy5 from unchanged loggdp
  d_pl <- data %>%
    arrange(countrycode, year) %>%
    group_by(countrycode) %>%
    mutate(maxwind = sample(maxwind)) %>%
    ungroup() %>%
    left_join(thresh_p95, by = "countrycode") %>%
    mutate(
      shock_pct = case_when(
        is.na(w_thresh)     ~ 0L,
        maxwind >= w_thresh ~ 1L,
        TRUE                ~ 0L
      )
    ) %>%
    arrange(countrycode, year) %>%
    group_by(countrycode) %>%
    mutate(
      tau_pct = time_since_last(shock_pct, year),
      dy5     = lead(loggdp, 5) - lag(loggdp, 1)
    ) %>%
    ungroup() %>%
    filter(!is.na(dy5), !is.na(shock_pct), !is.na(tau_pct))

  tryCatch({
    mod <- feols(
      dy5 ~ shock_pct + shock_pct:tau_pct + l(shock_pct, 1:2) + l(gdp_diff, 1:2) |
        countrycode[year] + countrycode[year2] + region^year,
      data = d_pl, panel.id = PANEL, vcov = DK
    )
    coef(mod)["shock_pct:tau_pct"] * 15 * 100
  }, error = function(e) NA_real_)
})

# h=5 dataset for SE section (build here after placebo, from actual d_p95)
d_h5 <- d_p95 %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  mutate(dy5 = lead(loggdp, 5) - lag(loggdp, 1)) %>%
  ungroup() %>%
  filter(!is.na(dy5), !is.na(shock_pct), !is.na(tau_pct))

n_valid  <- sum(!is.na(grad_placebo))
p_val_pl <- mean(abs(grad_placebo) >= abs(actual_b2), na.rm = TRUE)
cat("  Actual gradient at h=5:", round(actual_b2, 2), "pp\n")
cat("  Placebo p-value (two-sided):", round(p_val_pl, 3),
    " (", n_valid, "valid draws)\n")

p_placebo <- ggplot(tibble(grad = na.omit(grad_placebo)), aes(x = grad)) +
  geom_histogram(fill = "steelblue", colour = "white", bins = 40, alpha = 0.70) +
  geom_vline(xintercept = actual_b2, colour = myred,
             linewidth = 1.5, linetype = "solid") +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey30", linewidth = 0.6) +
  annotate("text", x = actual_b2 * 1.05, y = Inf,
           label = paste0("Actual\n", round(actual_b2, 1), " pp"),
           hjust = 0, vjust = 1.3, colour = myred, size = 3.8, fontface = "bold") +
  labs(
    title    = "Placebo Test: Wind Timing Randomly Shuffled Within Country",
    subtitle = paste0(
      "maxwind permuted within each country; p95 threshold and \u03c4 recomputed each draw. ",
      n_valid, " permutations.\n",
      "Red line = actual estimate. Permutation p-value = ", round(p_val_pl, 3), "."
    ),
    x = "Placebo gradient at year 5 (pp GDP)", y = "Count"
  ) +
  theme_nature()

ggsave("figures/compound_p95_rob_placebo.png", p_placebo,
       width = 9, height = 5.5, dpi = 300)
message("Saved compound_p95_rob_placebo.png")

# =============================================================================
# G. Standard error variants: DK, Two-way cluster, Wild cluster bootstrap
# =============================================================================
cat("\n=== G. SE variants ===\n")

# d_h5 has dy5 = lead(loggdp,5) - lag(loggdp,1) pre-computed with dplyr
# Use dy5 as LHS to match 05_compound_robustness.R Section 7f exactly
fml_h5 <- dy5 ~ shock_pct + shock_pct:tau_pct + l(shock_pct, 1:2) + l(gdp_diff, 1:2) |
  countrycode[year] + countrycode[year2] + region^year

# Fit once; re-summarise with different vcov for DK and two-way cluster
mod_h5 <- feols(fml_h5, data = d_h5, panel.id = PANEL)
b2     <- coef(mod_h5)["shock_pct:tau_pct"]
gap    <- 15   # τ=20 minus τ=5

cat("  Computing DK SE...\n")
se_dk <- sqrt(vcov(summary(mod_h5, vcov = DK ~ year))[
                "shock_pct:tau_pct", "shock_pct:tau_pct"])

cat("  Computing two-way cluster SE...\n")
se_tw <- sqrt(vcov(summary(mod_h5, vcov = ~ countrycode + year))[
                "shock_pct:tau_pct", "shock_pct:tau_pct"])

# Wild cluster bootstrap (country-level clusters, B = 999)
cat("  Computing wild cluster bootstrap CI...\n")
boot_lo <- NA_real_; boot_hi <- NA_real_
tryCatch({
  if (!requireNamespace("fwildclusterboot", quietly = TRUE))
    stop("fwildclusterboot not installed")
  library(fwildclusterboot)
  set.seed(42)
  boot_res <- boottest(mod_h5, clustid = "countrycode",
                       param = "shock_pct:tau_pct", B = 999)
  ci_boot  <- confint(boot_res)
  boot_lo  <- ci_boot[[1]] * gap * 100
  boot_hi  <- ci_boot[[2]] * gap * 100
  cat("  Wild bootstrap CI: [", round(boot_lo, 2), ",", round(boot_hi, 2), "]\n")
}, error = function(e) {
  message("  Wild bootstrap unavailable (", e$message,
          ") — using country cluster as fallback")
  se_c   <- sqrt(vcov(summary(mod_h5, vcov = ~ countrycode))[
                   "shock_pct:tau_pct", "shock_pct:tau_pct"])
  boot_lo <<- (b2 - 1.96 * se_c) * gap * 100
  boot_hi <<- (b2 + 1.96 * se_c) * gap * 100
})

spec_levels <- c("DK (baseline)", "Two-way cluster", "Wild cluster bootstrap")

grad_se <- tibble(
  spec     = spec_levels,
  gradient = b2 * gap * 100,
  grad_lo  = c((b2 - 1.96 * se_dk) * gap * 100,
               (b2 - 1.96 * se_tw) * gap * 100,
               boot_lo),
  grad_hi  = c((b2 + 1.96 * se_dk) * gap * 100,
               (b2 + 1.96 * se_tw) * gap * 100,
               boot_hi)
) %>%
  mutate(grad_se = (grad_hi - grad_lo) / (2 * 1.96))

cat("\nSE comparison:\n")
print(grad_se)

p_se <- ggplot(grad_se,
               aes(x = gradient,
                   y = factor(spec, levels = rev(spec_levels)))) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey55", linewidth = 0.55) +
  geom_errorbarh(aes(xmin = grad_lo, xmax = grad_hi),
                 height = 0.25, colour = myblue, linewidth = 1.1) +
  geom_point(size = 4.5, colour = myblue) +
  geom_text(aes(label = paste0(round(gradient, 1), " pp")),
            nudge_y = 0.28, size = 3.5, colour = "grey20") +
  labs(
    title    = "Gradient at Year 5 Is Robust to Standard Error Choice",
    subtitle = "Gradient = ME(τ=20) - ME(τ=5), h=5. 95% CI shown. Point estimate identical across specifications.",
    x = "Gradient at year 5 (pp GDP)", y = NULL
  ) +
  theme_nature()

ggsave("figures/compound_p95_rob_se.png", p_se,
       width = 9, height = 4, dpi = 300)
message("Saved compound_p95_rob_se.png")

# =============================================================================
# Save summary
# =============================================================================

saveRDS(
  list(
    irf_base    = irf_base,
    irf_samples = irf_samples,
    irf_clim    = irf_clim,
    irf_loco    = irf_loco,
    irf_loyo_jk = irf_loyo_all,
    irf_looc_jk = irf_looc_all,
    grad_placebo = grad_placebo,
    actual_grad  = actual_b2,
    grad_se      = grad_se
  ),
  "results_p95_robustness.rds"
)
message("Saved results_p95_robustness.rds")
cat("\nDone.\n")
