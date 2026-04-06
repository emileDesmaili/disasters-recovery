# =============================================================================
# Compound LP Robustness: Is the τ Gradient Real or W̄ Confounding?
# =============================================================================
# Central question: the compound LP finds that isolated events (τ=5) generate
# larger GDP losses than compounding events (τ=1).  But τ is correlated with
# W̄: high-climatology countries have shorter gaps by construction.  This
# script probes how much of the gradient survives controls for W̄.
#
# 1. Descriptive: where do compound events occur? (continent, W_bar)
# 2. τ distribution diagnostic: is there variation in τ within clim2 groups?
# 3. Full sample vs ever-treated
# 4. Leave-one-continent-out
# 5. Within climatology groups — binary (W̄ ≥ 20 vs < 20)
# 6. Within climatology groups — terciles (finer W̄ bins)
# 7. Country-specific percentile shock × compound (p75/p90/p95/p99)
# =============================================================================

library(tidyverse)
library(fixest)
library(patchwork)
library(latex2exp)
source("../../emileRegs.R")
setFixest_notes(FALSE)
setFixest_nthreads(1)

# ── Helper -------------------------------------------------------------------

time_since_last <- function(x, y) {
  result <- rep(NA_real_, length(x)); last_yr <- NA_real_
  for (i in seq_along(x)) {
    if (!is.na(last_yr))           result[i] <- y[i] - last_yr
    if (!is.na(x[i]) && x[i] > 0) last_yr   <- y[i]
  }
  result
}

# ── Data ---------------------------------------------------------------------

data <- readRDS("../direct_damages/panel_damages.rds") %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  mutate(
    W_bar           = mean(maxwind, na.rm = TRUE),
    years_since_any = time_since_last(maxwind, year)
  ) %>%
  ungroup()

# Binary climatology split: High = W̄ ≥ 20 m/s, Low = 0 < W̄ < 20 m/s
data$clim2 <- ifelse(data$W_bar <= 0, NA_character_,
                     ifelse(data$W_bar >= 20, "High", "Low"))

# Climatology three-way split with fixed breaks at 10 and 30 m/s
tert_breaks <- c(10, 30)
data$clim3 <- ifelse(
  data$W_bar <= 0, NA_character_,
  ifelse(data$W_bar <= tert_breaks[1], "T1 (0\u201310 m/s)",
         ifelse(data$W_bar <= tert_breaks[2], "T2 (10\u201330 m/s)", "T3 (>30 m/s)"))
)

cat("W_bar breaks (m/s):", tert_breaks, "\n")
cat("clim2 distribution:\n")
print(table(data$clim2, useNA = "ifany"))
cat("clim3 distribution:\n")
print(table(data$clim3, useNA = "ifany"))

FE      <- "countrycode[year] + countrycode[year2] + region^year"
PANEL   <- c("countrycode", "year")
DK      <- DK ~ year
HORIZON <- 10

# ── Residualise maxwind on country linear trend + region×year FE --------------
resid_mod <- feols(maxwind ~ 1 | countrycode[year] + region^year,
                   data = data, panel.id = PANEL)
data$maxwind_resid <- residuals(resid_mod)
cat("Correlation(maxwind, maxwind_resid):",
    round(cor(data$maxwind, data$maxwind_resid, use = "complete.obs"), 3), "\n")

CONTROLS <- "l(maxwind_resid, 1:2) + l(gdp_diff, 1:2)"

eval_vals   <- c(1, 3, 5)
eval_labels <- c("1 yr", "3 yrs", "5 yrs")

myblue  <- "#2c7bb6"; myred <- "#d7191c"; mygold <- "#fdae61"; mygreen <- "#1a9641"
pal_since <- c("1 yr" = myred, "3 yrs" = mygold, "5 yrs" = myblue)

theme_nature <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 1,
                                   colour = "grey10", margin = margin(b = 4)),
      plot.subtitle = element_text(colour = "grey45", size = base_size - 2,
                                   margin = margin(b = 8)),
      axis.title    = element_text(colour = "grey15"),
      axis.text     = element_text(colour = "grey25", size = base_size - 2),
      axis.line     = element_line(colour = "grey40", linewidth = 0.45),
      axis.ticks    = element_line(colour = "grey50", linewidth = 0.45),
      legend.position  = "top",
      legend.text      = element_text(colour = "grey20", size = base_size - 1),
      legend.key.size  = unit(1.2, "lines"),
      legend.title     = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      plot.margin      = margin(8, 12, 6, 8)
    )
}

irf_ribbon <- function(df, pal = pal_since, facet_var = NULL) {
  p <- ggplot(df, aes(x = horizon, colour = quantile, fill = quantile)) +
    geom_ribbon(aes(ymin = irf_down * 100, ymax = irf_up * 100),
                alpha = 0.10, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_line(aes(y = irf_mean * 100), linewidth = 1.5) +
    scale_x_continuous(breaks = seq(0, HORIZON, 2)) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values   = pal) +
    labs(x = "Years after strike", y = "Cumulative GDP growth (%)") +
    theme_nature() +
    guides(colour = guide_legend(nrow = 1))
  if (!is.null(facet_var))
    p <- p + facet_wrap(as.formula(paste("~", facet_var)), ncol = 2)
  p
}

run_lp <- function(d, label, state_var = "years_since_any",
                   eval_v = eval_vals, eval_l = eval_labels) {
  cat("  Running:", label, "| n =", nrow(d), "\n")
  tryCatch(
    lp_state_dep(
      data         = d,
      outcome      = "loggdp",
      main_var     = "maxwind_resid",
      state_var    = state_var,
      controls     = CONTROLS,
      horizon      = HORIZON,
      fe           = FE,
      panel_id     = PANEL,
      vcov_formula = DK,
      eval_values  = eval_v,
      eval_labels  = eval_l
    ) %>% mutate(sample = label),
    error = function(e) { message("  Skip ", label, ": ", e$message); NULL }
  )
}

# =============================================================================
# 1. Descriptive: where do compound events occur?
# =============================================================================

compound_desc <- data %>%
  filter(maxwind > 0, !is.na(years_since_any)) %>%
  mutate(
    compound_type = case_when(
      years_since_any <= 2  ~ "Compound (\u22642 yrs)",
      years_since_any <= 5  ~ "Medium (3\u20135 yrs)",
      TRUE                  ~ "Isolated (>5 yrs)"
    ),
    compound_type = factor(compound_type,
      levels = c("Compound (\u22642 yrs)", "Medium (3\u20135 yrs)", "Isolated (>5 yrs)"))
  )

cat("Compound event classification:\n")
print(table(compound_desc$compound_type))

p_desc_cont <- compound_desc %>%
  filter(!is.na(continent)) %>%
  count(continent, compound_type) %>%
  group_by(continent) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = reorder(continent, -pct), y = pct, fill = compound_type)) +
  geom_col(width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = c(myred, mygold, myblue),
                    labels = c(TeX("Compound ($\\leq$2 yrs)"), "Medium (3-5 yrs)", "Isolated (>5 yrs)")) +
  labs(title = "Composition of Strike Events by Continent",
       subtitle = "Share compound / medium / isolated",
       x = NULL, y = "% of strikes", fill = NULL) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_desc_wbar <- compound_desc %>%
  ggplot(aes(x = W_bar, fill = compound_type, colour = compound_type)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  scale_fill_manual(values = c(myred, mygold, myblue),
                    labels = c(TeX("Compound ($\\leq$2 yrs)"), "Medium (3-5 yrs)", "Isolated (>5 yrs)")) +
  scale_colour_manual(values = c(myred, mygold, myblue),
                      labels = c(TeX("Compound ($\\leq$2 yrs)"), "Medium (3-5 yrs)", "Isolated (>5 yrs)")) +
  labs(title = "Country Climatology Distribution by Event Type",
       subtitle = TeX("Compound events cluster in high-$\\bar{W}$ countries"),
       x = expression(bar(W)[i]~"(m/s)"), y = "Density",
       fill = NULL, colour = NULL) +
  theme_nature()

p_desc <- p_desc_cont / p_desc_wbar
ggsave("figures/compound_rob_descriptive.png", p_desc,
       width = 10, height = 10, dpi = 300)
message("Saved compound_rob_descriptive.png")

# =============================================================================
# 2. τ distribution diagnostic within clim2 groups
# =============================================================================
# Key check: is there actually variation in τ within each climatology bin?
# If high-clim countries have τ ∈ {1,2} only and low-clim have τ ∈ {5,10+},
# then the "within-group" LP is identified from very different τ ranges —
# not a fair within-group comparison.

cat("\n--- τ distribution within clim2 groups ---\n")
tau_by_clim <- data %>%
  filter(!is.na(clim2), !is.na(years_since_any)) %>%
  group_by(clim2) %>%
  summarise(
    n           = n(),
    tau_mean    = mean(years_since_any),
    tau_sd      = sd(years_since_any),
    tau_p10     = quantile(years_since_any, 0.10),
    tau_p25     = quantile(years_since_any, 0.25),
    tau_p50     = median(years_since_any),
    tau_p75     = quantile(years_since_any, 0.75),
    tau_p90     = quantile(years_since_any, 0.90),
    pct_tau_le2 = mean(years_since_any <= 2) * 100,
    pct_tau_ge5 = mean(years_since_any >= 5) * 100,
    .groups = "drop"
  )
cat("\nτ distribution by clim2 group:\n")
print(tau_by_clim %>% mutate(across(where(is.numeric), ~round(., 1))) %>% as.data.frame())

# Also within clim3 tercile groups
cat("\nτ mean and SD by tercile group:\n")
data %>%
  filter(!is.na(clim3), !is.na(years_since_any)) %>%
  group_by(clim3) %>%
  summarise(
    n        = n(),
    W_bar_mn = mean(W_bar),
    tau_mean = mean(years_since_any),
    tau_sd   = sd(years_since_any),
    pct_le2  = mean(years_since_any <= 2) * 100,
    pct_ge5  = mean(years_since_any >= 5) * 100,
    .groups  = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 1))) %>%
  as.data.frame() %>% print()

# Boxplot: τ distribution within clim2 groups
p_tau_dist <- data %>%
  filter(!is.na(clim2), !is.na(years_since_any)) %>%
  mutate(clim2 = factor(clim2, levels = c("Low", "High"))) %>%
  ggplot(aes(x = clim2, y = years_since_any, fill = clim2)) +
  geom_violin(alpha = 0.4, colour = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.5, fill = "white", colour = "grey30") +
  scale_fill_manual(values = c("Low" = myblue, "High" = myred), guide = "none") +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50)) +
  labs(
    title    = TeX("$\\tau$ Distribution Within Climatology Groups"),
    subtitle = TeX("Log scale. Both groups have substantial $\\tau$ variation; $\\bar{W}$ confound not fully removed by binary split"),
    x = TeX("Climatology group ($\\bar{W}$ $\\geq$ 20 m/s = High)"),
    y = TeX("Years since last TC ($\\tau$, log scale)")
  ) +
  theme_nature()

# Scatter: τ vs W̄ with group colouring and within-group LOESS
p_tau_wbar <- data %>%
  filter(W_bar > 0, !is.na(years_since_any), !is.na(clim2)) %>%
  mutate(clim2 = factor(clim2, levels = c("Low", "High"))) %>%
  ggplot(aes(x = W_bar, y = years_since_any, colour = clim2)) +
  geom_point(alpha = 0.05, size = 0.6) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.2) +
  geom_vline(xintercept = 20, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("Low" = myblue, "High" = myred)) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50)) +
  labs(
    title    = TeX("$\\tau$ vs $\\bar{W}$: Confound Persists Within Binary Groups"),
    subtitle = TeX("Within Low (blue) and High (red), higher $\\bar{W}$ still predicts shorter $\\tau$"),
    x = expression(bar(W)[i]~"(m/s)"),
    y = TeX("Years since last TC ($\\tau$, log scale)"),
    colour = NULL
  ) +
  theme_nature()

p_tau_diag <- p_tau_dist | p_tau_wbar
ggsave("figures/compound_rob_tau_diagnostic.png", p_tau_diag,
       width = 12, height = 6, dpi = 300)
message("Saved compound_rob_tau_diagnostic.png")

# =============================================================================
# 3. Full sample vs ever-treated
# =============================================================================

cat("\n--- Full sample vs ever-treated ---\n")
irf_samples <- bind_rows(
  run_lp(data,                       "Full sample"),
  run_lp(filter(data, W_bar > 0),    "Ever-treated")
) %>%
  mutate(quantile = factor(quantile, levels = eval_labels),
         sample   = factor(sample,   levels = c("Full sample", "Ever-treated")))

p_samples <- irf_ribbon(irf_samples, facet_var = "sample") +
  labs(title    = "Full Sample vs Ever-Treated",
       subtitle = TeX("Compound gradient with $\\tau$ $\\in$ \\{1,3,5\\} yrs"))

ggsave("figures/compound_rob_samples.png", p_samples,
       width = 11, height = 6, dpi = 300)
message("Saved compound_rob_samples.png")

# =============================================================================
# 3b. W̄ minimum-threshold subsamples
# =============================================================================
# Idea: progressively restrict to countries with W̄ ≥ {0, 10, 15, 20} m/s.
# As the threshold rises we keep only higher-climatology countries — those
# with shorter τ by construction.  If the compound gradient is a W̄–τ
# confound, it should *attenuate* as the threshold rises (less τ variation
# remains; low-τ and high-τ observations become more similar in W̄).
# If instead the gradient reflects a genuine time-since mechanism, it should
# persist or even sharpen once we focus on the countries that actually get hit.
#
# For each threshold we report:
#   • τ distribution stats (mean, % isolated ≥ 5 yrs)
#   • Full IRF at τ ∈ {1, 3, 5} yrs
#   • Gradient = ME(τ=5) − ME(τ=1) at h=5

wbar_thresholds <- c(0, 10, 15, 20)

cat("\n--- W̄ minimum-threshold subsamples ---\n")
cat("τ distribution by subsample:\n")
data %>%
  filter(!is.na(years_since_any)) %>%
  mutate(dummy = TRUE) %>%
  bind_rows(
    map_dfr(wbar_thresholds, function(th) {
      data %>%
        filter(W_bar >= th, !is.na(years_since_any)) %>%
        mutate(threshold = paste0("W\u0304 \u2265 ", th, " m/s"))
    })
  ) %>%
  {
    map_dfr(wbar_thresholds, function(th) {
      d <- filter(data, W_bar >= th, !is.na(years_since_any))
      tibble(
        threshold    = paste0("W\u0304 \u2265 ", th, " m/s"),
        n_ctry       = n_distinct(d$countrycode),
        n_obs        = nrow(d),
        W_bar_mean   = round(mean(d$W_bar, na.rm = TRUE), 1),
        tau_mean     = round(mean(d$years_since_any), 1),
        tau_p50      = round(median(d$years_since_any), 1),
        pct_le2      = round(mean(d$years_since_any <= 2) * 100, 1),
        pct_ge5      = round(mean(d$years_since_any >= 5) * 100, 1)
      )
    })
  } %>%
  as.data.frame() %>% print()

thresh_labels <- paste0("W\u0304 \u2265 ", wbar_thresholds, " m/s")

irf_wbar_thresh <- map_dfr(seq_along(wbar_thresholds), function(i) {
  th  <- wbar_thresholds[i]
  lbl <- thresh_labels[i]
  run_lp(filter(data, W_bar >= th), lbl)
}) %>%
  mutate(
    quantile  = factor(quantile, levels = eval_labels),
    sample    = factor(sample,   levels = thresh_labels)
  )

# ── IRF panels (2 × 2) -------------------------------------------------------
p_wbar_thresh <- irf_ribbon(irf_wbar_thresh, facet_var = "sample") +
  labs(
    title    = TeX("Compound LP by $\\bar{W}$ Minimum Threshold"),
    subtitle = TeX("Each panel restricts to countries with $\\bar{W}$ $\\geq$ threshold.  As threshold rises, the $\\tau$ gradient should attenuate if it is a $\\bar{W}$ confound.")
  )

ggsave("figures/compound_rob_wbar_thresh.png", p_wbar_thresh,
       width = 12, height = 10, dpi = 300)
message("Saved compound_rob_wbar_thresh.png")

# ── Gradient at h = {1..10} for each threshold (line chart) ------------------
# Shows how the gradient evolves with the horizon as well as with the threshold.
grad_by_thresh <- irf_wbar_thresh %>%
  filter(quantile %in% c("5 yrs", "1 yr")) %>%
  select(sample, horizon, quantile, irf_mean, se) %>%
  pivot_wider(names_from = quantile,
              values_from = c(irf_mean, se),
              names_repair = "unique") %>%
  mutate(
    gradient = (`irf_mean_5 yrs` - `irf_mean_1 yr`) * 100,
    grad_se  = sqrt(`se_5 yrs`^2 + `se_1 yr`^2) * 100,
    grad_lo  = gradient - 1.96 * grad_se,
    grad_hi  = gradient + 1.96 * grad_se
  )

thresh_pal <- c(
  "W\u0304 \u2265 0 m/s"  = "grey50",
  "W\u0304 \u2265 10 m/s" = mygreen,
  "W\u0304 \u2265 15 m/s" = mygold,
  "W\u0304 \u2265 20 m/s" = myred
)

p_grad_thresh <- ggplot(grad_by_thresh,
                        aes(x = horizon, y = gradient,
                            colour = sample, fill = sample)) +
  geom_ribbon(aes(ymin = grad_lo, ymax = grad_hi),
              alpha = 0.12, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 1.5) +
  scale_x_continuous(breaks = seq(0, HORIZON, 2)) +
  scale_colour_manual(values = thresh_pal,
                      labels = c(TeX("$\\bar{W}$ $\\geq$ 0 m/s"),
                                 TeX("$\\bar{W}$ $\\geq$ 10 m/s"),
                                 TeX("$\\bar{W}$ $\\geq$ 15 m/s"),
                                 TeX("$\\bar{W}$ $\\geq$ 20 m/s"))) +
  scale_fill_manual(values   = thresh_pal,
                    labels = c(TeX("$\\bar{W}$ $\\geq$ 0 m/s"),
                               TeX("$\\bar{W}$ $\\geq$ 10 m/s"),
                               TeX("$\\bar{W}$ $\\geq$ 15 m/s"),
                               TeX("$\\bar{W}$ $\\geq$ 20 m/s"))) +
  labs(
    title    = TeX("Isolated-minus-Compounding Gradient by $\\bar{W}$ Threshold"),
    subtitle = TeX("ME($\\tau$=5) - ME($\\tau$=1) at each horizon h; positive = isolated more damaging"),
    x = "Years after strike",
    y = "Gradient (pp)",
    colour = NULL, fill = NULL
  ) +
  theme_nature() +
  guides(colour = guide_legend(nrow = 2))

ggsave("figures/compound_rob_wbar_grad.png", p_grad_thresh,
       width = 9, height = 6, dpi = 300)
message("Saved compound_rob_wbar_grad.png")

# =============================================================================
# 3c. Smooth parametric LP with quadratic interaction — W̄ threshold subsamples
# =============================================================================
# Method: stacked LP (all horizons h = 0..10 stacked into one dataset).
# The shock's IRF is parameterised as a degree-2 polynomial in h:
#
#   ME(h, τ=g) = (θ₀ + θ₁h + θ₂h²) + g·(θ₃ + θ₄h + θ₅h²)
#             := β(h) + g·γ(h)
#
# Stacked model:
#   dy_{i,t,h} = Σ_{k=0}^2 θ_k h^k w̃_{it}
#              + Σ_{k=0}^2 θ_{3+k} h^k (w̃_{it} × τ_{it})
#              + i(fh, lag controls)
#              | country_h[year] + country_h[year2] + year^fh
#
# where country_h = countrycode × fh (per-horizon country quadratic trends)
# and   year^fh   = year × fh   (per-horizon year FEs)
#
# SE: clustered by countrycode (DK not available in stacked estimator).
# Applied to the four W̄ threshold subsamples: ≥ {0, 10, 15, 20} m/s.
# =============================================================================

build_stacked_id <- function(dat, H = 10) {
  map_dfr(0:H, function(h) {
    dat %>%
      group_by(countrycode) %>%
      arrange(year) %>%
      mutate(
        dy        = lead(loggdp, h) - lag(loggdp, 1),
        l_gdp_1   = lag(gdp_diff, 1),
        l_gdp_2   = lag(gdp_diff, 2),
        l_shock_1 = lag(maxwind_resid, 1),
        l_shock_2 = lag(maxwind_resid, 2)
      ) %>%
      ungroup() %>%
      filter(!is.na(dy), !is.na(l_gdp_1), !is.na(l_gdp_2)) %>%
      mutate(h = h)
  }) %>%
  mutate(
    fh         = factor(h),
    country_h  = interaction(countrycode, fh, drop = TRUE),
    # Polynomial shock terms (w̃ × h^k)
    shock_h0   = maxwind_resid,
    shock_h1   = maxwind_resid * h,
    shock_h2   = maxwind_resid * h^2,
    # Polynomial gap-interaction terms (w̃ × τ × h^k)
    shock_gap_h0 = maxwind_resid * years_since_any,
    shock_gap_h1 = maxwind_resid * years_since_any * h,
    shock_gap_h2 = maxwind_resid * years_since_any * h^2
  )
}

estimate_smooth <- function(stacked_df, label) {
  cat("  Smooth LP:", label, "| stacked rows =", nrow(stacked_df), "\n")
  sv <- c("shock_h0", "shock_h1", "shock_h2",
          "shock_gap_h0", "shock_gap_h1", "shock_gap_h2")

  fml <- as.formula(paste0(
    "dy ~ ", paste(sv, collapse = " + "),
    " + i(fh, l_gdp_1) + i(fh, l_gdp_2)",
    " + i(fh, l_shock_1) + i(fh, l_shock_2)",
    " | country_h[year] + country_h[year2] + year^fh"
  ))

  mod <- tryCatch(
    feols(fml, data = stacked_df, vcov = ~countrycode),
    error = function(e) { message("  Skip: ", e$message); NULL }
  )
  if (is.null(mod)) return(NULL)

  theta <- coef(mod)[sv]
  V     <- vcov(mod)[sv, sv]

  h_seq <- seq(0, HORIZON, by = 0.25)

  map_dfr(c(1, 3, 5), function(g) {
    map_dfr(h_seq, function(hh) {
      x       <- c(1, hh, hh^2, g, g * hh, g * hh^2)
      irf_val <- sum(x * theta)
      irf_se  <- as.numeric(sqrt(t(x) %*% V %*% x))
      tibble(
        horizon  = hh,
        gap      = g,
        irf_mean = irf_val * 100,   # to pp immediately
        se       = irf_se  * 100,
        irf_down = (irf_val - 1.96 * irf_se) * 100,
        irf_up   = (irf_val + 1.96 * irf_se) * 100,
        sample   = label
      )
    })
  })
}

cat("\n--- Smooth parametric LP by W\u0304 threshold ---\n")
irf_smooth_thresh <- map_dfr(seq_along(wbar_thresholds), function(i) {
  th  <- wbar_thresholds[i]
  lbl <- thresh_labels[i]
  d   <- filter(data, W_bar >= th)
  sk  <- build_stacked_id(d)
  estimate_smooth(sk, lbl)
}) %>%
  mutate(
    gap    = factor(gap, levels = c(1, 3, 5),
                    labels = c("1 yr (compounding)", "3 yrs", "5 yrs (isolated)")),
    sample = factor(sample, levels = thresh_labels)
  )

pal_smooth_gap <- c(
  "1 yr (compounding)" = myred,
  "3 yrs"              = mygold,
  "5 yrs (isolated)"   = myblue
)

p_smooth_thresh <- ggplot(irf_smooth_thresh,
                           aes(x = horizon, colour = gap, fill = gap)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.10, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(aes(y = irf_mean), linewidth = 1.4) +
  facet_wrap(~ sample, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, HORIZON, 2)) +
  scale_colour_manual(values = pal_smooth_gap) +
  scale_fill_manual(values   = pal_smooth_gap) +
  labs(
    title    = TeX("Smooth Compound LP by $\\bar{W}$ Threshold (Quadratic in h)"),
    subtitle = TeX("Stacked LP: ME(h, $\\tau$=g) = $\\beta$(h) + g$\\cdot$$\\gamma$(h), $\\beta$,$\\gamma$ quadratic in h.  Evaluated at $\\tau$ $\\in$ \\{1,3,5\\} yrs.  SE clustered by country."),
    x = "Years after strike", y = "Cumulative GDP growth (%)",
    colour = NULL, fill = NULL
  ) +
  theme_nature() +
  guides(colour = guide_legend(nrow = 1))

ggsave("figures/compound_smooth_wbar_thresh.png", p_smooth_thresh,
       width = 12, height = 10, dpi = 300)
message("Saved compound_smooth_wbar_thresh.png")

# ── Smooth gap coefficient γ(h) across thresholds ----------------------------
# γ(h) = θ₃ + θ₄h + θ₅h² is the interaction term: the extra GDP loss per
# additional year of τ at horizon h.  Positive γ(h) = isolated more damaging.

p_smooth_gamma <- irf_smooth_thresh %>%
  # Recover γ(h) = [ME(h,5) - ME(h,1)] / 4  (finite difference over 4 yrs)
  filter(gap %in% c("1 yr (compounding)", "5 yrs (isolated)")) %>%
  select(sample, horizon, gap, irf_mean, irf_down, irf_up) %>%
  pivot_wider(names_from = gap,
              values_from = c(irf_mean, irf_down, irf_up),
              names_repair = "unique") %>%
  mutate(
    gamma      = (`irf_mean_5 yrs (isolated)` - `irf_mean_1 yr (compounding)`) / 4,
    gamma_lo   = (`irf_down_5 yrs (isolated)` - `irf_up_1 yr (compounding)`)   / 4,
    gamma_hi   = (`irf_up_5 yrs (isolated)`   - `irf_down_1 yr (compounding)`) / 4
  ) %>%
  ggplot(aes(x = horizon, y = gamma, colour = sample, fill = sample)) +
  geom_ribbon(aes(ymin = gamma_lo, ymax = gamma_hi), alpha = 0.10, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 1.5) +
  scale_x_continuous(breaks = seq(0, HORIZON, 2)) +
  scale_colour_manual(values = thresh_pal,
                      labels = c(TeX("$\\bar{W}$ $\\geq$ 0 m/s"),
                                 TeX("$\\bar{W}$ $\\geq$ 10 m/s"),
                                 TeX("$\\bar{W}$ $\\geq$ 15 m/s"),
                                 TeX("$\\bar{W}$ $\\geq$ 20 m/s"))) +
  scale_fill_manual(values   = thresh_pal,
                    labels = c(TeX("$\\bar{W}$ $\\geq$ 0 m/s"),
                               TeX("$\\bar{W}$ $\\geq$ 10 m/s"),
                               TeX("$\\bar{W}$ $\\geq$ 15 m/s"),
                               TeX("$\\bar{W}$ $\\geq$ 20 m/s"))) +
  labs(
    title    = TeX("Gap-Interaction Coefficient $\\gamma$(h) by $\\bar{W}$ Threshold"),
    subtitle = TeX("[ME($\\tau$=5) - ME($\\tau$=1)] / 4 per yr \u2014 positive = isolated events more damaging"),
    x = "Years after strike", y = TeX("$\\gamma$(h) (pp per yr of $\\tau$)"),
    colour = NULL, fill = NULL
  ) +
  theme_nature() +
  guides(colour = guide_legend(nrow = 2))

ggsave("figures/compound_smooth_gamma.png", p_smooth_gamma,
       width = 9, height = 6, dpi = 300)
message("Saved compound_smooth_gamma.png")

# =============================================================================
# 4. Leave-one-continent-out
# =============================================================================

continents <- data %>%
  filter(!is.na(continent), maxwind > 0) %>%
  count(continent) %>%
  filter(n >= 20) %>%
  pull(continent)

cat("\n--- Leave-one-continent-out:", continents, "---\n")

irf_loco <- map_dfr(continents, function(ct) {
  run_lp(filter(data, continent != ct | is.na(continent)), paste0("Drop ", ct))
}) %>%
  mutate(quantile = factor(quantile, levels = eval_labels))

p_loco <- irf_ribbon(irf_loco, facet_var = "sample") +
  labs(title    = "Leave-One-Continent-Out",
       subtitle = TeX("Compound LP: $\\tau$ $\\in$ \\{1,3,5\\} yrs, dropping one continent at a time"))

ggsave("figures/compound_rob_loco.png", p_loco,
       width = 11, height = 8, dpi = 300)
message("Saved compound_rob_loco.png")

# =============================================================================
# 5. Within climatology groups — binary
# =============================================================================
# Interpretation: does the τ gradient persist within each W̄ bin?
# Limitation: within-group W̄ heterogeneity still exists, so this only
# partially controls for W̄ confounding.

cat("\n--- Within clim2 (binary) groups ---\n")
irf_clim2 <- bind_rows(
  run_lp(filter(data, clim2 == "Low"),  "Low (W\u0304 < 20 m/s)"),
  run_lp(filter(data, clim2 == "High"), "High (W\u0304 \u2265 20 m/s)")
) %>%
  mutate(quantile = factor(quantile, levels = eval_labels),
         sample   = factor(sample, levels = c("Low (W\u0304 < 20 m/s)",
                                              "High (W\u0304 \u2265 20 m/s)")))

p_clim2 <- irf_ribbon(irf_clim2, facet_var = "sample") +
  labs(title    = "Within Climatology Groups (Binary)",
       subtitle = TeX("Does the $\\tau$ gradient survive splitting at $\\bar{W}$ = 20 m/s?"))

ggsave("figures/compound_rob_clim_groups.png", p_clim2,
       width = 11, height = 6, dpi = 300)
message("Saved compound_rob_clim_groups.png")

# =============================================================================
# 6. Within climatology groups — terciles (finer W̄ bins)
# =============================================================================
# Finer bins reduce within-group W̄ spread, making the τ variation
# more orthogonal to W̄ level differences.

cat("\n--- Within clim3 (tercile) groups ---\n")
clim3_levels <- c("T1 (0\u201310 m/s)", "T2 (10\u201330 m/s)", "T3 (>30 m/s)")
irf_clim3 <- map_dfr(clim3_levels, function(g) {
    d <- filter(data, clim3 == g)
    if (nrow(d) < 100) { message("Skip group ", g, ": too few obs"); return(NULL) }
    run_lp(d, g)
  }
) %>%
  mutate(quantile = factor(quantile, levels = eval_labels),
         sample   = factor(sample, levels = clim3_levels))

if (nrow(irf_clim3) > 0) {
  p_clim3 <- irf_ribbon(irf_clim3, facet_var = "sample") +
    labs(title    = "Within Climatology Groups (Fixed Breaks: 10 and 30 m/s)",
         subtitle = "T1: 0\u201310 m/s | T2: 10\u201330 m/s | T3: >30 m/s")
  ggsave("figures/compound_rob_clim_terciles.png", p_clim3,
         width = 14, height = 6, dpi = 300)
  message("Saved compound_rob_clim_terciles.png")
}

# =============================================================================
# 7b. Country-specific percentile shock × compound: p50 / p75 / p90 / p95
# =============================================================================
# Motivation: absolute wind thresholds conflate intensity with climatology.
# Country-specific percentile thresholds define "extreme" relative to each
# country's own distribution of struck-year maxwind.
#
# Design:
#   For each percentile p ∈ {50, 75, 90, 95}:
#     shock_p_{it}  = 1{maxwind_{it} ≥ W̄_p_i}   (binary, country-specific)
#     τ_p_{it}      = years since last shock_p event
#
#   Regular LP with continuous interaction (horizon-by-horizon, DK SE):
#     Δloggdp_{i,t+h} ~ shock_p + shock_p × τ_p + controls | FE
#     ME(τ=c) = β₁ + β₂·c,  evaluated at c ∈ {3, 10} yrs
#
#   Faceted by percentile: 4 panels, two lines each (red = recent, blue = isolated)
# =============================================================================

cat("\n--- Country-specific percentile shock \u00d7 compound (p75/p90/p95/p99) ---\n")

pctile_grid <- tibble(
  pct     = c(0.75, 0.90, 0.95, 0.99),
  pct_lbl = c("p75", "p90", "p95", "p99")
)

eval_tau_pct  <- c(2, 10)
eval_labs_pct <- c("2 yrs (recent)", "10 yrs (isolated)")
pal_pct_tau   <- c("2 yrs (recent)"    = myred,
                   "10 yrs (isolated)" = myblue)

# ── Compute thresholds and run LP for each percentile -----------------------
irf_pctile_all <- map_dfr(seq_len(nrow(pctile_grid)), function(i) {
  pv  <- pctile_grid$pct[i]
  lbl <- pctile_grid$pct_lbl[i]

  # Country-specific threshold
  thresh_cty <- data %>%
    filter(W_bar > 0, maxwind > 0) %>%
    group_by(countrycode) %>%
    summarise(w_thresh = quantile(maxwind, pv, na.rm = TRUE), .groups = "drop")

  d_tmp <- data %>%
    left_join(thresh_cty, by = "countrycode") %>%
    mutate(
      shock_pct = case_when(
        is.na(w_thresh)       ~ 0L,
        maxwind >= w_thresh   ~ 1L,
        TRUE                  ~ 0L
      )
    ) %>%
    arrange(countrycode, year) %>%
    group_by(countrycode) %>%
    mutate(tau_pct = time_since_last(shock_pct, year)) %>%
    ungroup()

  n_ev  <- sum(d_tmp$shock_pct)
  r_cor <- round(cor(d_tmp$W_bar, d_tmp$tau_pct, use = "complete.obs"), 3)
  cat(sprintf("  %s: %d events, cor(W\u0304, \u03c4) = %.3f\n", lbl, n_ev, r_cor))

  tryCatch(
    lp_state_dep(
      data         = d_tmp,
      outcome      = "loggdp",
      main_var     = "shock_pct",
      state_var    = "tau_pct",
      controls     = "l(shock_pct, 1:2) + l(gdp_diff, 1:2)",
      horizon      = HORIZON,
      fe           = FE,
      panel_id     = PANEL,
      vcov_formula = DK,
      eval_values  = eval_tau_pct,
      eval_labels  = eval_labs_pct
    ) %>%
      mutate(
        pct_label = lbl,
        irf_mean  = irf_mean * 100,
        se        = se       * 100,
        irf_down  = irf_down * 100,
        irf_up    = irf_up   * 100
      ),
    error = function(e) { message("  Skip ", lbl, ": ", e$message); NULL }
  )
}) %>%
  mutate(
    quantile  = factor(quantile,  levels = eval_labs_pct),
    pct_label = factor(pct_label, levels = pctile_grid$pct_lbl)
  )

cat("\nGradient at h = 5 across percentiles:\n")
irf_pctile_all %>%
  filter(horizon == 5) %>%
  select(pct_label, quantile, irf_mean, se, irf_down, irf_up) %>%
  mutate(across(where(is.numeric), ~round(., 2))) %>%
  as.data.frame() %>% print()

# ── 4-panel IRF figure (one panel per percentile, free y) -------------------
p_p90_irf <- ggplot(irf_pctile_all,
                    aes(x = horizon, colour = quantile, fill = quantile)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(aes(y = irf_mean), linewidth = 1.6) +
  facet_wrap(~ pct_label, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, HORIZON, 2)) +
  scale_colour_manual(values = pal_pct_tau) +
  scale_fill_manual(values   = pal_pct_tau) +
  labs(
    title    = TeX("Compound LP: Country-Specific Percentile Shock $\\times$ $\\tau$"),
    subtitle = TeX("Shock = 1\\{maxwind $\\geq$ $\\bar{W}_{p,i}$\\} for p $\\in$ \\{75, 90, 95, 99\\}.  Regular LP (DK SE).  $\\tau$ evaluated at \\{2, 10\\} yrs.  Full panel."),
    x = "Years after strike", y = "Cumulative GDP growth (%)",
    colour = NULL, fill = NULL
  ) +
  theme_nature() +
  guides(colour = guide_legend(nrow = 1))

# ── Diagnostic: τ vs W̄ for each percentile (strip) -------------------------
p_tau_p90_wbar <- irf_pctile_all   # reuse name for combined plot below; skip separate diagnostic

ggsave("figures/compound_p90_irf.png", p_p90_irf,
       width = 12, height = 10, dpi = 300)
message("Saved compound_p90_irf.png")

# ── Gradient at h=5 as a function of percentile (summary) -------------------
grad_by_pct <- irf_pctile_all %>%
  filter(horizon == 5) %>%
  select(pct_label, quantile, irf_mean, se) %>%
  group_by(pct_label) %>%
  summarise(
    grad    = diff(irf_mean),          # ME(10yr) - ME(3yr)
    grad_se = sqrt(sum(se^2)),         # approx SE
    .groups = "drop"
  ) %>%
  mutate(
    grad_lo = grad - 1.96 * grad_se,
    grad_hi = grad + 1.96 * grad_se
  )

cat("\nGradient [ME(\u03c4=10) \u2212 ME(\u03c4=3)] at h=5 by percentile:\n")
grad_by_pct %>% mutate(across(where(is.numeric), ~round(., 2))) %>%
  as.data.frame() %>% print()

p_p90_combined <- p_p90_irf   # single figure; no separate diagnostic

ggsave("figures/compound_p90_combined.png", p_p90_irf,
       width = 12, height = 10, dpi = 300)
message("Saved compound_p90_combined.png")

irf_p90 <- irf_pctile_all   # for saveRDS compatibility

# =============================================================================
# 7c. Smooth parametric LP version of the percentile IRF figure
# =============================================================================
# Same design as 7b but replaces horizon-by-horizon LP with the stacked
# quadratic smooth LP: ME(h, τ=g) = β(h) + g·γ(h), β,γ quadratic in h.
# Evaluated at τ ∈ {2, 10} yrs.  SE clustered by country.
# =============================================================================

cat("\n--- Smooth LP: country-specific percentile shock (p75/p90/p95/p99) ---\n")

build_stacked_pct <- function(dat, H = 10) {
  map_dfr(0:H, function(h) {
    dat %>%
      group_by(countrycode) %>%
      arrange(year) %>%
      mutate(
        dy        = lead(loggdp, h) - lag(loggdp, 1),
        l_gdp_1   = lag(gdp_diff, 1),
        l_gdp_2   = lag(gdp_diff, 2),
        l_shock_1 = lag(shock_pct, 1),
        l_shock_2 = lag(shock_pct, 2)
      ) %>%
      ungroup() %>%
      filter(!is.na(dy), !is.na(l_gdp_1), !is.na(l_gdp_2),
             !is.na(shock_pct), !is.na(tau_pct)) %>%
      mutate(h = h)
  }) %>%
  mutate(
    fh           = factor(h),
    country_h    = interaction(countrycode, fh, drop = TRUE),
    shock_h0     = shock_pct,
    shock_h1     = shock_pct * h,
    shock_h2     = shock_pct * h^2,
    shock_gap_h0 = shock_pct * tau_pct,
    shock_gap_h1 = shock_pct * tau_pct * h,
    shock_gap_h2 = shock_pct * tau_pct * h^2
  )
}

estimate_smooth_pct <- function(stacked_df, label,
                                eval_g   = eval_tau_pct,
                                eval_lbl = eval_labs_pct) {
  cat("  Smooth LP:", label, "| stacked rows =", nrow(stacked_df), "\n")
  sv <- c("shock_h0", "shock_h1", "shock_h2",
          "shock_gap_h0", "shock_gap_h1", "shock_gap_h2")

  fml <- as.formula(paste0(
    "dy ~ ", paste(sv, collapse = " + "),
    " + i(fh, l_gdp_1) + i(fh, l_gdp_2)",
    " + i(fh, l_shock_1) + i(fh, l_shock_2)",
    " | country_h[year] + country_h[year2] + year^fh"
  ))

  mod <- tryCatch(
    feols(fml, data = stacked_df, vcov = ~countrycode),
    error = function(e) { message("  Skip: ", e$message); NULL }
  )
  if (is.null(mod)) return(NULL)

  theta <- coef(mod)[sv]
  V     <- vcov(mod)[sv, sv]
  h_seq <- seq(0, HORIZON, by = 0.25)

  map_dfr(seq_along(eval_g), function(j) {
    g   <- eval_g[j]
    lbl <- eval_lbl[j]
    map_dfr(h_seq, function(hh) {
      x       <- c(1, hh, hh^2, g, g * hh, g * hh^2)
      irf_val <- sum(x * theta)
      irf_se  <- as.numeric(sqrt(t(x) %*% V %*% x))
      tibble(
        horizon  = hh,
        quantile = lbl,
        irf_mean = irf_val * 100,
        se       = irf_se  * 100,
        irf_down = (irf_val - 1.96 * irf_se) * 100,
        irf_up   = (irf_val + 1.96 * irf_se) * 100,
        sample   = label
      )
    })
  })
}

irf_smooth_pct <- map_dfr(seq_len(nrow(pctile_grid)), function(i) {
  pv  <- pctile_grid$pct[i]
  lbl <- pctile_grid$pct_lbl[i]

  thresh_cty <- data %>%
    filter(W_bar > 0, maxwind > 0) %>%
    group_by(countrycode) %>%
    summarise(w_thresh = quantile(maxwind, pv, na.rm = TRUE), .groups = "drop")

  d_tmp <- data %>%
    left_join(thresh_cty, by = "countrycode") %>%
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

  sk <- build_stacked_pct(d_tmp)
  estimate_smooth_pct(sk, lbl)
}) %>%
  mutate(
    quantile  = factor(quantile,  levels = eval_labs_pct),
    pct_label = factor(sample,    levels = pctile_grid$pct_lbl)
  )

p_smooth_pct <- ggplot(irf_smooth_pct,
                       aes(x = horizon, colour = quantile, fill = quantile)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(aes(y = irf_mean), linewidth = 1.6) +
  facet_wrap(~ pct_label, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, HORIZON, 2)) +
  scale_colour_manual(values = pal_pct_tau) +
  scale_fill_manual(values   = pal_pct_tau) +
  labs(
    title    = TeX("Smooth Compound LP: Country-Specific Percentile Shock $\\times$ $\\tau$"),
    subtitle = TeX("Stacked LP: ME(h, $\\tau$=g) = $\\beta$(h) + g$\\cdot$$\\gamma$(h), quadratic in h.  Shock = 1\\{maxwind $\\geq$ $\\bar{W}_{p,i}$\\} for p $\\in$ \\{75, 90, 95, 99\\}.  $\\tau$ $\\in$ \\{2, 10\\} yrs.  SE clustered by country."),
    x = "Years after strike", y = "Cumulative GDP growth (%)",
    colour = NULL, fill = NULL
  ) +
  theme_nature() +
  guides(colour = guide_legend(nrow = 1))

ggsave("figures/compound_p90_smooth_irf.png", p_smooth_pct,
       width = 12, height = 10, dpi = 300)
message("Saved compound_p90_smooth_irf.png")

# =============================================================================
# 7d. p95 sensitivity to evaluation point pairs: regular LP, quadratic, cubic
# =============================================================================
# Three versions of the same figure, each with four panels for the pairs:
#   {2,10}  {1,5}  {2,15}  {5,20}
# Red = short τ, blue = long τ (exact values shown in panel strip).
# =============================================================================

cat("\n--- p95 eval-point sensitivity (regular LP + quadratic + cubic) ---\n")

# ── Shared data ---------------------------------------------------------------
thresh_p95 <- data %>%
  filter(W_bar > 0, maxwind > 0) %>%
  group_by(countrycode) %>%
  summarise(w_thresh = quantile(maxwind, 0.95, na.rm = TRUE), .groups = "drop")

d_p95 <- data %>%
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

eval_pairs <- list(
  list(vals = c(15, 20), panel = "15 vs 20 yrs since last strike"),
  list(vals = c(10, 20), panel = "10 vs 20 yrs since last strike"),
  list(vals = c(5,  20), panel = "5 vs 20 yrs since last strike"),
  list(vals = c(1,  20), panel = "1 vs 20 yrs since last strike")
)
panel_levels <- map_chr(eval_pairs, "panel")

# ── Helper: extract IRF curves from a fitted model at one eval pair -----------
eval_pair_irf <- function(theta, V, ep, h_seq, degree) {
  map_dfr(seq_along(ep$vals), function(j) {
    g    <- ep$vals[j]
    type <- if (j == 1) "Recent strike" else "Long recovery"
    map_dfr(h_seq, function(hh) {
      x <- switch(as.character(degree),
        "2" = c(1, hh, hh^2, g, g*hh, g*hh^2),
        "3" = c(1, hh, hh^2, hh^3, g, g*hh, g*hh^2, g*hh^3)
      )
      irf_val <- sum(x * theta)
      irf_se  <- as.numeric(sqrt(t(x) %*% V %*% x))
      tibble(
        horizon  = hh,
        type     = type,
        tau_val  = g,
        irf_mean = irf_val * 100,
        irf_down = (irf_val - 1.96 * irf_se) * 100,
        irf_up   = (irf_val + 1.96 * irf_se) * 100,
        panel    = ep$panel
      )
    })
  })
}

pal_type <- c("Recent strike" = myred, "Long recovery" = myblue)
h_seq    <- seq(0, HORIZON, by = 0.25)

make_sens_plot <- function(irf_df, title_str, subtitle_str) {
  ggplot(irf_df, aes(x = horizon, colour = type, fill = type)) +
    geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_line(aes(y = irf_mean), linewidth = 1.6) +
    # label each line with its τ value at the right edge
    geom_text(
      data = irf_df %>%
        group_by(panel, type) %>%
        filter(horizon == max(horizon)),
      aes(x = horizon, y = irf_mean,
          label = paste0("tau=", tau_val)),
      hjust = -0.15, size = 3.2, show.legend = FALSE
    ) +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, HORIZON, 2),
                       expand = expansion(mult = c(0.02, 0.12))) +
    scale_colour_manual(
      values = pal_type,
      labels = c("Recent strike", "Long recovery")
    ) +
    scale_fill_manual(values = pal_type, guide = "none") +
    labs(title = title_str, subtitle = subtitle_str,
         x = "Years after strike", y = "Cumulative GDP growth (%)",
         colour = NULL) +
    theme_nature() +
    guides(colour = guide_legend(nrow = 1,
                                 override.aes = list(linewidth = 2,
                                                     fill = NA)))
}

# ── (A) Regular LP -----------------------------------------------------------
cat("  Fitting regular LP for all eval points...\n")
all_vals <- unique(unlist(map(eval_pairs, "vals")))
irf_reg_all <- lp_state_dep(
  data         = d_p95,
  outcome      = "loggdp",
  main_var     = "shock_pct",
  state_var    = "tau_pct",
  controls     = "l(shock_pct, 1:2) + l(gdp_diff, 1:2)",
  horizon      = HORIZON,
  fe           = FE,
  panel_id     = PANEL,
  vcov_formula = DK,
  eval_values  = sort(all_vals),
  eval_labels  = as.character(sort(all_vals))
)

irf_reg_pairs <- map_dfr(eval_pairs, function(ep) {
  bind_rows(
    irf_reg_all %>%
      filter(quantile == as.character(ep$vals[1])) %>%
      mutate(type = "Recent strike", tau_val = ep$vals[1], panel = ep$panel),
    irf_reg_all %>%
      filter(quantile == as.character(ep$vals[2])) %>%
      mutate(type = "Long recovery", tau_val = ep$vals[2], panel = ep$panel)
  )
}) %>%
  mutate(
    irf_mean = irf_mean * 100,
    irf_down = irf_down * 100,
    irf_up   = irf_up   * 100,
    panel    = factor(panel, levels = panel_levels),
    type     = factor(type,  levels = c("Recent strike", "Long recovery"))
  )

p_reg_sens <- make_sens_plot(
  irf_reg_pairs,
  "Recovery Time Matters: Standard LP (fixed long gap = 20 yrs)",
  "Four short-gap values vs fixed long gap of 20 yrs. Shock = country p95 wind indicator. DK standard errors."
)
ggsave("figures/compound_p95_reg_sensitivity.png", p_reg_sens,
       width = 12, height = 10, dpi = 300)
message("Saved compound_p95_reg_sensitivity.png")

# ── (B) Quadratic smooth LP --------------------------------------------------
sv2 <- c("shock_h0","shock_h1","shock_h2",
         "shock_gap_h0","shock_gap_h1","shock_gap_h2")
sk_p95 <- build_stacked_pct(d_p95)
cat("  Fitting quadratic smooth LP | stacked rows =", nrow(sk_p95), "\n")
mod_q <- feols(
  as.formula(paste0(
    "dy ~ ", paste(sv2, collapse = " + "),
    " + i(fh, l_gdp_1) + i(fh, l_gdp_2) + i(fh, l_shock_1) + i(fh, l_shock_2)",
    " | country_h[year] + country_h[year2] + year^fh"
  )),
  data = sk_p95, vcov = ~countrycode
)
theta_q <- coef(mod_q)[sv2]
V_q     <- vcov(mod_q)[sv2, sv2]

irf_quad_pairs <- map_dfr(eval_pairs, function(ep) {
  eval_pair_irf(theta_q, V_q, ep, h_seq, degree = 2) %>%
    mutate(panel = factor(panel, levels = panel_levels),
           type  = factor(type,  levels = c("Recent strike", "Long recovery")))
})
p_quad_sens <- make_sens_plot(
  irf_quad_pairs,
  "Recovery Time Matters: Quadratic Smooth LP (fixed long gap = 20 yrs)",
  "Stacked LP with quadratic polynomial in horizon (6 params). Shock = country p95 wind indicator. SE clustered by country."
)
ggsave("figures/compound_p95_eval_sensitivity.png", p_quad_sens,
       width = 12, height = 10, dpi = 300)
message("Saved compound_p95_eval_sensitivity.png")

# ── (C) Cubic smooth LP ------------------------------------------------------
# Extend stacked dataset with h^3 polynomial terms
sk_p95_c <- sk_p95 %>%
  mutate(
    shock_h3     = shock_pct * h^3,
    shock_gap_h3 = shock_pct * tau_pct * h^3
  )
sv3 <- c("shock_h0", "shock_h1", "shock_h2", "shock_h3",
         "shock_gap_h0", "shock_gap_h1", "shock_gap_h2", "shock_gap_h3")
cat("  Fitting cubic smooth LP | stacked rows =", nrow(sk_p95_c), "\n")
mod_c <- feols(
  as.formula(paste0(
    "dy ~ ", paste(sv3, collapse = " + "),
    " + i(fh, l_gdp_1) + i(fh, l_gdp_2) + i(fh, l_shock_1) + i(fh, l_shock_2)",
    " | country_h[year] + country_h[year2] + year^fh"
  )),
  data = sk_p95_c, vcov = ~countrycode
)
theta_c <- coef(mod_c)[sv3]
V_c     <- vcov(mod_c)[sv3, sv3]

irf_cubic_pairs <- map_dfr(eval_pairs, function(ep) {
  eval_pair_irf(theta_c, V_c, ep, h_seq, degree = 3) %>%
    mutate(panel = factor(panel, levels = panel_levels),
           type  = factor(type,  levels = c("Recent strike", "Long recovery")))
})
p_cubic_sens <- make_sens_plot(
  irf_cubic_pairs,
  "Recovery Time Matters: Cubic Smooth LP (fixed long gap = 20 yrs)",
  "Stacked LP with cubic polynomial in horizon (8 params). Shock = country p95 wind indicator. SE clustered by country."
)
ggsave("figures/compound_p95_cubic_sensitivity.png", p_cubic_sens,
       width = 12, height = 10, dpi = 300)
message("Saved compound_p95_cubic_sensitivity.png")

# =============================================================================
# 7e. compound_p95_irf: fixed long τ = 20, short τ ∈ {10, 5, 3, 1}, shared y
# =============================================================================

cat("\n--- compound_p95_irf (long=20, short varies, shared y) ---\n")

pal_fixed <- c("Recent strike" = myred, "Long recovery" = myblue)

irf_pairs_fixed <- list(
  list(vals = c(15, 20), panel = "15 vs 20 yrs since last strike"),
  list(vals = c(10, 20), panel = "10 vs 20 yrs since last strike"),
  list(vals = c( 5, 20), panel = "5 vs 20 yrs since last strike"),
  list(vals = c( 1, 20), panel = "1 vs 20 yrs since last strike")
)
panel_levels_fixed <- map_chr(irf_pairs_fixed, "panel")

# ── (A) Regular LP -----------------------------------------------------------
all_vals_fixed <- c(1, 5, 10, 15, 20)
irf_reg_fixed_all <- lp_state_dep(
  data         = d_p95,
  outcome      = "loggdp",
  main_var     = "shock_pct",
  state_var    = "tau_pct",
  controls     = "l(shock_pct, 1:2) + l(gdp_diff, 1:2)",
  horizon      = HORIZON,
  fe           = FE,
  panel_id     = PANEL,
  vcov_formula = DK,
  eval_values  = all_vals_fixed,
  eval_labels  = as.character(all_vals_fixed)
)

irf_reg_fixed <- map_dfr(irf_pairs_fixed, function(ep) {
  bind_rows(
    irf_reg_fixed_all %>%
      filter(quantile == as.character(ep$vals[1])) %>%
      mutate(type = "Recent strike", tau_val = ep$vals[1],
             panel = as.character(ep$panel)),
    irf_reg_fixed_all %>%
      filter(quantile == as.character(ep$vals[2])) %>%
      mutate(type = "Long recovery", tau_val = ep$vals[2],
             panel = as.character(ep$panel))
  )
}) %>%
  mutate(
    irf_mean = irf_mean * 100,
    irf_down = irf_down * 100,
    irf_up   = irf_up   * 100,
    panel    = factor(panel, levels = panel_levels_fixed),
    type     = factor(type,  levels = c("Recent strike", "Long recovery"))
  )

p_reg_fixed <- ggplot(irf_reg_fixed,
                      aes(x = horizon, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  geom_line(aes(y = irf_mean), linewidth = 1.6) +
  geom_text(
    data = irf_reg_fixed %>%
      group_by(panel, type) %>% filter(horizon == max(horizon)),
    aes(x = horizon, y = irf_mean,
        label = paste0(tau_val, " yrs")),
    hjust = -0.15, size = 3.2, show.legend = FALSE
  ) +
  facet_wrap(~ panel, ncol = 2, scales = "fixed") +
  scale_x_continuous(breaks = seq(0, HORIZON, 2),
                     expand = expansion(mult = c(0.02, 0.12))) +
  scale_colour_manual(values = pal_fixed,
                      labels = c("Recent strike", "Long recovery")) +
  scale_fill_manual(values = pal_fixed, guide = "none") +
  labs(
    title    = "GDP Response After a Major Cyclone: Does Recovery Time Matter?",
    subtitle = "Each panel compares countries hit recently (red) vs. after a long recovery window (blue).\nShock defined at country-specific 95th percentile of wind speed. Shared y-axis.",
    x = "Years after strike", y = "Cumulative GDP growth (%)",
    colour = NULL
  ) +
  theme_nature() +
  guides(colour = guide_legend(nrow = 1,
                               override.aes = list(linewidth = 2, fill = NA)))

ggsave("figures/compound_p95_irf.png", p_reg_fixed,
       width = 12, height = 10, dpi = 300)
message("Saved compound_p95_irf.png")

# ── (B) Quadratic smooth LP --------------------------------------------------
irf_quad_fixed <- map_dfr(irf_pairs_fixed, function(ep) {
  eval_pair_irf(theta_q, V_q, ep, h_seq, degree = 2) %>%
    mutate(
      panel = factor(panel, levels = panel_levels_fixed),
      type  = factor(type,  levels = c("Recent strike", "Long recovery"))
    )
})

p_quad_fixed <- ggplot(irf_quad_fixed,
                       aes(x = horizon, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  geom_line(aes(y = irf_mean), linewidth = 1.6) +
  facet_wrap(~ panel, ncol = 2, scales = "fixed") +
  scale_x_continuous(breaks = seq(0, HORIZON, 2)) +
  scale_colour_manual(values = pal_fixed) +
  scale_fill_manual(values = pal_fixed, guide = "none") +
  labs(
    title    = "GDP Response After a Major Cyclone: Does Recovery Time Matter? (Smooth LP)",
    subtitle = "Each panel compares countries hit recently (red) vs. after a long recovery window (blue).\nSmooth LP with quadratic time trend. Shock at country-specific 95th percentile. Shared y-axis.",
    x = "Years after strike", y = "Cumulative GDP growth (%)",
    colour = NULL
  ) +
  theme_nature() +
  guides(colour = guide_legend(nrow = 1,
                               override.aes = list(linewidth = 2, fill = NA)))

ggsave("figures/compound_p95_irf_smooth.png", p_quad_fixed,
       width = 12, height = 10, dpi = 300)
message("Saved compound_p95_irf_smooth.png")

# =============================================================================
# Section 7f: Formal gradient significance test at h = 5
# =============================================================================
# Gradient = ME(tau=20) - ME(tau=c) = beta2 * (20 - c)
# SE       = |20 - c| * SE(beta2)   [delta method; beta1 cancels in difference]
# Run h=5 regression directly on d_p95 to extract beta2 and its DK SE.
# =============================================================================

cat("\n--- Gradient significance test at h=5 ---\n")

d_h5 <- d_p95 %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  mutate(dy5 = lead(loggdp, 5) - lag(loggdp, 1)) %>%
  ungroup() %>%
  filter(!is.na(dy5), !is.na(shock_pct), !is.na(tau_pct))

mod_h5 <- feols(
  dy5 ~ shock_pct + shock_pct:tau_pct + l(shock_pct, 1:2) + l(gdp_diff, 1:2) |
    countrycode[year] + countrycode[year2] + region^year,
  data     = d_h5,
  panel.id = PANEL,
  vcov     = DK
)

b2    <- coef(mod_h5)["shock_pct:tau_pct"]
se_b2 <- sqrt(vcov(mod_h5)["shock_pct:tau_pct", "shock_pct:tau_pct"])
cat("beta2 =", round(b2 * 100, 3), "pp/yr  SE =", round(se_b2 * 100, 3), "\n")

grad_test <- tibble(
  tau_short = c(15, 10, 5, 1),
  label     = paste0(tau_short, " vs 20 yrs")
) %>%
  mutate(
    gap      = 20 - tau_short,
    gradient = b2   * gap * 100,
    grad_se  = se_b2 * gap * 100,
    # 90 % CI
    lo90     = gradient - 1.645 * grad_se,
    hi90     = gradient + 1.645 * grad_se,
    # 95 % CI
    lo95     = gradient - 1.960 * grad_se,
    hi95     = gradient + 1.960 * grad_se,
    t_stat   = gradient / grad_se,
    p_val    = round(2 * pnorm(-abs(t_stat)), 3),
    sig      = case_when(p_val < .01 ~ "***", p_val < .05 ~ "**",
                         p_val < .10 ~ "*",   TRUE ~ "n.s.")
  )

cat("Gradient test results:\n")
print(grad_test %>%
  select(label, gradient, grad_se, lo90, hi90, lo95, hi95, p_val, sig) %>%
  mutate(across(c(gradient, grad_se, lo90, hi90, lo95, hi95), ~ round(., 2))))

# Nature-style colour: deep teal-blue for significant rows, muted grey otherwise
col_sig   <- "#0D6B8C"   # deep teal, significant at 10%
col_insig <- "#9DAFBA"   # muted blue-grey, for reference

y_levels <- paste0(c(15, 10, 5, 1), " vs 20 yrs")

p_grad <- ggplot(grad_test,
                 aes(y = fct_rev(factor(label, levels = y_levels)))) +
  # Null line
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.55) +
  # 95 % CI — thinner outer whisker
  geom_errorbarh(aes(xmin = lo95, xmax = hi95),
                 height = 0, colour = col_sig, linewidth = 0.75, alpha = 0.55) +
  # 90 % CI — thicker inner whisker
  geom_errorbarh(aes(xmin = lo90, xmax = hi90),
                 height = 0, colour = col_sig, linewidth = 2.0) +
  # Point estimate
  geom_point(aes(x = gradient),
             size = 4.5, colour = col_sig, shape = 21,
             fill = "white", stroke = 1.6) +
  # Annotation
  geom_text(aes(x = hi95,
                label = paste0(" ", round(gradient, 1), " pp  [", sig, "]")),
            hjust = 0, size = 3.7, colour = "grey25") +
  scale_x_continuous(
    name   = "Gradient at year 5: long-recovery minus recent-strike (pp GDP)",
    expand = expansion(mult = c(0.05, 0.25)),
    breaks = scales::pretty_breaks(5)
  ) +
  labs(
    title    = "Recovery Time Matters: Recent Strikes Are More Damaging",
    subtitle = paste0(
      "Countries struck 20 years ago suffer smaller GDP losses than those struck recently.\n",
      "Country-specific p95 threshold. Thick bar = 90% CI, thin bar = 95% CI (DK SE).\n",
      "Significant at the 10% level for all comparisons (p = 0.069)."
    ),
    y = "Recovery gap compared to 20 yrs"
  ) +
  theme_nature(base_size = 14) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

ggsave("figures/compound_p95_gradient_test.png", p_grad,
       width = 9, height = 5, dpi = 300)
message("Saved compound_p95_gradient_test.png")

# =============================================================================
# 7g. Categorical τ bins — shock × bin interaction in a single LP
# =============================================================================
# Analogue of lp_state_dep but with a categorical state variable.
# Formula per horizon:
#   f(loggdp,h) - l(loggdp,1) ~ shock_pct:i(tau_bin) + controls | FE
#
# Using i(tau_bin) without a reference level means fixest omits the
# standalone shock_pct term and estimates one coefficient per bin directly —
# each coefficient IS the ME for that bin (no delta method needed, no ref).
# This mirrors how lp_state_dep recovers ME = β_main + β_inter×c, except
# the "state" is categorical so each bin gets its own β.
#
# Fixed bins (years since last country-specific p95 event):
#   Bin 1 : 1–5 yrs
#   Bin 2 : 5–10 yrs
#   Bin 3 : 10–15 yrs
#   Bin 4 : 15–20 yrs
# =============================================================================

cat("\n--- Categorical τ bins: shock × bin interaction LP (all 4 bins shown) ---\n")

bin_labels <- c("1\u20135 yrs", "5\u201310 yrs", "10\u201315 yrs", "15\u201320 yrs")
pal_bins   <- c("1\u20135 yrs"   = myred,
                "5\u201310 yrs"  = mygold,
                "10\u201315 yrs" = mygreen,
                "15\u201320 yrs" = myblue)

d_p95_cat <- d_p95 %>%
  mutate(
    tau_bin = case_when(
      is.na(tau_pct)                 ~ NA_character_,
      tau_pct >= 1  & tau_pct <= 5   ~ "1\u20135 yrs",
      tau_pct >  5  & tau_pct <= 10  ~ "5\u201310 yrs",
      tau_pct >  10 & tau_pct <= 15  ~ "10\u201315 yrs",
      tau_pct >  15 & tau_pct <= 20  ~ "15\u201320 yrs",
      TRUE                           ~ NA_character_
    ),
    tau_bin = factor(tau_bin, levels = bin_labels)
  )

cat("  Events per bin:\n")
print(table(d_p95_cat$tau_bin[d_p95_cat$shock_pct == 1], useNA = "ifany"))

# One LP per horizon: explicit shock_pct × bin dummies — each coefficient is
# directly the ME for that bin (same logic as lp_state_dep but categorical).
# Using explicit dummies avoids fixest misinterpreting i() as a varying slope
# when the FE spec contains countrycode[year] slope terms.
d_p95_cat <- d_p95_cat %>%
  mutate(
    s_b1 = shock_pct * (tau_bin == "1\u20135 yrs"   & !is.na(tau_bin)),
    s_b2 = shock_pct * (tau_bin == "5\u201310 yrs"  & !is.na(tau_bin)),
    s_b3 = shock_pct * (tau_bin == "10\u201315 yrs" & !is.na(tau_bin)),
    s_b4 = shock_pct * (tau_bin == "15\u201320 yrs" & !is.na(tau_bin))
  )

irf_bins <- map_dfr(0:HORIZON, function(h) {
  tryCatch({
    mod <- feols(
      as.formula(paste0(
        "f(loggdp, ", h, ") - l(loggdp, 1) ~ ",
        "s_b1 + s_b2 + s_b3 + s_b4 + ",
        "l(shock_pct, 1:2) + l(gdp_diff, 1:2) | ", FE
      )),
      data = d_p95_cat, panel.id = PANEL, vcov = DK
    )
    V  <- vcov(mod)
    cf <- coef(mod)
    map_dfr(
      list(c("s_b1", "1\u20135 yrs"), c("s_b2", "5\u201310 yrs"),
           c("s_b3", "10\u201315 yrs"), c("s_b4", "15\u201320 yrs")),
      function(x) {
        nm <- x[1]; bl <- x[2]
        if (!(nm %in% names(cf))) return(NULL)
        est <- cf[nm]; se <- sqrt(V[nm, nm])
        tibble(horizon  = h, tau_bin  = bl,
               irf_mean = est * 100,
               irf_down = (est - 1.96 * se) * 100,
               irf_up   = (est + 1.96 * se) * 100)
      })
  }, error = function(e) { message("h=", h, ": ", e$message); NULL })
}) %>%
  mutate(tau_bin = factor(tau_bin, levels = bin_labels))

p_bins <- ggplot(irf_bins,
                 aes(x = horizon, colour = tau_bin, fill = tau_bin)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.10, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(aes(y = irf_mean), linewidth = 1.6) +
  scale_x_continuous(breaks = seq(0, HORIZON, 2)) +
  scale_colour_manual(values = pal_bins, name = NULL) +
  scale_fill_manual(values   = pal_bins, guide = "none") +
  labs(
    title    = "GDP Response by Recovery Gap at Strike Time",
    subtitle = paste0(
      "Single LP per horizon: shock_pct \u00d7 \u03c4-bin categorical interaction (95% CI, DK SE).\n",
      "Each coefficient = ME for that bin. Country-specific p95 threshold."
    ),
    x = "Years after strike",
    y = "Cumulative GDP change (pp)"
  ) +
  theme_nature() +
  guides(colour = guide_legend(nrow = 1,
                               override.aes = list(linewidth = 2, fill = NA)))

ggsave("figures/compound_p95_irf_bins.png", p_bins,
       width = 8, height = 6.5, dpi = 300)
message("Saved compound_p95_irf_bins.png")

# ── Save ---------------------------------------------------------------------

saveRDS(
  list(irf_samples        = irf_samples,
       irf_wbar_thresh    = irf_wbar_thresh,
       irf_smooth_thresh  = irf_smooth_thresh,
       irf_loco           = irf_loco,
       irf_clim2          = irf_clim2,
       irf_clim3          = irf_clim3,
       irf_p90            = irf_p90,
       irf_smooth_pct     = irf_smooth_pct,
       irf_reg_pairs      = irf_reg_pairs,
       irf_quad_pairs     = irf_quad_pairs,
       irf_cubic_pairs    = irf_cubic_pairs,
       irf_reg_fixed      = irf_reg_fixed,
       irf_quad_fixed     = irf_quad_fixed,
       grad_test          = grad_test,
       irf_bins           = irf_bins),
  "results_compound_robustness.rds"
)
message("Saved results_compound_robustness.rds")
cat("\nDone.\n")
