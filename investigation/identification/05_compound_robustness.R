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
# 7. Residualised τ: project years_since_any on W̄, use residual as state var
#    — this is the direct non-parametric control for the W̄–τ correlation
# 8. Gradient summary: ME(5yr)−ME(1yr) at h=5 across all cuts
# =============================================================================

library(tidyverse)
library(fixest)
library(patchwork)
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

# ── Residualise years_since_any on W̄ ------------------------------------------
# This is the key remedy: remove the cross-country W̄–τ correlation.
# Within-country year-to-year variation in τ is preserved; only the
# systematic "high-W̄ → shorter τ" level difference is removed.
tau_lm <- lm(years_since_any ~ W_bar,
             data = filter(data, W_bar > 0, !is.na(years_since_any)))
cat("\nτ ~ W̄ regression (cross-country confound):\n")
print(summary(tau_lm)$coefficients)

data$tau_resid <- NA_real_
idx <- data$W_bar > 0 & !is.na(data$years_since_any)
data$tau_resid[idx] <- residuals(lm(years_since_any ~ W_bar,
                                    data = data[idx, ]))

cat("\nτ residual (after removing W̄ trend) quantiles:\n")
print(round(quantile(data$tau_resid, c(.1,.25,.5,.75,.9), na.rm = TRUE), 2))
cat("Correlation(years_since_any, tau_resid):",
    round(cor(data$years_since_any, data$tau_resid, use = "complete.obs"), 3), "\n")
cat("Correlation(W_bar,          tau_resid):",
    round(cor(data$W_bar, data$tau_resid, use = "complete.obs"), 3), "\n\n")

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
  scale_fill_manual(values = c(myred, mygold, myblue)) +
  labs(title = "Composition of Strike Events by Continent",
       subtitle = "Share compound / medium / isolated",
       x = NULL, y = "% of strikes", fill = NULL) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_desc_wbar <- compound_desc %>%
  ggplot(aes(x = W_bar, fill = compound_type, colour = compound_type)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  scale_fill_manual(values = c(myred, mygold, myblue)) +
  scale_colour_manual(values = c(myred, mygold, myblue)) +
  labs(title = "Country Climatology Distribution by Event Type",
       subtitle = "Compound events cluster in high-W\u0304 countries",
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
    title    = "\u03c4 Distribution Within Climatology Groups",
    subtitle = "Log scale. Both groups have substantial \u03c4 variation; W\u0304 confound not fully removed by binary split",
    x = "Climatology group (W\u0304 \u2265 20 m/s = High)",
    y = "Years since last TC (\u03c4, log scale)"
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
    title    = "\u03c4 vs W\u0304: Confound Persists Within Binary Groups",
    subtitle = "Within Low (blue) and High (red), higher W\u0304 still predicts shorter \u03c4",
    x = expression(bar(W)[i]~"(m/s)"),
    y = "Years since last TC (\u03c4, log scale)",
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
       subtitle = "Compound gradient with \u03c4 \u2208 {1,3,5} yrs")

ggsave("figures/compound_rob_samples.png", p_samples,
       width = 11, height = 6, dpi = 300)
message("Saved compound_rob_samples.png")

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

loco_h5 <- irf_loco %>%
  filter(horizon == 5) %>%
  select(sample, quantile, irf_mean, irf_down, irf_up) %>%
  mutate(across(c(irf_mean, irf_down, irf_up), ~ . * 100))

cat("\nLOCO gradient at h=5:\n")
loco_h5 %>% mutate(across(where(is.numeric), ~round(., 2))) %>%
  as.data.frame() %>% print()

p_loco <- ggplot(loco_h5,
                 aes(x = quantile, y = irf_mean, colour = quantile, shape = quantile)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(ymin = irf_down, ymax = irf_up), width = 0.18, linewidth = 0.8) +
  geom_point(size = 3.5) +
  facet_wrap(~ sample, nrow = 2) +
  scale_colour_manual(values = pal_since, guide = "none") +
  scale_shape_manual(values  = c(19, 17, 15), guide = "none") +
  labs(title    = "Leave-One-Continent-Out at h\u00a0=\u00a05",
       subtitle = "Point estimate \u00b1 95% CI by \u03c4 evaluation point",
       x = "\u03c4 (yrs since last TC)", y = "Cumulative GDP growth (%) at h=5") +
  theme_nature()

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

cat("\nBinary-group gradient at h=5:\n")
irf_clim2 %>% filter(horizon == 5) %>%
  select(sample, quantile, irf_mean, se, irf_down, irf_up) %>%
  mutate(across(where(is.numeric), ~round(. * 100, 2))) %>%
  as.data.frame() %>% print()

p_clim2 <- irf_ribbon(irf_clim2, facet_var = "sample") +
  labs(title    = "Within Climatology Groups (Binary)",
       subtitle = "Does the \u03c4 gradient survive splitting at W\u0304 = 20 m/s?")

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

cat("\nTercile-group gradient at h=5:\n")
irf_clim3 %>% filter(horizon == 5) %>%
  select(sample, quantile, irf_mean, se, irf_down, irf_up) %>%
  mutate(across(where(is.numeric), ~round(. * 100, 2))) %>%
  as.data.frame() %>% print()

if (nrow(irf_clim3) > 0) {
  p_clim3 <- irf_ribbon(irf_clim3, facet_var = "sample") +
    labs(title    = "Within Climatology Groups (Fixed Breaks: 10 and 30 m/s)",
         subtitle = "T1: 0\u201310 m/s | T2: 10\u201330 m/s | T3: >30 m/s")
  ggsave("figures/compound_rob_clim_terciles.png", p_clim3,
         width = 14, height = 6, dpi = 300)
  message("Saved compound_rob_clim_terciles.png")
}

# =============================================================================
# 7. Residualised τ: project years_since_any on W̄, use residual
# =============================================================================
# This is the most direct remedy for the W̄–τ confound.
# τ_resid = years_since_any − E[years_since_any | W̄]
# Removes the cross-sectional level correlation; residual variation is
# idiosyncratic timing within what the country's climatology would predict.
#
# Evaluation points on the τ_resid scale: 10th, 50th, 90th percentiles
# among non-NA τ_resid values, reflecting "short relative to W̄",
# "average", and "long relative to W̄".

cat("\n--- Residualised τ (τ | W̄ removed) ---\n")

tau_resid_qtiles <- quantile(data$tau_resid, c(0.10, 0.50, 0.90), na.rm = TRUE)
cat("τ_resid evaluation quantiles (p10/p50/p90):", round(tau_resid_qtiles, 2), "\n")

eval_vals_r   <- as.numeric(tau_resid_qtiles)
eval_labels_r <- c("\u03c4_r p10 (short)", "\u03c4_r p50 (avg)", "\u03c4_r p90 (long)")

pal_resid <- c("\u03c4_r p10 (short)" = myred,
               "\u03c4_r p50 (avg)"   = mygold,
               "\u03c4_r p90 (long)"  = myblue)

irf_tau_resid <- bind_rows(
  run_lp(filter(data, W_bar > 0),  "Ever-treated, \u03c4_resid",
         state_var = "tau_resid",
         eval_v = eval_vals_r, eval_l = eval_labels_r),
  run_lp(data,                     "Full sample, \u03c4_resid",
         state_var = "tau_resid",
         eval_v = eval_vals_r, eval_l = eval_labels_r)
) %>%
  mutate(quantile = factor(quantile, levels = eval_labels_r),
         sample   = factor(sample, levels = c("Full sample, \u03c4_resid",
                                              "Ever-treated, \u03c4_resid")))

cat("\nResidualised-τ gradient at h=5:\n")
irf_tau_resid %>% filter(horizon == 5) %>%
  select(sample, quantile, irf_mean, se, irf_down, irf_up) %>%
  mutate(across(where(is.numeric), ~round(. * 100, 2))) %>%
  as.data.frame() %>% print()

p_tau_resid <- irf_ribbon(irf_tau_resid, pal = pal_resid, facet_var = "sample") +
  labs(
    title    = "Compound LP with Residualised \u03c4",
    subtitle = paste0(
      "State = \u03c4 \u2212 E[\u03c4|W\u0304].  ",
      "Evaluated at p10/p50/p90 of residual (\u2248",
      paste(round(tau_resid_qtiles, 1), collapse = "/"), " yrs).  ",
      "Red = short relative to W\u0304, blue = long relative to W\u0304"
    )
  )

ggsave("figures/compound_rob_tau_resid.png", p_tau_resid,
       width = 11, height = 6, dpi = 300)
message("Saved compound_rob_tau_resid.png")

# =============================================================================
# 8. Gradient summary: ME(long τ) − ME(short τ) at h=5 across all cuts
# =============================================================================
# Compute the isolated-minus-compounding gap and its approximate SE at h=5.
# For original τ: long = "5 yrs", short = "1 yr"
# For τ_resid:    long = p90,     short = p10

gradient_h5 <- function(df, long_label, short_label, cut_name) {
  df %>%
    filter(horizon == 5, quantile %in% c(long_label, short_label)) %>%
    select(sample, quantile, irf_mean, se) %>%
    pivot_wider(names_from = quantile, values_from = c(irf_mean, se),
                names_repair = "unique") %>%
    mutate(
      cut      = cut_name,
      gradient = (.data[[paste0("irf_mean_", long_label)]] -
                  .data[[paste0("irf_mean_", short_label)]]) * 100,
      # Delta-method SE assuming independence across evaluation points
      grad_se  = sqrt((.data[[paste0("se_", long_label)]]^2 +
                       .data[[paste0("se_", short_label)]]^2)) * 100
    ) %>%
    select(cut, sample, gradient, grad_se)
}

grad_all <- bind_rows(
  gradient_h5(irf_samples,   "5 yrs", "1 yr",   "Sample restriction"),
  gradient_h5(irf_loco,      "5 yrs", "1 yr",   "LOCO"),
  gradient_h5(irf_clim2,     "5 yrs", "1 yr",   "Binary W\u0304 split"),
  gradient_h5(irf_clim3,     "5 yrs", "1 yr",   "W\u0304 split (10/30 m/s)"),
  gradient_h5(irf_tau_resid,
              eval_labels_r[3], eval_labels_r[1], "\u03c4 residualised on W\u0304")
) %>%
  mutate(
    grad_lo = gradient - 1.96 * grad_se,
    grad_hi = gradient + 1.96 * grad_se,
    cut = factor(cut, levels = c("Sample restriction", "LOCO",
                                 "Binary W\u0304 split", "W\u0304 split (10/30 m/s)",
                                 "\u03c4 residualised on W\u0304"))
  )

# Drop rows where the gradient is not identified (degenerate τ support:
# 100% of obs at τ ≤ 2, so evaluation at τ=5 is pure extrapolation).
# Flagged by implausibly large SE (> 10 pp).
grad_all_plot <- grad_all %>% filter(abs(grad_se) <= 10)

cat("\n--- Gradient summary (isolated − compounding, pp at h=5) ---\n")
cat("Full table (including degenerate extrapolations):\n")
grad_all %>% mutate(across(where(is.numeric), ~round(., 2))) %>%
  as.data.frame() %>% print()
cat("\nFiltered table (SE ≤ 10 pp; excludes degenerate τ groups):\n")
grad_all_plot %>% mutate(across(where(is.numeric), ~round(., 2))) %>%
  as.data.frame() %>% print()

p_grad <- ggplot(grad_all_plot, aes(y = reorder(sample, gradient),
                                x = gradient, colour = cut)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = grad_lo, xmax = grad_hi),
                 height = 0.3, linewidth = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ cut, ncol = 1, scales = "free_y") +
  scale_colour_manual(
    values = c("Sample restriction"           = myblue,
               "LOCO"                          = mygold,
               "Binary W\u0304 split"          = myred,
               "W\u0304 split (10/30 m/s)"     = mygreen,
               "\u03c4 residualised on W\u0304" = "grey20"),
    guide = "none"
  ) +
  labs(
    title    = "Isolated-minus-Compounding Gradient at h\u00a0=\u00a05",
    subtitle = "ME(\u03c4=long) \u2212 ME(\u03c4=short), pp; positive = isolated events more damaging",
    x = "GDP growth difference (pp)", y = NULL
  ) +
  theme_nature()

ggsave("figures/compound_rob_gradient.png", p_grad,
       width = 10, height = 12, dpi = 300)
message("Saved compound_rob_gradient.png")

# ── Save ---------------------------------------------------------------------

saveRDS(
  list(irf_samples   = irf_samples,
       irf_loco       = irf_loco,
       irf_clim2      = irf_clim2,
       irf_clim3      = irf_clim3,
       irf_tau_resid  = irf_tau_resid,
       grad_all        = grad_all),
  "results_compound_robustness.rds"
)
message("Saved results_compound_robustness.rds")
cat("\nDone.\n")
