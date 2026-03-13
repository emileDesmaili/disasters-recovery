# investigation_globalization.R
# Tests whether globalization/trade integration explains the post-1990
# divergence in cyclone-GDP IRFs.
#
# Downloads additional PWT variables (trade openness, investment share),
# merges onto the panel, and runs:
#   1. State-dependent LP with trade openness interaction
#   2. State-dependent LP with investment share interaction
#   3. Period-interacted LP within trade-openness terciles

library(tidyverse)
library(fixest)
library(haven)
source("../emileRegs.R")

# ==================================================================
# 1. Load and prepare base data (same pipeline as lp_gdp_shocks.R)
# ==================================================================
pwt <- read_stata("../raw_data/pwt_clean.dta")
tcs <- read_stata("../raw_data/ibtracs_clean.dta")
data <- pwt %>% left_join(tcs, by = c("year", "countrycode"))

knot_to_ms <- 0.514444

data <- data %>%
  mutate(
    maxwind = replace_na(max_ann_wind_i_nob, 0) * knot_to_ms,
    energy  = replace_na(sum_ann_energy_i_nob, 0),
    nlands  = replace_na(sum_a_lands_nob, 0),
    loggdp   = ln_real_gdp_usd_pc,
    gdp_diff = growth_real_gdp_usd_pc,
    year2    = year^2,
    period   = ifelse(year <= 1990, "Pre-1990", "Post-1990")
  ) %>%
  filter(year >= 1969)

# Global 95th-percentile threshold
data <- data %>%
  mutate(
    maxwind_p95 = quantile(maxwind[year >= 1970 & year <= 2014], 0.95, na.rm = TRUE)
  ) %>%
  mutate(
    maxwind_95 = as.integer(maxwind_p95 > 0 & maxwind >= maxwind_p95)
  )

# Standardize continuous shock measure
data <- data %>%
  mutate(
    maxwind_sd = maxwind / sd(maxwind[maxwind > 0], na.rm = TRUE)
  )

# LP settings
outcome  <- "loggdp"
horizon  <- 10
make_controls <- function(shock_var) {
  paste0("l(gdp_diff,1:2) + l(", shock_var, ",1:2)")
}
fe       <- "countrycode[year] + countrycode[year2] + year"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

# Plotting theme
irf_theme <- theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# ==================================================================
# 2. Download full PWT data with trade/investment variables
# ==================================================================
cat("\n============================================================\n")
cat("Step 1: Loading PWT trade and investment variables\n")
cat("============================================================\n\n")

pwt_full_path <- "../raw_data/pwt_full.rds"
if (file.exists(pwt_full_path)) {
  cat("Loading cached PWT full data from pwt_full.rds...\n")
  pwt_full <- readRDS(pwt_full_path)
} else {
  # Try the pwt10 package approach
  if (!requireNamespace("pwt10", quietly = TRUE)) {
    cat("Installing pwt10 package...\n")
    install.packages("pwt10", repos = "https://cloud.r-project.org")
  }

  if (requireNamespace("pwt10", quietly = TRUE)) {
    cat("Loading PWT 10.01 from pwt10 package...\n")
    pwt_full <- pwt10::pwt10.01
    saveRDS(pwt_full, pwt_full_path)
    cat(sprintf("Saved PWT full data to %s (%d rows, %d columns)\n",
                pwt_full_path, nrow(pwt_full), ncol(pwt_full)))
  } else {
    # Fallback: try downloading CSV from PWT website
    cat("pwt10 package unavailable, trying CSV download...\n")
    pwt_url <- "https://dataverse.nl/api/access/datafile/354098"
    temp_file <- tempfile(fileext = ".xlsx")
    tryCatch({
      download.file(pwt_url, temp_file, mode = "wb", quiet = TRUE)
      if (!requireNamespace("readxl", quietly = TRUE)) {
        install.packages("readxl", repos = "https://cloud.r-project.org")
      }
      pwt_full <- readxl::read_excel(temp_file, sheet = "Data")
      saveRDS(pwt_full, pwt_full_path)
      cat("Downloaded and saved PWT data from Dataverse.\n")
    }, error = function(e) {
      stop("Could not load PWT data. Please install the pwt10 package:\n",
           "  install.packages('pwt10')\n",
           "Error: ", e$message)
    })
  }
}

# Extract the variables we need
# PWT uses 'isocode' as the country identifier (ISO 3-letter)
cat("\nExtracting trade and investment variables from PWT...\n")

# Identify the country code column
if ("isocode" %in% names(pwt_full)) {
  cc_col <- "isocode"
} else if ("countrycode" %in% names(pwt_full)) {
  cc_col <- "countrycode"
} else {
  # Try to find it
  cc_col <- grep("iso|country.*code", names(pwt_full), value = TRUE, ignore.case = TRUE)[1]
  cat(sprintf("  Using '%s' as country code column\n", cc_col))
}

pwt_trade <- pwt_full %>%
  as.data.frame() %>%
  select(
    countrycode = all_of(cc_col),
    year,
    csh_x,    # share of merchandise exports at current PPPs
    csh_m,    # share of merchandise imports at current PPPs
    csh_i,    # share of gross capital formation at current PPPs
    any_of(c("ctfp", "labsh"))  # TFP, labor share (if available)
  ) %>%
  mutate(
    year = as.integer(year),
    countrycode = as.character(countrycode),
    # Trade openness: exports + abs(imports)
    # csh_m is typically negative in PWT (imports are negative in expenditure decomposition)
    trade_open = csh_x + abs(csh_m),
    inv_share  = csh_i
  )

cat(sprintf("  PWT trade data: %d rows, %d countries, years %d-%d\n",
            nrow(pwt_trade), n_distinct(pwt_trade$countrycode),
            min(pwt_trade$year), max(pwt_trade$year)))
cat(sprintf("  Non-missing trade_open: %d (%.1f%%)\n",
            sum(!is.na(pwt_trade$trade_open)),
            100 * mean(!is.na(pwt_trade$trade_open))))
cat(sprintf("  Non-missing inv_share: %d (%.1f%%)\n",
            sum(!is.na(pwt_trade$inv_share)),
            100 * mean(!is.na(pwt_trade$inv_share))))

# ==================================================================
# 3. Merge onto main data
# ==================================================================
cat("\n============================================================\n")
cat("Step 2: Merging trade/investment variables onto panel\n")
cat("============================================================\n\n")

data <- data %>%
  left_join(
    pwt_trade %>% select(countrycode, year, trade_open, inv_share, csh_x, csh_m, csh_i),
    by = c("countrycode", "year")
  )

cat(sprintf("After merge: %d rows\n", nrow(data)))
cat(sprintf("  Non-missing trade_open: %d (%.1f%%)\n",
            sum(!is.na(data$trade_open)),
            100 * mean(!is.na(data$trade_open))))
cat(sprintf("  Non-missing inv_share: %d (%.1f%%)\n",
            sum(!is.na(data$inv_share)),
            100 * mean(!is.na(data$inv_share))))

# Demean state variables (for interaction models)
data <- data %>%
  mutate(
    trade_open_dm = trade_open - mean(trade_open, na.rm = TRUE),
    inv_share_dm  = inv_share  - mean(inv_share, na.rm = TRUE)
  )

# ==================================================================
# Descriptive statistics: trade openness and investment by period
# ==================================================================
cat("\n--- Trade Openness and Investment Share by Period ---\n")
desc_trade <- data %>%
  filter(!is.na(trade_open) | !is.na(inv_share)) %>%
  group_by(period) %>%
  summarise(
    n = n(),
    trade_open_mean   = mean(trade_open, na.rm = TRUE),
    trade_open_median = median(trade_open, na.rm = TRUE),
    trade_open_sd     = sd(trade_open, na.rm = TRUE),
    inv_share_mean    = mean(inv_share, na.rm = TRUE),
    inv_share_median  = median(inv_share, na.rm = TRUE),
    inv_share_sd      = sd(inv_share, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(desc_trade), digits = 4, row.names = FALSE)

cat("\n--- Trade Openness by Period for Cyclone-Exposed Countries ---\n")
desc_trade_exposed <- data %>%
  filter(maxwind > 0, !is.na(trade_open)) %>%
  group_by(period) %>%
  summarise(
    n_exposed = n(),
    trade_open_mean   = mean(trade_open, na.rm = TRUE),
    trade_open_median = median(trade_open, na.rm = TRUE),
    inv_share_mean    = mean(inv_share, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(desc_trade_exposed), digits = 4, row.names = FALSE)

# ==================================================================
# 4. State-dependent LP function
# ==================================================================

lp_state_dep <- function(data, outcome, main_var, state_var,
                         controls = NULL, horizon = 10,
                         fe = "countrycode + year",
                         panel_id = c("countrycode", "year"),
                         vcov_formula = DK ~ year) {
  rhs_controls <- if (!is.null(controls)) paste0(" + ", controls) else ""
  irf_list <- list()

  for (h in 0:horizon) {
    fml <- as.formula(paste0(
      "f(", outcome, ", ", h, ") - l(", outcome, ", 1) ~ ",
      main_var, " + ", main_var, ":", state_var,
      rhs_controls, " | ", fe
    ))
    mod <- feols(fml, data = data, panel.id = panel_id, vcov = vcov_formula)

    # Evaluate marginal effect at different state values
    beta_main  <- coef(mod)[main_var]
    inter_name <- grep(paste0(main_var, ":", state_var), names(coef(mod)), value = TRUE)
    if (length(inter_name) == 0)
      inter_name <- grep(paste0(state_var, ":", main_var), names(coef(mod)), value = TRUE)
    beta_inter <- coef(mod)[inter_name]

    # Evaluate at 10th, 50th, 90th percentile of state_var
    state_vals <- quantile(data[[state_var]], c(0.10, 0.50, 0.90), na.rm = TRUE)
    for (j in seq_along(state_vals)) {
      me <- beta_main + beta_inter * state_vals[j]
      # Delta method SE
      grad <- c(1, state_vals[j])
      se <- sqrt(t(grad) %*% vcov(mod)[c(main_var, inter_name),
                                        c(main_var, inter_name)] %*% grad)
      irf_list[[length(irf_list) + 1]] <- data.frame(
        horizon = h, irf_mean = me, se = as.numeric(se),
        irf_down = me - 1.96 * as.numeric(se),
        irf_up = me + 1.96 * as.numeric(se),
        quantile = paste0("p", c(10, 50, 90)[j]),
        state_val = state_vals[j]
      )
    }
  }
  bind_rows(irf_list)
}

# ==================================================================
# ANALYSIS 1: Trade Openness Interaction
# ==================================================================
cat("\n\n============================================================\n")
cat("ANALYSIS 1: State-Dependent LP -- Trade Openness Interaction\n")
cat("============================================================\n\n")

cat("Running state-dependent LP: maxwind_sd x trade_open_dm...\n")

irfs_trade <- lp_state_dep(
  data = data,
  outcome = outcome,
  main_var = "maxwind_sd",
  state_var = "trade_open_dm",
  controls = make_controls("maxwind_sd"),
  horizon = horizon,
  fe = fe,
  panel_id = panel_id,
  vcov_formula = vcov_fm
) %>%
  mutate(
    quantile = factor(quantile, levels = c("p10", "p50", "p90"),
                      labels = c("Low openness (p10)",
                                 "Median openness (p50)",
                                 "High openness (p90)"))
  )

# Print summary
cat("\n--- IRF at h=5 and h=10, by trade openness quantile ---\n")
irfs_trade %>%
  filter(horizon %in% c(5, 10)) %>%
  select(horizon, quantile, state_val, irf_mean, se, irf_down, irf_up) %>%
  arrange(horizon, quantile) %>%
  as.data.frame() %>%
  print(digits = 4, row.names = FALSE)

# Plot
fig_trade <- ggplot(irfs_trade, aes(x = horizon, y = irf_mean,
                                     color = quantile, fill = quantile)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Low openness (p10)" = "steelblue",
                                "Median openness (p50)" = "grey40",
                                "High openness (p90)" = "firebrick")) +
  scale_fill_manual(values = c("Low openness (p10)" = "steelblue",
                               "Median openness (p50)" = "grey40",
                               "High openness (p90)" = "firebrick")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "Horizon (years)",
    y = "GDP response (per 1 SD wind)",
    title = "Cyclone IRF by Trade Openness",
    subtitle = "Marginal effect at p10, p50, p90 of trade openness (demeaned)",
    color = "Trade openness", fill = "Trade openness"
  ) +
  irf_theme +
  theme(legend.position = "top")

ggsave("figures/irf_trade_openness.png", fig_trade,
       width = 8, height = 5.5, dpi = 300)
cat("  Saved: figures/irf_trade_openness.png\n")

# ==================================================================
# ANALYSIS 2: Investment Share Interaction
# ==================================================================
cat("\n\n============================================================\n")
cat("ANALYSIS 2: State-Dependent LP -- Investment Share Interaction\n")
cat("============================================================\n\n")

cat("Running state-dependent LP: maxwind_sd x inv_share_dm...\n")

irfs_inv <- lp_state_dep(
  data = data,
  outcome = outcome,
  main_var = "maxwind_sd",
  state_var = "inv_share_dm",
  controls = make_controls("maxwind_sd"),
  horizon = horizon,
  fe = fe,
  panel_id = panel_id,
  vcov_formula = vcov_fm
) %>%
  mutate(
    quantile = factor(quantile, levels = c("p10", "p50", "p90"),
                      labels = c("Low investment (p10)",
                                 "Median investment (p50)",
                                 "High investment (p90)"))
  )

# Print summary
cat("\n--- IRF at h=5 and h=10, by investment share quantile ---\n")
irfs_inv %>%
  filter(horizon %in% c(5, 10)) %>%
  select(horizon, quantile, state_val, irf_mean, se, irf_down, irf_up) %>%
  arrange(horizon, quantile) %>%
  as.data.frame() %>%
  print(digits = 4, row.names = FALSE)

# Plot
fig_inv <- ggplot(irfs_inv, aes(x = horizon, y = irf_mean,
                                 color = quantile, fill = quantile)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Low investment (p10)" = "steelblue",
                                "Median investment (p50)" = "grey40",
                                "High investment (p90)" = "firebrick")) +
  scale_fill_manual(values = c("Low investment (p10)" = "steelblue",
                               "Median investment (p50)" = "grey40",
                               "High investment (p90)" = "firebrick")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "Horizon (years)",
    y = "GDP response (per 1 SD wind)",
    title = "Cyclone IRF by Investment Share",
    subtitle = "Marginal effect at p10, p50, p90 of investment share (demeaned)",
    color = "Investment share", fill = "Investment share"
  ) +
  irf_theme +
  theme(legend.position = "top")

ggsave("figures/irf_investment_share.png", fig_inv,
       width = 8, height = 5.5, dpi = 300)
cat("  Saved: figures/irf_investment_share.png\n")

# ==================================================================
# ANALYSIS 3: Trade Openness Terciles x Period
# ==================================================================
cat("\n\n============================================================\n")
cat("ANALYSIS 3: Period-Interacted LP by Trade Openness Terciles\n")
cat("============================================================\n\n")

# Compute country-level mean trade openness
country_trade <- data %>%
  filter(!is.na(trade_open)) %>%
  group_by(countrycode) %>%
  summarise(
    mean_trade_open = mean(trade_open, na.rm = TRUE),
    mean_inv_share  = mean(inv_share, na.rm = TRUE),
    n_years         = n(),
    n_shocks_95     = sum(maxwind_95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    trade_tercile = ntile(mean_trade_open, 3),
    trade_tercile_label = factor(
      trade_tercile,
      levels = 1:3,
      labels = c("T1 (Least Open)", "T2 (Middle)", "T3 (Most Open)")
    )
  )

data <- data %>%
  left_join(
    country_trade %>% select(countrycode, mean_trade_open, trade_tercile,
                              trade_tercile_label),
    by = "countrycode"
  )

# Print tercile composition
cat("--- Trade Openness Tercile Composition ---\n")
tercile_summary <- country_trade %>%
  group_by(trade_tercile_label) %>%
  summarise(
    n_countries    = n(),
    min_trade_open = min(mean_trade_open),
    max_trade_open = max(mean_trade_open),
    mean_trade_open_grp = mean(mean_trade_open),
    n_with_shocks  = sum(n_shocks_95 > 0),
    total_shocks   = sum(n_shocks_95),
    .groups = "drop"
  )
print(as.data.frame(tercile_summary), digits = 4, row.names = FALSE)

# Run period-interacted LP within each tercile
cat("\nRunning period-interacted LP within each trade-openness tercile...\n")

irfs_tercile <- map_dfr(1:3, function(t) {
  label <- levels(country_trade$trade_tercile_label)[t]
  cat(sprintf("  Estimating LP for %s...\n", label))
  sub_data <- data %>% filter(trade_tercile == t)

  tryCatch({
    lp_panel_inter(
      data = sub_data,
      outcome = outcome,
      main_var = "maxwind_95",
      interact_var = "period",
      controls = make_controls("maxwind_95"),
      horizon = horizon,
      fe = fe,
      panel_id = panel_id,
      vcov_formula = vcov_fm
    ) %>%
      mutate(trade_group = label)
  }, error = function(e) {
    cat(sprintf("    Error for %s: %s\n", label, e$message))
    return(NULL)
  })
})

# Clean category names
irfs_tercile <- irfs_tercile %>%
  mutate(category = ifelse(grepl("Pre", category), "Pre-1990", "Post-1990"))

# Print summary at h=5
cat("\n--- IRF at h=5 by trade-openness tercile and period ---\n")
irfs_tercile %>%
  filter(horizon == 5) %>%
  select(trade_group, category, irf_mean, se, irf_down, irf_up) %>%
  arrange(trade_group, category) %>%
  as.data.frame() %>%
  print(digits = 4, row.names = FALSE)

cat("\n--- IRF at h=10 by trade-openness tercile and period ---\n")
irfs_tercile %>%
  filter(horizon == 10) %>%
  select(trade_group, category, irf_mean, se, irf_down, irf_up) %>%
  arrange(trade_group, category) %>%
  as.data.frame() %>%
  print(digits = 4, row.names = FALSE)

# Plot: 3-panel figure
fig_tercile <- irfs_tercile %>%
  mutate(
    category = factor(category, levels = c("Pre-1990", "Post-1990")),
    trade_group = factor(trade_group,
                         levels = c("T1 (Least Open)", "T2 (Middle)", "T3 (Most Open)"))
  ) %>%
  ggplot(aes(x = horizon, y = irf_mean, color = category, fill = category)) +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~trade_group, ncol = 3) +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_fill_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(
    x = "Horizon (years)", y = "GDP response",
    title = "Cyclone IRF by Period, Split by Trade Openness Tercile",
    subtitle = "maxwind_95 shock, period-interacted LP",
    color = "Period", fill = "Period"
  ) +
  irf_theme +
  theme(legend.position = "top")

ggsave("figures/irf_trade_terciles_by_period.png", fig_tercile,
       width = 13, height = 5, dpi = 300)
cat("  Saved: figures/irf_trade_terciles_by_period.png\n")

# ==================================================================
# SUMMARY
# ==================================================================
cat("\n\n============================================================\n")
cat("SUMMARY OF FINDINGS\n")
cat("============================================================\n\n")

cat("1. TRADE OPENNESS:\n")
cat(sprintf("   Global mean trade openness (exports + |imports|)/GDP: %.3f\n",
            mean(data$trade_open, na.rm = TRUE)))
cat(sprintf("   Pre-1990 mean: %.3f,  Post-1990 mean: %.3f\n",
            desc_trade$trade_open_mean[desc_trade$period == "Pre-1990"],
            desc_trade$trade_open_mean[desc_trade$period == "Post-1990"]))

cat("\n   State-dependent IRFs (h=10):\n")
for (q in levels(irfs_trade$quantile)) {
  val <- irfs_trade %>% filter(horizon == 10, quantile == q)
  cat(sprintf("     %s: %.4f (SE=%.4f, [%.4f, %.4f])\n",
              q, val$irf_mean, val$se, val$irf_down, val$irf_up))
}

cat("\n2. INVESTMENT SHARE:\n")
cat(sprintf("   Global mean investment share: %.3f\n",
            mean(data$inv_share, na.rm = TRUE)))
cat(sprintf("   Pre-1990 mean: %.3f,  Post-1990 mean: %.3f\n",
            desc_trade$inv_share_mean[desc_trade$period == "Pre-1990"],
            desc_trade$inv_share_mean[desc_trade$period == "Post-1990"]))

cat("\n   State-dependent IRFs (h=10):\n")
for (q in levels(irfs_inv$quantile)) {
  val <- irfs_inv %>% filter(horizon == 10, quantile == q)
  cat(sprintf("     %s: %.4f (SE=%.4f, [%.4f, %.4f])\n",
              q, val$irf_mean, val$se, val$irf_down, val$irf_up))
}

cat("\n3. TRADE OPENNESS TERCILES x PERIOD:\n")
for (tg in levels(factor(irfs_tercile$trade_group,
                         levels = c("T1 (Least Open)", "T2 (Middle)", "T3 (Most Open)")))) {
  cat(sprintf("   %s:\n", tg))
  for (per in c("Pre-1990", "Post-1990")) {
    val <- irfs_tercile %>% filter(trade_group == tg, category == per, horizon == 10)
    if (nrow(val) > 0) {
      cat(sprintf("     %s h=10: %.4f (SE=%.4f)\n", per, val$irf_mean, val$se))
    }
  }
}

cat("\n\nKEY QUESTION: Does trade openness or investment share explain\n")
cat("why cyclone damage got worse post-1990?\n\n")

# Evaluate the evidence
trade_h10_low  <- irfs_trade %>% filter(horizon == 10, quantile == "Low openness (p10)")
trade_h10_high <- irfs_trade %>% filter(horizon == 10, quantile == "High openness (p90)")
inv_h10_low    <- irfs_inv %>% filter(horizon == 10, quantile == "Low investment (p10)")
inv_h10_high   <- irfs_inv %>% filter(horizon == 10, quantile == "High investment (p90)")

cat("EVIDENCE ASSESSMENT:\n")
cat(sprintf("  Trade openness: High-open h=10 IRF = %.4f vs Low-open = %.4f\n",
            trade_h10_high$irf_mean, trade_h10_low$irf_mean))
cat(sprintf("    => Difference = %.4f (more negative = more damage for open economies)\n",
            trade_h10_high$irf_mean - trade_h10_low$irf_mean))

cat(sprintf("  Investment share: High-inv h=10 IRF = %.4f vs Low-inv = %.4f\n",
            inv_h10_high$irf_mean, inv_h10_low$irf_mean))
cat(sprintf("    => Difference = %.4f\n",
            inv_h10_high$irf_mean - inv_h10_low$irf_mean))

if (trade_h10_high$irf_mean < trade_h10_low$irf_mean) {
  cat("\n  ==> More trade-open countries suffer LARGER GDP losses from cyclones.\n")
  cat("      Since trade openness increased post-1990, this could partially explain\n")
  cat("      the worsening of cyclone damage in the later period.\n")
} else {
  cat("\n  ==> More trade-open countries do NOT suffer larger GDP losses from cyclones.\n")
  cat("      Trade openness does NOT appear to explain the post-1990 divergence.\n")
}

if (inv_h10_high$irf_mean < inv_h10_low$irf_mean) {
  cat("\n  ==> Higher investment-share countries suffer LARGER GDP losses from cyclones.\n")
  cat("      If investment shares changed post-1990, this could contribute.\n")
} else {
  cat("\n  ==> Higher investment-share countries do NOT suffer larger GDP losses.\n")
  cat("      Investment share does NOT appear to explain the post-1990 divergence.\n")
}

cat("\n\nFigures saved:\n")
cat("  - figures/irf_trade_openness.png\n")
cat("  - figures/irf_investment_share.png\n")
cat("  - figures/irf_trade_terciles_by_period.png\n")

cat("\nDone.\n")
