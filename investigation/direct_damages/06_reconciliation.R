# =============================================================================
# Direct Damages — Reconciling State-Dependent LP vs Vulnerability-Tercile LP
# =============================================================================
#
# The state-dependent LP (04_heterogeneity.R) and the vulnerability-tercile
# sub-sample LP produce different pictures of damage heterogeneity.  This
# script systematically tests four hypotheses for why they diverge:
#
#  H1. SAMPLE SELECTION
#      State-dep LP only uses obs where ln_damage_lagGDP is defined
#      (damage_mil > 0).  Tercile LP uses ALL country-years in each tercile.
#      → Test: run tercile LP restricted to storm-years (maxwind > 0).
#
#  H2. IDENTIFICATION (within vs. between)
#      State-dep exploits within-country variation in annual damage.
#      Tercile splits exploit permanent between-country differences.
#      → Test: run state-dep LP using country-mean ln_damage_lagGDP
#              (time-invariant) as the state variable.
#
#  H3. FUNCTIONAL FORM
#      Linear interaction may miss discrete jumps captured by binary splits.
#      → Test: re-evaluate state-dep LP at binary-boundary quantiles
#              (p25 / p75) rather than p10 / p50 / p90.
#
# EDA:
#   EDA-1.  Distribution of contemporaneous ln_damage_lagGDP by binary group
#           (storm-years with damage > 0 only).  Tests whether the permanent
#           index aligns with the contemporaneous state variable.
#   EDA-2.  Country composition (income, continent) within each group.
#
# Figures (6):
#   reconcile_damage_dist.png       EDA-1
#   reconcile_composition.png       EDA-2
#   reconcile_overlay.png           Direct comparison: both approaches
#   reconcile_stormyears.png        H1: binary LP, storm-years only
#   reconcile_timeinvariant.png     H2: state-dep with country-mean state var
#   reconcile_tercile_cuts.png      H3: state-dep evaluated at p25/p75
# =============================================================================

library(tidyverse)
library(fixest)
library(patchwork)
source("../../emileRegs.R")
setFixest_notes(FALSE)

data    <- readRDS("panel_damages.rds")
res_het <- readRDS("results_heterogeneity.rds")   # irf_state, irf_by_vuln

# Binary indicator for High vulnerability (used in interaction LPs)
data <- data %>%
  mutate(damage_high_bin = as.integer(damage_binary == "High vulnerability"))

# ── Palette & theme -----------------------------------------------------------

myblue   <- "#1f78b4"
myred    <- "#e31a1c"
myorange <- "#f28e2b"
mygreen  <- "#2ca25f"
mypurple <- "#756bb1"

pal_tercile <- c(
  "Low damage\nvulnerability"  = myblue,
  "Medium"                     = myorange,
  "High damage\nvulnerability" = myred
)

pal_binary <- c(
  "Low vulnerability"  = myblue,
  "High vulnerability" = myred
)

FE      <- "countrycode[year] + countrycode[year2] + region^year"
PANEL   <- c("countrycode", "year")
DK      <- DK ~ year
HORIZON <- 8

theme_nature <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 1,
                                   colour = "grey10", margin = margin(b = 4)),
      plot.subtitle = element_text(colour = "grey45", size = base_size - 2,
                                   margin = margin(b = 8)),
      axis.title    = element_text(colour = "grey15", size = base_size),
      axis.text     = element_text(colour = "grey25", size = base_size - 1),
      axis.line     = element_line(colour = "grey40", linewidth = 0.45),
      axis.ticks    = element_line(colour = "grey50", linewidth = 0.45),
      legend.position  = "top",
      legend.text      = element_text(colour = "grey20", size = base_size - 1),
      legend.key.size  = unit(1.3, "lines"),
      legend.title     = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(12, 16, 10, 12)
    )
}

irf_line <- function(df, pal, lt = NULL, sub = NULL) {
  p <- ggplot(df, aes(x = horizon, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = irf_down * 100, ymax = irf_up * 100),
                alpha = 0.10, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.5) +
    geom_line(aes(y = irf_mean * 100), linewidth = 1.4) +
    scale_x_continuous(breaks = 0:HORIZON) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_fill_manual(values   = pal, name = NULL) +
    labs(x = "Years after strike", y = "Cumulative GDP growth (%)") +
    theme_nature() +
    guides(colour = guide_legend(nrow = 1))
  if (!is.null(sub)) p <- p + labs(subtitle = sub)
  if (!is.null(lt))  p <- p + scale_linetype_manual(values = lt, name = NULL)
  p
}

# ─────────────────────────────────────────────────────────────────────────────
# DIAGNOSTIC TABLES
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=========================================================\n")
cat("DIAGNOSTIC: state variable vs vulnerability index alignment\n")
cat("=========================================================\n")

# Global quantiles of the contemporaneous state variable
storm_obs <- data %>% filter(!is.na(ln_damage_lagGDP))
q_global  <- quantile(storm_obs$ln_damage_lagGDP,
                      c(0.10, 0.25, 0.50, 0.75, 0.90),
                      na.rm = TRUE)
cat("\nGlobal quantiles of ln_damage_lagGDP (storm-years, damage > 0):\n")
print(round(q_global, 3))

# Per-tercile summary of ln_damage_lagGDP (contemporaneous)
tab_dist <- data %>%
  filter(!is.na(ln_damage_lagGDP), !is.na(damage_binary)) %>%
  group_by(damage_binary) %>%
  summarise(
    n_obs       = n(),
    n_countries = n_distinct(countrycode),
    mean        = mean(ln_damage_lagGDP),
    sd          = sd(ln_damage_lagGDP),
    p10         = quantile(ln_damage_lagGDP, 0.10),
    p50         = quantile(ln_damage_lagGDP, 0.50),
    p90         = quantile(ln_damage_lagGDP, 0.90),
    .groups     = "drop"
  )
cat("\nln_damage_lagGDP summary by vulnerability tercile (storm-years):\n")
print(tab_dist, digits = 3)

# Damage-positive obs per binary group (state-dep LP sample vs binary LP sample)
tab_sample <- data %>%
  filter(!is.na(damage_binary)) %>%
  group_by(damage_binary) %>%
  summarise(
    total_obs     = n(),
    damage_pos    = sum(!is.na(ln_damage_lagGDP)),
    storm_years   = sum(maxwind > 0, na.rm = TRUE),
    pct_dam_pos   = round(100 * damage_pos / total_obs, 1),
    .groups = "drop"
  )
cat("\nSample sizes by binary group — total vs damage-positive vs storm-year obs:\n")
print(tab_sample)

# Within-country SD of ln_damage_lagGDP by binary group
# (how much within-country variation does the state-dep LP actually use?)
tab_within <- data %>%
  filter(!is.na(ln_damage_lagGDP), !is.na(damage_binary)) %>%
  group_by(countrycode, damage_binary) %>%
  summarise(within_sd = sd(ln_damage_lagGDP), n = n(), .groups = "drop") %>%
  group_by(damage_binary) %>%
  summarise(
    n_countries    = n(),
    med_within_sd  = median(within_sd, na.rm = TRUE),
    mean_within_sd = mean(within_sd,   na.rm = TRUE),
    .groups = "drop"
  )
cat("\nWithin-country SD of ln_damage_lagGDP by binary group:\n")
print(tab_within, digits = 3)

# ─────────────────────────────────────────────────────────────────────────────
# EDA-1: Distribution of ln_damage_lagGDP by vulnerability binary group
# ─────────────────────────────────────────────────────────────────────────────

eda1 <- data %>%
  filter(!is.na(ln_damage_lagGDP), !is.na(damage_binary))

p_eda1 <- ggplot(eda1, aes(x = damage_binary, y = ln_damage_lagGDP,
                            fill = damage_binary, colour = damage_binary)) +
  geom_violin(alpha = 0.22, linewidth = 0.55, trim = TRUE) +
  geom_boxplot(width = 0.16, outlier.size = 0.7, outlier.alpha = 0.35,
               fill = "white", linewidth = 0.5) +
  # p25/p75 reference lines matching state-dep LP evaluation points
  geom_hline(yintercept = q_global["25%"], linetype = "dotted",
             colour = myblue, linewidth = 0.9) +
  geom_hline(yintercept = q_global["75%"], linetype = "dotted",
             colour = myred, linewidth = 0.9) +
  annotate("text", x = 0.55, y = q_global["25%"] + 0.2,
           label = "p25", colour = myblue, size = 5.5, hjust = 0) +
  annotate("text", x = 0.55, y = q_global["75%"] + 0.2,
           label = "p75", colour = myred, size = 5.5, hjust = 0) +
  scale_fill_manual(values = pal_binary, guide = "none") +
  scale_colour_manual(values = pal_binary, guide = "none") +
  labs(
    title    = "Contemporaneous Damage Distribution by Vulnerability Group",
    subtitle = "Dotted lines = global p25 / p75 evaluation points used in state-dep LP",
    x        = NULL,
    y        = "log(damage / GDP[t-1])"
  ) +
  theme_nature()

ggsave("figures/reconcile_damage_dist.png", p_eda1,
       width = 8.5, height = 6.5, dpi = 300)
message("Saved reconcile_damage_dist.png")

# ─────────────────────────────────────────────────────────────────────────────
# EDA-2: Regional composition by vulnerability group
# ─────────────────────────────────────────────────────────────────────────────

comp_cont <- data %>%
  distinct(countrycode, damage_binary, continent) %>%
  filter(!is.na(damage_binary), !is.na(continent)) %>%
  count(damage_binary, continent) %>%
  group_by(damage_binary) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

p_comp <- ggplot(comp_cont,
                 aes(x = damage_binary, y = pct, fill = continent)) +
  geom_col(position = "stack", alpha = 0.85, width = 0.65) +
  scale_y_continuous(labels = scales::percent_format(1)) +
  labs(
    title    = "Regional Composition by Vulnerability Group (Binary Split)",
    subtitle = "Systematic concentration in certain regions signals between-country identification",
    x = NULL, y = "Share of countries"
  ) +
  theme_nature() +
  guides(fill = guide_legend(nrow = 2))

p_comp <- p_comp +
  theme(plot.title    = element_text(face = "bold", size = 15, colour = "grey10"),
        plot.subtitle = element_text(colour = "grey45", size = 12,
                                     margin = margin(b = 6)))

ggsave("figures/reconcile_composition.png", p_comp,
       width = 13, height = 6.5, dpi = 300)
message("Saved reconcile_composition.png")

# ─────────────────────────────────────────────────────────────────────────────
# EDA-3: World map of vulnerability tercile assignment
# ─────────────────────────────────────────────────────────────────────────────

library(maps)

vuln_map <- data %>%
  distinct(countrycode, damage_binary) %>%
  mutate(vuln_lab = case_when(
    damage_binary == "Low vulnerability"  ~ "Low",
    damage_binary == "High vulnerability" ~ "High",
    TRUE                                  ~ "No data"
  )) %>%
  mutate(vuln_lab = factor(vuln_lab, levels = c("Low", "High", "No data")))

world <- map_data("world") %>%
  mutate(iso3 = countrycode::countrycode(region, "country.name", "iso3c",
                                         warn = FALSE))

map_vuln_df <- world %>%
  left_join(vuln_map, by = c("iso3" = "countrycode")) %>%
  mutate(vuln_lab = replace_na(as.character(vuln_lab), "No data"),
         vuln_lab = factor(vuln_lab, levels = c("Low", "High", "No data")))

pal_map <- c("Low"     = "#6baed6",
             "High"    = "#cb181d",
             "No data" = "grey88")

p_map_vuln <- ggplot(map_vuln_df,
                     aes(x = long, y = lat, group = group,
                         fill = vuln_lab)) +
  geom_polygon(colour = "white", linewidth = 0.06) +
  scale_fill_manual(values = pal_map, name = "Vulnerability") +
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-60, 85)) +
  labs(title = "Vulnerability Binary Split by Country",
       subtitle = "Above/below median mean log(damage / GDP[t-1]) in positive-damage storm-years") +
  theme_void(base_size = 16) +
  theme(
    plot.title    = element_text(face = "bold", colour = "grey10",
                                 hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(colour = "grey40", hjust = 0.5,
                                 size = 13, margin = margin(b = 6)),
    legend.position = "bottom",
    legend.text     = element_text(size = 13),
    legend.title    = element_text(size = 13),
    plot.background = element_rect(fill = "white", colour = NA)
  )

ggsave("figures/reconcile_vuln_map.png", p_map_vuln,
       width = 11, height = 6, dpi = 300)
message("Saved reconcile_vuln_map.png")

# ─────────────────────────────────────────────────────────────────────────────
# EDA-4: Balance check — is binary assignment "random"?
# ─────────────────────────────────────────────────────────────────────────────
# Compare pre-determined country characteristics across binary groups:
#   - Log GDP per capita (1970-1990 avg)
#   - Log population
#   - Log land area (sq km)
#   - Storm frequency (fraction of years with maxwind > 0)
#   - Conditional wind intensity (mean maxwind among strike-years)
#   - Capital per worker (log rkna / pop)
# If the binary split is quasi-randomly assigned, distributions should overlap.

cat("\n=== EDA-4: Balance check across vulnerability groups (binary) ===\n")

# Land area (km²) from Natural Earth — try sf, fall back to terra
land_area_df <- tryCatch({
  suppressMessages({
    w <- rnaturalearth::ne_countries(returnclass = "sf")
    sf::sf_use_s2(FALSE)
    w$area_sqkm <- as.numeric(sf::st_area(w)) / 1e6
    as.data.frame(w) %>%
      transmute(countrycode = iso_a3,
                log_land_area = log(area_sqkm + 1)) %>%
      filter(!is.na(countrycode), countrycode != "-99")
  })
}, error = function(e) {
  message("sf area failed (", conditionMessage(e), "); trying terra...")
  tryCatch({
    suppressMessages({
      w <- rnaturalearth::ne_countries(returnclass = "sv")
      areas <- terra::expanse(w, unit = "km")
      data.frame(
        countrycode   = w$iso_a3,
        log_land_area = log(areas + 1)
      ) %>% filter(!is.na(countrycode), countrycode != "-99")
    })
  }, error = function(e2) {
    message("terra area also failed; log_land_area will be NA")
    data.frame(countrycode = character(0), log_land_area = numeric(0))
  })
})

balance_vars <- data %>%
  filter(!is.na(damage_binary)) %>%
  group_by(countrycode, damage_binary) %>%
  summarise(
    log_gdppc_pre   = mean(log(real_gdp_usd_pc)[year <= 1990],  na.rm = TRUE),
    log_pop         = mean(log(pop),                             na.rm = TRUE),
    storm_freq      = mean(maxwind > 0,                          na.rm = TRUE),
    cond_wind       = mean(maxwind[maxwind > 0],                 na.rm = TRUE),
    log_cap_worker  = mean(log(rkna / pop)[year <= 1990],        na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(land_area_df, by = "countrycode")

# Print means by binary group
cat("\nMeans by binary vulnerability group:\n")
balance_vars %>%
  group_by(damage_binary) %>%
  summarise(across(where(is.numeric), ~round(mean(., na.rm = TRUE), 3)),
            n = n(), .groups = "drop") %>%
  print()

# Density plots for each covariate
bal_long <- balance_vars %>%
  select(countrycode, damage_binary, log_gdppc_pre, log_pop, log_land_area,
         storm_freq, cond_wind, log_cap_worker) %>%
  pivot_longer(cols = -c(countrycode, damage_binary),
               names_to = "variable", values_to = "value") %>%
  mutate(variable = recode(variable,
    log_gdppc_pre  = "Log GDP/capita (pre-1990)",
    log_pop        = "Log population",
    log_land_area  = "Log land area (km\u00b2)",
    storm_freq     = "Storm frequency",
    cond_wind      = "Cond. wind (m/s)",
    log_cap_worker = "Log capital/worker (pre-1990)"
  ))

pal_bal <- c(
  "Low vulnerability"  = myblue,
  "High vulnerability" = myred
)

p_balance <- ggplot(bal_long %>% filter(!is.na(value)),
                    aes(x = value, colour = damage_binary,
                        fill = damage_binary)) +
  geom_density(alpha = 0.15, linewidth = 0.9) +
  geom_rug(aes(colour = damage_binary), alpha = 0.4, linewidth = 0.4,
           sides = "b") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  scale_colour_manual(values = pal_bal, name = NULL) +
  scale_fill_manual(values   = pal_bal, name = NULL) +
  labs(
    title    = "Balance Check: Country Characteristics by Vulnerability Group (Binary)",
    subtitle = "Non-overlap indicates binary group correlates with structural country features",
    x = NULL, y = "Density"
  ) +
  theme_nature(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 13, colour = "grey15")
  ) +
  guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))

ggsave("figures/reconcile_balance.png", p_balance,
       width = 13, height = 10, dpi = 300)
message("Saved reconcile_balance.png")

# ─────────────────────────────────────────────────────────────────────────────
# OVERLAY: direct visual comparison of both approaches
# ─────────────────────────────────────────────────────────────────────────────

state_overlay <- res_het$irf_state %>%
  mutate(
    group  = case_when(
      as.character(quantile) %in% c("Low damage (p10)", "p10") ~ "Low damage (p25)",
      as.character(quantile) %in% c("High damage (p90)", "p90") ~ "High damage (p75)",
      TRUE ~ as.character(quantile)
    ),
    ltype  = "State-dep (contemporaneous)",
    source = "State-dep"
  )

vuln_overlay <- res_het$irf_by_vuln %>%
  mutate(
    ltype  = "Binary sub-sample",
    source = "Binary"
  ) %>%
  rename(group = group)

pal_ov <- c(
  "Low damage (p25)"   = "#6baed6",
  "High damage (p75)"  = "#cb181d",
  "Low vulnerability"  = "#08519c",
  "High vulnerability" = "#67000d"
)
lt_ov <- c(
  "Low damage (p25)"   = "dashed",
  "High damage (p75)"  = "dashed",
  "Low vulnerability"  = "solid",
  "High vulnerability" = "solid"
)

overlay_df <- bind_rows(state_overlay, vuln_overlay)

p_overlay <- ggplot(overlay_df, aes(x = horizon, colour = group,
                                    linetype = group)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  geom_line(aes(y = irf_mean * 100), linewidth = 1.4) +
  scale_x_continuous(breaks = 0:HORIZON) +
  scale_colour_manual(values = pal_ov, name = NULL) +
  scale_linetype_manual(values = lt_ov, name = NULL) +
  labs(
    title    = "State-Dep LP vs Binary Sub-Sample LP",
    subtitle = "Dashed = state-dep at p25/p75 | Solid = binary sub-sample",
    x = "Years after strike", y = "Cumulative GDP growth (%)"
  ) +
  theme_nature() +
  guides(colour   = guide_legend(nrow = 2),
         linetype = guide_legend(nrow = 2))

ggsave("figures/reconcile_overlay.png", p_overlay,
       width = 9, height = 7, dpi = 300)
message("Saved reconcile_overlay.png")

# ─────────────────────────────────────────────────────────────────────────────
# H1: Full-panel interaction LP — discrete vs continuous
# ─────────────────────────────────────────────────────────────────────────────
# Direct comparison: lp_panel_inter (discrete damage_binary) on full data —
# same as 04_heterogeneity.R section 2.  This is the baseline for comparing
# H2 (time-invariant state) and H3 (continuous at p25/p75).

cat("\n=== H1: Interaction LP — full panel (discrete binary) ===\n")

irf_h1 <- lp_panel_inter(
  data         = data,
  outcome      = "loggdp",
  main_var     = "wind_cat1plus",
  interact_var = "damage_binary",
  controls     = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2)",
  horizon      = HORIZON,
  fe           = FE,
  panel_id     = PANEL,
  vcov_formula = DK
) %>% rename(group = category)

if (!is.null(irf_h1) && nrow(irf_h1) > 0) {
  cat("\nH1 results at h=5 (interaction LP, full panel):\n")
  irf_h1 %>% filter(horizon == 5) %>%
    mutate(across(c(irf_mean, se, irf_down, irf_up), ~round(. * 100, 2))) %>%
    print()

  p_h1 <- irf_line(irf_h1 %>% mutate(group = factor(group, levels = names(pal_binary))),
                   pal_binary,
                   sub = "H1: Discrete interaction LP on full panel (baseline for H2/H3)") +
    labs(title = "Interaction LP by Binary Vulnerability (Full Panel)")

  ggsave("figures/reconcile_stormyears.png", p_h1,
         width = 7.5, height = 7, dpi = 300)
  message("Saved reconcile_stormyears.png")
}

# ─────────────────────────────────────────────────────────────────────────────
# H2: State-dep LP with time-invariant (country-mean) state variable
# ─────────────────────────────────────────────────────────────────────────────
# Replace contemporaneous ln_damage_lagGDP with each country's long-run mean
# over storm-years.  A time-invariant state variable is conceptually closer
# to the vulnerability index used in the tercile analysis.

data <- data %>%
  group_by(countrycode) %>%
  mutate(
    ln_damage_mean = mean(ln_damage_lagGDP[maxwind > 0], na.rm = TRUE)
  ) %>%
  ungroup()

cat("\n=== H2: State-dep LP — time-invariant country mean ===\n")
cat("Country-mean ln_damage quantiles:\n")
print(round(quantile(
  distinct(data, countrycode, ln_damage_mean)$ln_damage_mean,
  c(0.10, 0.33, 0.50, 0.67, 0.90), na.rm = TRUE
), 3))

irf_h2 <- lp_state_dep(
  data           = data,
  outcome        = "loggdp",
  main_var       = "wind_cat1plus",
  state_var      = "ln_damage_mean",
  controls       = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2)",
  horizon        = HORIZON,
  fe             = FE,
  panel_id       = PANEL,
  vcov_formula   = DK,
  eval_quantiles = c(0.10, 0.50, 0.90),
  eval_labels    = c("Low (p10)", "Median (p50)", "High (p90)")
) %>% rename(group = quantile)

if (!is.null(irf_h2) && nrow(irf_h2) > 0) {
  cat("\nH2 results at h=5:\n")
  irf_h2 %>% filter(horizon == 5) %>%
    mutate(across(c(irf_mean, se, irf_down, irf_up), ~round(. * 100, 2))) %>%
    print()

  pal_h2 <- c("Low (p10)" = myblue, "Median (p50)" = myorange,
              "High (p90)" = myred)

  p_h2 <- irf_line(irf_h2, pal_h2,
                   sub = "H2: Country-mean log(damage/GDP[t-1]) — between-country, time-invariant") +
    labs(title = "State-Dep LP: Time-Invariant State Variable")

  ggsave("figures/reconcile_timeinvariant.png", p_h2,
         width = 7.5, height = 7, dpi = 300)
  message("Saved reconcile_timeinvariant.png")
}

# ─────────────────────────────────────────────────────────────────────────────
# H3: State-dep LP evaluated at tercile boundary quantiles (p33 / p67)
# ─────────────────────────────────────────────────────────────────────────────
# Re-runs the original state-dep LP (contemporaneous state var) but evaluates
# the ME at the tercile cut points.  If the two approaches agree on the ranking
# but disagree on the level, they share the same functional form.

cat("\n=== H3: State-dep LP at binary cut points (p25 / p75) ===\n")

irf_h3 <- lp_state_dep(
  data           = data,
  outcome        = "loggdp",
  main_var       = "wind_cat1plus",
  state_var      = "ln_damage_lagGDP",
  controls       = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2)",
  horizon        = HORIZON,
  fe             = FE,
  panel_id       = PANEL,
  vcov_formula   = DK,
  eval_quantiles = c(0.25, 0.75),
  eval_labels    = c("p25 (low/med cut)", "p75 (med/high cut)")
) %>% rename(group = quantile)

if (!is.null(irf_h3) && nrow(irf_h3) > 0) {
  cat("\nH3 results at h=5:\n")
  irf_h3 %>% filter(horizon == 5) %>%
    mutate(across(c(irf_mean, se, irf_down, irf_up), ~round(. * 100, 2))) %>%
    print()

  pal_h3 <- c("p25 (low/med cut)"  = myblue,
              "p75 (med/high cut)" = myred)

  p_h3 <- irf_line(irf_h3, pal_h3,
                   sub = "H3: Same model as state-dep LP; evaluated at binary cut points (p25/p75)") +
    labs(title = "State-Dep LP Evaluated at Binary Boundaries (p25/p75)")

  ggsave("figures/reconcile_tercile_cuts.png", p_h3,
         width = 7.5, height = 7, dpi = 300)
  message("Saved reconcile_tercile_cuts.png")
}

# ─────────────────────────────────────────────────────────────────────────────
# ROBUSTNESS: Ever-treated countries (at least one year with maxwind > 0)
# ─────────────────────────────────────────────────────────────────────────────
# Repeat the two main heterogeneity approaches from 04_heterogeneity.R
# restricted to countries that were ever struck — excludes never-struck
# countries that contribute to FE precision but never enter as treated.

cat("\n=== ROBUSTNESS: Ever-treated countries ===\n")

ever_treated_countries <- data %>%
  filter(maxwind > 0) %>%
  pull(countrycode) %>%
  unique() %>%
  sort()

data_rob <- data %>% filter(countrycode %in% ever_treated_countries)
cat("Ever-treated countries:", length(ever_treated_countries),
    "| full panel:", n_distinct(data$countrycode), "\n")

# (a) Continuous interaction (state-dep LP) on ever-treated sample
irf_rob_state <- lp_state_dep(
  data           = data_rob,
  outcome        = "loggdp",
  main_var       = "wind_cat1plus",
  state_var      = "ln_damage_lagGDP",
  controls       = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2)",
  horizon        = HORIZON,
  fe             = FE,
  panel_id       = PANEL,
  vcov_formula   = DK,
  eval_quantiles = c(0.25, 0.75),
  eval_labels    = c("Low damage (p25)", "High damage (p75)")
)

# (b) Discrete interaction (lp_panel_inter) on ever-treated sample
irf_rob_vuln <- lp_panel_inter(
  data         = data_rob,
  outcome      = "loggdp",
  main_var     = "wind_cat1plus",
  interact_var = "damage_binary",
  controls     = "l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2)",
  horizon      = HORIZON,
  fe           = FE,
  panel_id     = PANEL,
  vcov_formula = DK
) %>% rename(group = category)

pal_state <- c("Low damage (p25)" = myblue, "High damage (p75)" = myred)

if (nrow(irf_rob_state) > 0 && nrow(irf_rob_vuln) > 0) {
  cat("\nRobustness state-dep at h=5:\n")
  irf_rob_state %>% filter(horizon == 5) %>%
    mutate(across(c(irf_mean, se, irf_down, irf_up), ~round(. * 100, 2))) %>%
    print()

  p_rob_state <- ggplot(
    irf_rob_state %>% mutate(quantile = factor(quantile, levels = names(pal_state))),
    aes(x = horizon, colour = quantile, fill = quantile)) +
    geom_ribbon(aes(ymin = irf_down * 100, ymax = irf_up * 100),
                alpha = 0.10, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60",
               linewidth = 0.5) +
    geom_line(aes(y = irf_mean * 100), linewidth = 1.4) +
    scale_x_continuous(breaks = 0:HORIZON) +
    scale_colour_manual(values = pal_state, name = NULL) +
    scale_fill_manual(values   = pal_state, name = NULL) +
    labs(
      title    = "State-Dep LP (Ever-Treated Countries)",
      subtitle = "Robustness: excludes never-struck countries",
      x = "Years after strike", y = "Cumulative GDP growth (%)"
    ) +
    theme_nature() + guides(colour = guide_legend(nrow = 1))

  p_rob_vuln <- irf_line(
    irf_rob_vuln %>% mutate(group = factor(group, levels = names(pal_binary))),
    pal_binary,
    sub = "Robustness: excludes never-struck countries"
  ) + labs(title = "Interaction LP by Vulnerability (Ever-Treated Countries)")

  p_rob <- (p_rob_state | p_rob_vuln) +
    plot_annotation(
      title = "Robustness: Ever-Treated Countries Only",
      theme = theme(plot.title = element_text(face = "bold", size = 17,
                                              colour = "grey10", hjust = 0.5))
    )

  ggsave("figures/reconcile_ever_treated.png", p_rob,
         width = 13, height = 7, dpi = 300)
  message("Saved reconcile_ever_treated.png")
}

# ─────────────────────────────────────────────────────────────────────────────
# Summary table at h = 5
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=========================================================\n")
cat("SUMMARY: All approaches at h = 5 (GDP growth, pp)\n")
cat("=========================================================\n")

summarise_at_h5 <- function(df, label) {
  df %>%
    filter(horizon == 5) %>%
    transmute(
      approach  = label,
      group     = group,
      est_pct   = round(irf_mean * 100, 2),
      se_pct    = round(se * 100, 2),
      ci90_lo   = round((irf_mean - 1.645 * se) * 100, 2),
      ci90_hi   = round((irf_mean + 1.645 * se) * 100, 2)
    )
}

bind_rows(
  summarise_at_h5(
    res_het$irf_state %>% rename(group = quantile), "State-dep LP (continuous)"),
  summarise_at_h5(res_het$irf_by_vuln, "Discrete interaction LP"),
  if (!is.null(irf_h1) && nrow(irf_h1) > 0)
    summarise_at_h5(irf_h1, "H1: discrete interaction (full)") else NULL,
  if (!is.null(irf_h2) && nrow(irf_h2) > 0)
    summarise_at_h5(irf_h2, "H2: time-invariant state") else NULL,
  if (!is.null(irf_h3) && nrow(irf_h3) > 0)
    summarise_at_h5(irf_h3, "H3: at binary cuts") else NULL
) %>%
  arrange(approach, group) %>%
  as.data.frame() %>%
  print()

# ─────────────────────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────────────────────

saveRDS(list(
  irf_h1          = irf_h1,
  irf_h2          = irf_h2,
  irf_h3          = irf_h3,
  irf_rob_state   = irf_rob_state,
  irf_rob_vuln    = irf_rob_vuln,
  overlay_df      = overlay_df
), "results_reconciliation.rds")
message("Saved results_reconciliation.rds")
cat("\nDone.\n")
