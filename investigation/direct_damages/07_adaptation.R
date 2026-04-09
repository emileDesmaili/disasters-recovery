# =============================================================================
# Direct Damages Investigation — Adaptation Heterogeneity
# =============================================================================
# Measure: within-country wind-damage elasticity.
#
# Regression: log(damage/GDP_{t-1}) ~ maxwind | country FE
# Estimated as a linear mixed-effects model where the slope of maxwind
# is a *random* effect (country-specific slope shrunk toward pooled mean).
# Partial pooling solves the small-sample OLS noise problem: countries with
# few damage observations get slopes pulled toward the global average, while
# countries with many observations keep their estimated slopes.
#
# Interpretation:
#   High slope (steep) = every extra m/s causes proportionally more damage
#                        = not adapted (infrastructure destroyed by intensity)
#   Low  slope (flat)  = damage roughly constant across wind speeds
#                        = adapted (infrastructure absorbs wind variation)
#
# Steps:
#  1. Mixed-effects model on damage-years: lmer(ln_damage_lagGDP ~ maxwind +
#     (maxwind | countrycode)); extract BLUP slopes per country
#  2. Describe terciles and show scatter of fitted lines
#  3. LP IRF heterogeneity by adaptation tercile
# =============================================================================

library(tidyverse)
library(fixest)
library(lme4)
library(patchwork)
library(ggrepel)
source("../../emileRegs.R")
setFixest_notes(FALSE)
setFixest_nthreads(1)

data <- readRDS("panel_damages.rds")

myblue   <- "#1f78b4"
myred    <- "#e31a1c"
myorange <- "#f28e2b"

FE      <- "countrycode[year] + countrycode[year2] + region^year"
PANEL   <- c("countrycode", "year")
DK      <- DK ~ year
HORIZON <- 8

theme_nature <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      colour = "grey10", margin = margin(b = 5)),
      plot.subtitle    = element_text(colour = "grey35", size = base_size - 1,
                                      margin = margin(b = 8)),
      axis.title       = element_text(colour = "grey15", size = base_size),
      axis.text        = element_text(colour = "grey25", size = base_size - 1),
      axis.line        = element_line(colour = "grey40", linewidth = 0.45),
      axis.ticks       = element_line(colour = "grey50", linewidth = 0.45),
      legend.position  = "top",
      legend.text      = element_text(colour = "grey20", size = base_size - 1),
      legend.key.size  = unit(1.2, "lines"),
      legend.title     = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(12, 16, 10, 12),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", colour = "grey15",
                                      size = base_size - 1)
    )
}

# Country name lookup
cnames <- data %>%
  distinct(countrycode) %>%
  mutate(country_name = countrycode::countrycode(countrycode, "iso3c",
                                                  "country.name", warn = FALSE)) %>%
  mutate(country_name = coalesce(country_name, countrycode))

# =============================================================================
# 1. Mixed-effects model: country-specific damage-wind slopes with shrinkage
# =============================================================================
cat("=== 1. Mixed-effects damage-wind model ===\n")

# Damage-years only: ln_damage_lagGDP is NA unless EMDAT damage > 0
damage_data <- data %>%
  filter(!is.na(ln_damage_lagGDP), maxwind > 0) %>%
  select(countrycode, year, ln_damage_lagGDP, maxwind,
         loggdp, gdp_diff, wind_cat1plus, region)

cat("  Damage-year obs:", nrow(damage_data), "\n")
cat("  Countries:", n_distinct(damage_data$countrycode), "\n")
obs_per <- count(damage_data, countrycode)
cat("  Obs distribution (n per country):\n")
print(table(cut(obs_per$n, c(0, 2, 5, 10, 20, 50))))

# Mixed-effects: random intercept + random slope for maxwind by country
# Random slope = the country-specific wind-damage elasticity
set.seed(42)
me_mod <- lmer(
  ln_damage_lagGDP ~ maxwind + (1 + maxwind | countrycode),
  data    = damage_data,
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl   = list(maxfun = 2e5))
)
cat("\nFixed effects (pooled slope):\n")
print(fixef(me_mod))

# BLUP (Best Linear Unbiased Predictors) = shrunk country-specific slopes
blup     <- ranef(me_mod)$countrycode
pooled_b <- fixef(me_mod)["maxwind"]

country_slopes <- blup %>%
  rownames_to_column("countrycode") %>%
  rename(rand_int = `(Intercept)`, rand_slope = maxwind) %>%
  mutate(
    # Total slope = pooled slope + country-specific deviation
    adapt_slope = pooled_b + rand_slope,
    n_damage    = obs_per$n[match(countrycode, obs_per$countrycode)]
  ) %>%
  left_join(cnames, by = "countrycode")

cat("\nAdaptation slopes (sorted, most adapted first):\n")
tmp <- country_slopes %>%
  select(country_name, adapt_slope, n_damage) %>%
  arrange(adapt_slope) %>%
  mutate(adapt_slope = round(adapt_slope, 4))
write.csv(tmp, row.names = FALSE)
rm(tmp)

# Tercile classification
country_slopes <- country_slopes %>%
  mutate(
    adapt_tercile = ntile(adapt_slope, 3),
    adapt_label   = factor(adapt_tercile,
                           levels = 1:3,
                           labels = c("High adaptation",
                                      "Mid adaptation",
                                      "Low adaptation"))
  )

pal_tercile <- c(
  "High adaptation" = myblue,
  "Mid adaptation"  = myorange,
  "Low adaptation"  = myred
)

cat("\nAdaptation tercile summary:\n")
print(country_slopes %>%
        group_by(adapt_label) %>%
        summarise(
          n             = n(),
          slope_min     = round(min(adapt_slope), 4),
          slope_median  = round(median(adapt_slope), 4),
          slope_max     = round(max(adapt_slope), 4),
          countries     = paste(sort(country_name), collapse = ", "),
          .groups = "drop"
        ))

# =============================================================================
# 2. Visualise: country-level damage-wind scatter with fitted lines per tercile
# =============================================================================
cat("\n=== 2. Figures ===\n")

# Fitted line endpoints per country (from BLUP intercept + slope)
fitted_lines <- country_slopes %>%
  mutate(
    intercept = fixef(me_mod)["(Intercept)"] + rand_int,
    wind_lo   = damage_data %>%
                  group_by(countrycode) %>%
                  summarise(lo = min(maxwind), .groups = "drop") %>%
                  pull(lo) %>% .[match(countrycode,
                    damage_data %>% distinct(countrycode) %>% pull(countrycode))],
    wind_hi   = damage_data %>%
                  group_by(countrycode) %>%
                  summarise(hi = max(maxwind), .groups = "drop") %>%
                  pull(hi) %>% .[match(countrycode,
                    damage_data %>% distinct(countrycode) %>% pull(countrycode))],
    y_lo      = intercept + adapt_slope * wind_lo,
    y_hi      = intercept + adapt_slope * wind_hi
  )

# ── Ranked slopes ─────────────────────────────────────────────────────────────
p_rank <- ggplot(country_slopes,
                 aes(x = adapt_slope,
                     y = reorder(country_name, -adapt_slope),
                     colour = adapt_label)) +
  geom_vline(xintercept = pooled_b, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_point(aes(size = n_damage), alpha = 0.85) +
  scale_colour_manual(values = pal_tercile) +
  scale_size_continuous(name = "Damage obs.", range = c(1.5, 5),
                        breaks = c(3, 10, 20, 40)) +
  annotate("text", x = pooled_b + 0.001, y = 1,
           label = "Pooled slope", colour = "grey30",
           hjust = 0, size = 3.5) +
  labs(
    title    = "Country-Specific Wind-Damage Elasticities",
    subtitle = paste0(
      "Mixed-effects model (random slopes shrunk toward pooled average).\n",
      "Low slope = damage rises little with wind intensity = adapted."
    ),
    x = "Wind-damage slope (log damage/GDP per m/s, shrunk BLUP)",
    y = NULL
  ) +
  theme_nature(base_size = 11) +
  theme(axis.text.y = element_text(size = 8.5)) +
  guides(colour = guide_legend(override.aes = list(size = 3.5)))

ggsave("figures/adapt_rank.png", p_rank,
       width = 9, height = max(7, nrow(country_slopes) * 0.28), dpi = 300,
       limitsize = FALSE)
message("Saved adapt_rank.png")

# ── Damage vs maxwind scatter, coloured by adaptation tercile ─────────────────
p_scatter <- ggplot(damage_data %>%
                      left_join(country_slopes %>%
                                  select(countrycode, adapt_label),
                                by = "countrycode"),
                    aes(x = maxwind, y = ln_damage_lagGDP,
                        colour = adapt_label)) +
  geom_point(alpha = 0.30, size = 1.2) +
  # Draw pooled line
  geom_abline(intercept = fixef(me_mod)["(Intercept)"],
              slope     = pooled_b,
              colour    = "grey20", linewidth = 0.9, linetype = "solid") +
  scale_colour_manual(values = pal_tercile, na.value = "grey80") +
  labs(
    title    = "log(Damage/GDP) vs Wind Intensity by Adaptation Tercile",
    subtitle = "Each point = one country-year with positive EMDAT damage.\nBlack line = pooled within-country slope.",
    x        = "Max wind speed in strike year (m/s)",
    y        = "log(damage / GDP)"
  ) +
  theme_nature() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 0.8)))

ggsave("figures/adapt_scatter.png", p_scatter,
       width = 9, height = 6.5, dpi = 300)
message("Saved adapt_scatter.png")

# =============================================================================
# 3. Merge adaptation tercile into full panel for LP
# =============================================================================
cat("\n=== 3. Merging adaptation group into full panel ===\n")

panel_adapt <- data %>%
  left_join(country_slopes %>%
              select(countrycode, adapt_tercile, adapt_label, adapt_slope),
            by = "countrycode")

cat("  Panel obs by tercile:\n")
print(count(panel_adapt, adapt_label))

# =============================================================================
# 4. LP IRF by adaptation tercile
# =============================================================================
cat("\n=== 4. LP IRF by adaptation tercile ===\n")

run_lp_irf <- function(dat, label) {
  cat("  Fitting LP for:", label, "\n")
  map_dfr(0:HORIZON, function(h) {
    tryCatch({
      mod <- feols(
        as.formula(paste0(
          "lead(loggdp, ", h, ") - lag(loggdp, 1) ~ ",
          "wind_cat1plus + l(wind_cat1plus, 1:2) + l(gdp_diff, 1:2) | ", FE
        )),
        data = dat, panel.id = PANEL, vcov = DK
      )
      v   <- vcov(summary(mod, vcov = DK))["wind_cat1plus", "wind_cat1plus"]
      est <- coef(mod)["wind_cat1plus"]
      tibble(horizon = h, group = label,
             irf_mean = est, irf_se = sqrt(v),
             irf_down = est - 1.645 * sqrt(v),
             irf_up   = est + 1.645 * sqrt(v))
    }, error = function(e) NULL)
  })
}

irf_adapt <- imap_dfr(
  setNames(1:3, c("High adaptation", "Mid adaptation", "Low adaptation")),
  function(tercile, lbl) {
    panel_adapt %>% filter(adapt_tercile == tercile) %>% run_lp_irf(lbl)
  }
) %>%
  mutate(group = factor(group, levels = c("High adaptation", "Mid adaptation",
                                          "Low adaptation")))

p_irf_adapt <- ggplot(irf_adapt,
                      aes(x = horizon, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = irf_down * 100, ymax = irf_up * 100),
              alpha = 0.10, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60",
             linewidth = 0.45) +
  geom_line(aes(y = irf_mean * 100), linewidth = 1.4) +
  scale_x_continuous(breaks = 0:HORIZON) +
  scale_colour_manual(values = pal_tercile) +
  scale_fill_manual(values   = pal_tercile) +
  labs(
    title    = "GDP Response to Cat 1+ Cyclone by Adaptation Level",
    subtitle = "Adaptation = country wind-damage elasticity (mixed-effects slopes).\n90% CI (DK SE). Blue = most adapted countries.",
    x        = "Years after strike",
    y        = "Cumulative GDP change (pp)"
  ) +
  theme_nature() +
  guides(colour = guide_legend(nrow = 1,
                               override.aes = list(linewidth = 2, fill = NA)))

ggsave("figures/irf_by_adaptation.png", p_irf_adapt,
       width = 8, height = 6.5, dpi = 300)
message("Saved irf_by_adaptation.png")

pal_hl <- c("High adaptation" = myblue, "Low adaptation" = myred)

p_irf_hl <- irf_adapt %>%
  filter(group != "Mid adaptation") %>% droplevels() %>%
  ggplot(aes(x = horizon, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = irf_down * 100, ymax = irf_up * 100),
              alpha = 0.12, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60",
             linewidth = 0.45) +
  geom_line(aes(y = irf_mean * 100), linewidth = 1.6) +
  scale_x_continuous(breaks = 0:HORIZON) +
  scale_colour_manual(values = pal_hl) +
  scale_fill_manual(values   = pal_hl) +
  labs(
    title    = "Adapted Countries Recover Faster",
    subtitle = "High adaptation (blue) vs low adaptation (red). 90% CI (DK SE).",
    x        = "Years after strike",
    y        = "Cumulative GDP change (pp)"
  ) +
  theme_nature() +
  guides(colour = guide_legend(nrow = 1,
                               override.aes = list(linewidth = 2, fill = NA)))

ggsave("figures/irf_by_adaptation_hl.png", p_irf_hl,
       width = 7.5, height = 6.5, dpi = 300)
message("Saved irf_by_adaptation_hl.png")

# =============================================================================
saveRDS(list(country_slopes = country_slopes,
             me_mod         = me_mod,
             irf_adapt      = irf_adapt),
        "results_adaptation.rds")
message("Saved results_adaptation.rds\nDone.")
