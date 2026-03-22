# =============================================================================
# Direct Damages Investigation — Data Construction
# =============================================================================
# Builds a self-contained panel merging PWT + IBTrACS + EMDAT.
# Key outputs:
#   - frac_dest_k   : EMDAT physical damages / capital stock (rkna)
#   - damage_pct_gdp: EMDAT damages / real GDP  (comparable to growth rate)
#   - frac_dest_h   : cyclone deaths / population
#   - ln_tfp_pop_hc : log TFP (Cobb-Douglas residual)
#   - loggdp, gdp_diff, maxwind, maxwind_sqkm, etc.
# Saved to: panel_damages.rds
# =============================================================================

library(tidyverse)
library(haven)
library(readxl)
library(countrycode)

knot_to_ms <- 0.514444

# ── 1. Load raw data ----------------------------------------------------------

pwt    <- read_stata("data/pwt_clean.dta")
tcs    <- read_stata("data/ibtracs_clean.dta")
wdi    <- read_stata("data/wdi_clean_2018.dta")
capidx <- read_stata("data/pwt_capitalprices.dta")

# ── 2. EMDAT: aggregate TC damages to country-year ---------------------------

emdat_raw <- read_excel("data/public_emdat_custom_request.xlsx",
                        sheet = "EM-DAT Data")

emdat <- emdat_raw %>%
  rename(
    countrycode  = ISO,
    year         = `Start Year`,
    total_deaths = `Total Deaths`,
    damage_adj   = `Total Damage, Adjusted ('000 US$)`  # already CPI-adjusted
  ) %>%
  mutate(
    year         = as.integer(year),
    countrycode  = as.character(countrycode),
    total_deaths = as.numeric(total_deaths),
    damage_adj   = as.numeric(damage_adj)
  ) %>%
  filter(year >= 1970, year <= 2015) %>%
  group_by(countrycode, year) %>%
  summarise(
    total_deaths    = sum(total_deaths, na.rm = TRUE),
    damage_adj_000  = sum(damage_adj,   na.rm = TRUE),   # thousands of $2011
    .groups = "drop"
  ) %>%
  mutate(
    damage_mil = damage_adj_000 / 1000   # → millions of $2011 (matches PWT units)
  )

# ── 3. Merge -----------------------------------------------------------------

data <- pwt %>%
  mutate(year = as.integer(year), countrycode = as.character(countrycode)) %>%
  left_join(
    wdi %>% mutate(year = as.integer(year), countrycode = as.character(countrycode)),
    by = c("countrycode", "year")
  ) %>%
  left_join(
    tcs %>% mutate(year = as.integer(year), countrycode = as.character(countrycode)),
    by = c("countrycode", "year")
  ) %>%
  left_join(emdat, by = c("countrycode", "year")) %>%
  filter(year >= 1970, year <= 2015)

# ── 4. TC shock variables (fill zeros) --------------------------------------

tc_zero_vars <- c("max_a_wind_i_sqkm_nob", "sum_a_lands_i_sqkm_nob",
                  "max_ann_wind_i_nob", "sum_ann_energy_i_sqkm_nob")
data <- data %>%
  mutate(across(all_of(tc_zero_vars), ~ replace_na(., 0)),
         total_deaths = replace_na(total_deaths, 0),
         damage_mil   = replace_na(damage_mil,   0))

# ── 5. Wind variables --------------------------------------------------------

data <- data %>%
  mutate(
    maxwind       = max_ann_wind_i_nob * knot_to_ms,
    maxwind_sqkm  = max_a_wind_i_sqkm_nob * knot_to_ms,
    nlands_sqkm   = sum_a_lands_i_sqkm_nob
  )

# ── 6. GDP outcome and lags --------------------------------------------------

data <- data %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  mutate(
    real_gdp_usd_pc        = if_else(!is.na(real_gdp_usd_pc),
                                     real_gdp_usd_pc, real_gdp_usd / pop),
    loggdp                 = log(real_gdp_usd_pc),
    lag_loggdp             = lag(loggdp),
    gdp_diff               = loggdp - lag_loggdp,
    lag_ln_real_gdp_usd_pc = lag(loggdp),
    lag_real_gdp_usd       = lag(real_gdp_usd, 1)
  ) %>%
  ungroup()

# ── 7. TFP ------------------------------------------------------------------

data <- data %>%
  mutate(
    ln_real_gdp_usd = log(real_gdp_usd),
    ln_rkna         = log(rkna),
    ln_hc           = log(hc),
    ln_pop          = log(pop),
    ln_tfp_pop_hc   = ln_real_gdp_usd - 0.33 * ln_rkna - 0.67 * (ln_hc + ln_pop)
  )

# ── 8. Direct damage variables -----------------------------------------------
# frac_dest_k: physical capital share destroyed (EMDAT damages / capital stock)
# damage_pct_gdp: damages as % of real GDP  -- comparable unit to GDP growth rate
# ln_frac_dest_k: log version for elasticity regressions

data <- data %>%
  mutate(
    dkil_annual     = total_deaths / 1e6,          # deaths → millions (PWT pop units)
    frac_dest_k     = damage_mil / rkna,            # both in millions $2011
    frac_dest_h     = dkil_annual / pop,
    damage_pct_gdp  = damage_mil / real_gdp_usd,   # fraction of GDP (both in mil $2011)
    # Log versions (storms only — NA otherwise)
    ln_frac_dest_k  = if_else(damage_mil > 0 & rkna > 0,
                              log(frac_dest_k), NA_real_),
    ln_frac_dest_h  = if_else(total_deaths > 0 & pop > 0,
                              log(frac_dest_h), NA_real_),
    ln_damage_pct_gdp = if_else(damage_mil > 0 & real_gdp_usd > 0,
                                log(damage_pct_gdp), NA_real_),
    ln_damage_lagGDP  = if_else(damage_mil > 0 & lag_real_gdp_usd > 0,
                                log(damage_mil / lag_real_gdp_usd), NA_real_)
  )

# ── 9. Shock indicators (p95 binary + Saffir-Simpson categories) ------------

data <- data %>%
  mutate(
    # Global p95 binary shock
    maxwind_p95   = quantile(maxwind[year >= 1970 & year <= 2014], 0.95, na.rm = TRUE),
    maxwind_95    = as.integer(maxwind_p95 > 0 & maxwind >= maxwind_p95),
    # Saffir-Simpson categorical binary indicators (m/s thresholds, no-storm baseline)
    wind_cat1plus = as.integer(maxwind >= 33),    # Cat 1+  (~64 knots)
    wind_cat2plus = as.integer(maxwind >= 43),    # Cat 2+  (~83 knots)
    wind_cat3plus = as.integer(maxwind >= 50),    # Cat 3+  (~96 knots)
    # Continuous area-normalised shock (knots/km²) — same as main Bakkensen-Barrage regressor
    wind_cont     = maxwind_sqkm
  )

# ── 10. Country characteristics (for interaction analysis) -------------------
# Vulnerability = country-level mean of log(damage / GDP_{t-1}) over storm-years
# with positive EMDAT damage.  Only countries with at least one positive-damage
# storm-year receive a non-NA vulnerability score.

data <- data %>%
  group_by(countrycode) %>%
  mutate(
    # Mean log(damage/GDP_{t-1}) across storm-years with positive reported damage
    damage_vuln = mean(
      ln_damage_lagGDP[maxwind > 0 & !is.na(ln_damage_lagGDP)],
      na.rm = TRUE
    ),
  ) %>%
  ungroup() %>%
  mutate(
    # Binary split: above / below median of damage_vuln
    # (only countries with at least one positive-damage storm receive a value)
    damage_binary = factor(
      ntile(damage_vuln, 2),
      levels = 1:2,
      labels = c("Low vulnerability", "High vulnerability")
    ),
    # Period indicator
    period  = if_else(year <= 1990, "Pre-1990", "Post-1990"),
    year2   = year^2,
    # Continent / region
    continent = countrycode(countrycode, "iso3c", "continent"),
    region    = countrycode(countrycode, "iso3c", "region")
  )

# ── 11. Contemporaneous damage residual (for the continuous interaction) -----
# damage_resid = damage_pct_gdp purged of wind (i.e., unexplained damage conditional on wind)
# Captures: country-specific vulnerability holding wind constant

data <- data %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  mutate(
    lag_damage_pct_gdp = lag(damage_pct_gdp, 1),
    any_damage         = as.integer(damage_mil > 0)
  ) %>%
  ungroup()

# ── 12. ENSO / climate controls (temperature, precipitation, Niño-3.4) -------

enso_raw <- read.csv("../../raw_data/ENSO_Growth_Panel_1960-2019.csv",
                     stringsAsFactors = FALSE)

enso <- enso_raw %>%
  select(countrycode = iso, year, t, p, nino34) %>%
  mutate(year = as.integer(year))

data <- data %>%
  left_join(enso, by = c("countrycode", "year"))

message("ENSO merge: t non-NA = ", sum(!is.na(data$t)),
        " | p non-NA = ", sum(!is.na(data$p)),
        " | nino34 non-NA = ", sum(!is.na(data$nino34)))

# ── 13. Save -----------------------------------------------------------------

saveRDS(data, "panel_damages.rds")
message("Saved panel_damages.rds: ", nrow(data), " rows, ",
        n_distinct(data$countrycode), " countries")
message("Storm-years with positive damage: ",
        sum(data$damage_mil > 0, na.rm = TRUE))
message("p95 threshold (global): ", round(data$maxwind_p95[1], 1), " m/s")
