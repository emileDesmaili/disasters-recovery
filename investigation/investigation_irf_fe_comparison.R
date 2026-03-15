# investigation_irf_fe_comparison.R
# Compare period-interacted IRFs across three FE specifications:
#   1. Baseline: countrycode[year] + countrycode[year2] + year
#   2. Region x Year + Country Polynomials: countrycode[year] + countrycode[year2] + region^year
#   3. Region x Year Only: region^year

library(tidyverse)
library(fixest)
library(haven)
library(countrycode)

source("../emileRegs.R")

setFixest_notes(FALSE)

# ------------------------------------------------------------------
# 1. Load and construct data
# ------------------------------------------------------------------
pwt  <- read_stata("../raw_data/pwt_clean.dta")
tcs  <- read_stata("../raw_data/ibtracs_clean.dta")
data <- pwt %>% left_join(tcs, by = c("year", "countrycode"))

knot_to_ms <- 0.514444

data <- data %>%
  mutate(
    maxwind  = replace_na(max_ann_wind_i_nob, 0) * knot_to_ms,
    loggdp   = ln_real_gdp_usd_pc,
    gdp_diff = growth_real_gdp_usd_pc,
    year2    = year^2,
    period   = ifelse(year <= 1990, "Pre-1990", "Post-1990")
  ) %>%
  filter(year >= 1969)

data <- data %>%
  mutate(
    maxwind_p95 = quantile(maxwind[year >= 1970 & year <= 2014], 0.95, na.rm = TRUE)
  ) %>%
  mutate(
    maxwind_95 = as.integer(maxwind_p95 > 0 & maxwind >= maxwind_p95)
  )

data <- data %>%
  mutate(
    region = countrycode(countrycode, "iso3c", "region"),
    region = ifelse(is.na(region), "Other", region)
  )

# ------------------------------------------------------------------
# 2. Three FE specifications
# ------------------------------------------------------------------
fe_specs <- list(
  "Baseline (Year FE)"              = "countrycode[year] + countrycode[year2] + year",
  "Region x Year + Country Poly"    = "countrycode[year] + countrycode[year2] + region^year",
  "Region x Year Only"              = "region^year"
)

outcome  <- "loggdp"
horizon  <- 10
controls <- "l(gdp_diff,1:2) + l(maxwind_95,1:2)"
panel_id <- c("countrycode", "year")
vcov_fm  <- DK ~ year

# ------------------------------------------------------------------
# 3. Estimate IRFs
# ------------------------------------------------------------------
all_irfs <- list()

for (spec_name in names(fe_specs)) {
  fe <- fe_specs[[spec_name]]
  cat(sprintf("Estimating: %s\n", spec_name))

  irf <- lp_panel_inter(
    data = data, outcome = outcome, main_var = "maxwind_95",
    interact_var = "period", controls = controls,
    horizon = horizon, fe = fe, panel_id = panel_id, vcov_formula = vcov_fm
  )

  irf <- irf %>%
    mutate(
      category = str_replace_all(category, "_", "-"),
      category = case_when(
        str_detect(category, "Pre")  ~ "Pre-1990",
        str_detect(category, "Post") ~ "Post-1990",
        TRUE ~ category
      ),
      spec = spec_name
    )

  all_irfs[[spec_name]] <- irf
}

plot_df <- bind_rows(all_irfs) %>%
  mutate(spec = factor(spec, levels = names(fe_specs)))

# ------------------------------------------------------------------
# 4. Plot
# ------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

p <- ggplot(plot_df, aes(x = horizon, y = irf_mean, color = category, fill = category)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~spec, ncol = 3) +
  scale_color_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_fill_manual(values = c("Pre-1990" = "steelblue", "Post-1990" = "firebrick")) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = 0:10) +
  labs(
    x = "Horizon (years)",
    y = "Cumulative effect on log GDP per capita",
    title = "Period-Interacted IRF: Comparing FE Specifications",
    subtitle = "Maxwind p95 binary shock | 95% CI (Driscoll-Kraay)",
    color = NULL, fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "bottom"
  )

ggsave("figures/irf_fe_comparison.png", p,
       width = 14, height = 5.5, dpi = 300)
cat("\nSaved: figures/irf_fe_comparison.png\n")

# Print h=5 gaps
cat("\nPost-Pre gaps at h=5:\n")
for (spec_name in names(fe_specs)) {
  h5 <- plot_df %>% filter(spec == spec_name, horizon == 5)
  post <- h5 %>% filter(category == "Post-1990") %>% pull(irf_mean)
  pre  <- h5 %>% filter(category == "Pre-1990")  %>% pull(irf_mean)
  cat(sprintf("  %s: %.2f pp\n", spec_name, (post - pre) * 100))
}
