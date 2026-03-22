# =============================================================================
# Direct Damages Investigation — Contemporaneous Damage ~ Wind
# =============================================================================

library(tidyverse)
library(fixest)
library(patchwork)
setFixest_notes(FALSE)

data <- readRDS("panel_damages.rds")

myblue   <- "#1f78b4"
myred    <- "#e31a1c"
myorange <- "#f28e2b"

FE    <- "countrycode[year] + countrycode[year2] + region^year"
PANEL <- c("countrycode", "year")
DK    <- DK ~ year

theme_nature <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      colour = "grey10", margin = margin(b = 5)),
      axis.title       = element_text(colour = "grey15", size = base_size),
      axis.text        = element_text(colour = "grey25", size = base_size - 1),
      axis.line        = element_line(colour = "grey40", linewidth = 0.45),
      axis.ticks       = element_line(colour = "grey50", linewidth = 0.45),
      legend.position  = "top",
      legend.text      = element_text(colour = "grey20", size = base_size - 1),
      legend.key.size  = unit(1.3, "lines"),
      legend.title     = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(12, 16, 10, 12)
    )
}

# ── 1. Scatter: log(damage/GDP) vs maxwind — two panels ----------------------

storm_data <- data %>% filter(maxwind > 0, !is.na(ln_damage_lagGDP))

make_scatter <- function(df, xvar, xlabel, panel_title) {
  ct  <- cor.test(df[[xvar]], df$ln_damage_lagGDP, method = "pearson")
  r   <- round(ct$estimate, 3)
  pv  <- ct$p.value
  plab <- if (pv < 0.001) "p < 0.001" else paste0("p = ", round(pv, 3))
  ann <- paste0("r = ", r, ",  ", plab)

  ggplot(df, aes(x = .data[[xvar]], y = ln_damage_lagGDP)) +
    geom_point(alpha = 0.22, size = 1.1, colour = myblue) +
    geom_smooth(method = "lm", colour = myred, fill = myred,
                alpha = 0.10, linewidth = 1.0, se = TRUE,
                linetype = "dashed") +
    annotate("text",
             x = -Inf, y = Inf,
             hjust = -0.1, vjust = 1.5,
             label = ann,
             size = 4.2, colour = "grey20") +
    labs(
      title = panel_title,
      x     = xlabel,
      y     = "log(damage / GDP[t-1])"
    ) +
    theme_nature() +
    theme(plot.title = element_text(size = 13, face = "plain",
                                    colour = "grey30", hjust = 0.5))
}

p_left <- make_scatter(storm_data, "maxwind",
                       "Max wind speed (m/s)",
                       "Max wind speed (m/s)")

storm_data_sqkm <- data %>%
  filter(maxwind_sqkm > 0, !is.na(ln_damage_lagGDP))

p_right <- make_scatter(storm_data_sqkm, "maxwind_sqkm",
                        "Max wind speed per km² (knots/km²)",
                        "Max wind speed per km²")

p_scatter <- (p_left | p_right) +
  plot_annotation(
    title = "Cyclone Damage vs Wind Speed",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, colour = "grey10",
                                hjust = 0.5, margin = margin(b = 8))
    )
  )

ggsave("figures/damage_vs_wind_scatter.png", p_scatter,
       width = 13, height = 6.5, dpi = 300)
message("Saved damage_vs_wind_scatter.png")

# ── 2. Point-whisker: category-specific coefficients per outcome -------------
# Cat 1 only: 33–43 m/s  |  Cat 2 only: 43–50 m/s  |  Cat 3+: ≥50 m/s

data <- data %>%
  mutate(
    wind_cat1only = as.integer(maxwind >= 33 & maxwind < 43),
    wind_cat2only = as.integer(maxwind >= 43 & maxwind < 50)
    # wind_cat3plus already in data
  )

cat_specs <- tribble(
  ~shock,          ~cat_label,
  "wind_cat1only", "Cat 1\n(33–43 m/s)",
  "wind_cat2only", "Cat 2\n(43–50 m/s)",
  "wind_cat3plus", "Cat 3+\n(≥ 50 m/s)"
)

outcome_specs <- tribble(
  ~var,                ~out_label,
  "ln_damage_lagGDP", "log(damage/GDP[t-1])",
  "ln_frac_dest_h",    "log(deaths/pop)"
)

coef_df <- cross_join(cat_specs, outcome_specs) %>%
  pmap_dfr(function(shock, cat_label, var, out_label) {
    m <- tryCatch(
      feols(as.formula(paste0(var, " ~ ", shock, " | ", FE)),
            data     = data %>% filter(!is.na(.data[[var]])),
            panel.id = PANEL,
            vcov     = ~countrycode),
      error = function(e) { message("Skip ", shock, "/", var, ": ", e$message); NULL }
    )
    if (is.null(m)) return(NULL)
    tibble(
      cat_label = cat_label,
      outcome   = out_label,
      est       = coef(m)[shock],
      se        = se(m)[shock],
      lo95      = est - 1.96  * se,
      hi95      = est + 1.96  * se,
      lo90      = est - 1.645 * se,
      hi90      = est + 1.645 * se
    )
  }) %>%
  mutate(
    cat_label = factor(cat_label, levels = cat_specs$cat_label),
    outcome   = factor(outcome,   levels = outcome_specs$out_label)
  )

pal_cat2 <- c(myblue, myorange, myred)

p_coef <- ggplot(coef_df, aes(x = cat_label, y = est, colour = cat_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60",
             linewidth = 0.5) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 0.9) +
  geom_linerange(aes(ymin = lo90, ymax = hi90), linewidth = 2.2) +
  geom_point(size = 4, shape = 21, fill = "white", stroke = 1.5) +
  scale_colour_manual(values = setNames(pal_cat2, levels(coef_df$cat_label)),
                      guide = "none") +
  facet_wrap(~ outcome, scales = "free_y") +
  labs(
    title = "Direct Damage by Cyclone Category",
    x     = NULL,
    y     = "Coefficient (log-points)"
  ) +
  theme_nature() +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 13, colour = "grey15"),
    axis.text.x      = element_text(size = 11, colour = "grey20")
  )

ggsave("figures/direct_damage_coef.png", p_coef,
       width = 11, height = 6, dpi = 300)
message("Saved direct_damage_coef.png")

# ── 3. feols marginal-effect curves: both log outcomes ~ maxwind -------------

outcomes_cont <- tribble(
  ~var,                ~label,              ~colour,
  "ln_damage_lagGDP", "log(damage/GDP[t-1])",   myred,
  "ln_frac_dest_h",    "log(deaths/pop)",   myblue
)

wind_grid <- seq(0, max(data$maxwind, na.rm = TRUE), length.out = 200)

pred_df <- pmap_dfr(outcomes_cont, function(var, label, colour) {
  df_fit <- data %>% filter(!is.na(.data[[var]]), maxwind > 0)
  m <- tryCatch(
    feols(as.formula(paste0(var, " ~ maxwind | ", FE)),
          data     = df_fit,
          panel.id = PANEL,
          vcov     = ~countrycode),
    error = function(e) { message("Skip ", var, ": ", e$message); NULL }
  )
  if (is.null(m)) return(NULL)
  b       <- coef(m)["maxwind"]
  se_b    <- se(m)["maxwind"]
  mu_wind <- mean(df_fit$maxwind, na.rm = TRUE)
  mu_y    <- mean(df_fit[[var]],  na.rm = TRUE)
  tibble(
    outcome = label,
    colour  = colour,
    maxwind = wind_grid,
    fit     = mu_y + b * (wind_grid - mu_wind),
    lo95    = mu_y + (b - 1.96 * se_b) * (wind_grid - mu_wind),
    hi95    = mu_y + (b + 1.96 * se_b) * (wind_grid - mu_wind)
  )
})

pal_out <- setNames(unique(pred_df$colour), unique(pred_df$outcome))

p_marg <- ggplot(pred_df, aes(x = maxwind, colour = outcome, fill = outcome)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.12, colour = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  scale_colour_manual(values = pal_out, name = NULL) +
  scale_fill_manual(values   = pal_out, name = NULL) +
  labs(
    title = "Damage & Mortality vs Wind Speed",
    x     = "Max wind speed (m/s)",
    y     = "ln(normalized loss)"
  ) +
  theme_nature()

ggsave("figures/damage_deaths_vs_maxwind_feols.png", p_marg,
       width = 7.5, height = 7, dpi = 300)
message("Saved damage_deaths_vs_maxwind_feols.png")

saveRDS(list(coef_df = coef_df, pred_df = pred_df),
        "results_direct_damage.rds")
message("Saved results_direct_damage.rds")
