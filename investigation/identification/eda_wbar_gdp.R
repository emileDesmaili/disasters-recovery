# =============================================================================
# EDA: Country-average wind exposure vs log GDP per capita (1970)
# Two panels:
#   Left  — ever-struck countries only (W_bar > 0)
#   Right — all countries
# =============================================================================

library(tidyverse)
library(patchwork)
library(latex2exp)
library(ggrepel)

source("../../emileRegs.R")

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
      legend.position  = "none",
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(10, 14, 8, 10)
    )
}

# ── Data ----------------------------------------------------------------------

data <- readRDS("../direct_damages/panel_damages.rds") %>%
  arrange(countrycode, year) %>%
  group_by(countrycode) %>%
  mutate(W_bar = mean(maxwind, na.rm = TRUE)) %>%
  ungroup()

# Country-level: W_bar + loggdp in 1970
cty <- data %>%
  group_by(countrycode) %>%
  summarise(
    W_bar    = first(W_bar),
    loggdp70 = loggdp[year == 1970][1],
    .groups  = "drop"
  ) %>%
  filter(!is.na(loggdp70))

ever_struck <- cty %>% filter(W_bar > 0)

# ── Correlation labels --------------------------------------------------------

r_all    <- cor(cty$W_bar,       cty$loggdp70,       use = "complete.obs")
r_struck <- cor(ever_struck$W_bar, ever_struck$loggdp70, use = "complete.obs")

lab_all    <- sprintf("italic(r) == '%.2f'~~(italic(n) == %d)", r_all,    nrow(cty %>% filter(!is.na(W_bar))))
lab_struck <- sprintf("italic(r) == '%.2f'~~(italic(n) == %d)", r_struck, nrow(ever_struck))

myblue <- "#2c7bb6"
myred  <- "#d7191c"

# ── Panel A: ever-struck -------------------------------------------------------

label_countries <- c("USA", "PHL", "JPN", "DOM", "MEX", "CHN", "IDN")
label_names     <- c("USA" = "United States", "PHL" = "Philippines",
                     "JPN" = "Japan",         "DOM" = "Dom. Republic",
                     "MEX" = "Mexico",        "CHN" = "China",
                     "IDN" = "Indonesia")

ever_struck_labels <- ever_struck %>%
  mutate(label = ifelse(countrycode %in% label_countries,
                        label_names[countrycode], NA_character_))

pA <- ggplot(ever_struck_labels, aes(x = W_bar, y = loggdp70)) +
  geom_point(colour = myblue, alpha = 0.65, size = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = myred,
              fill = myred, alpha = 0.15, linewidth = 1) +
  geom_label_repel(
    data = ever_struck_labels %>% filter(!is.na(label)),
    aes(label = label),
    size          = 3.8,
    colour        = "grey10",
    fill          = "white",
    label.size    = 0.25,
    label.padding = unit(0.2, "lines"),
    segment.colour = "grey40",
    segment.size  = 0.4,
    arrow         = arrow(length = unit(0.012, "npc"), type = "closed"),
    box.padding   = 0.6,
    point.padding = 0.3,
    min.segment.length = 0,
    seed          = 42
  ) +
  annotate("text", x = Inf, y = Inf, label = lab_struck,
           hjust = 1.05, vjust = 1.8, size = 4.5, colour = "grey20",
           parse = TRUE) +
  labs(
    title    = "Ever-Struck Countries",
    subtitle = TeX("$\\bar{W} > 0$ m/s"),
    x        = TeX("Country-average max wind speed, $\\bar{W}$ (m/s)"),
    y        = TeX("Log GDP per capita, 1970")
  ) +
  theme_nature()

# ── Panel B: all countries ----------------------------------------------------

pB <- ggplot(cty %>% filter(!is.na(W_bar)), aes(x = W_bar, y = loggdp70)) +
  geom_point(colour = myblue, alpha = 0.55, size = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = myred,
              fill = myred, alpha = 0.15, linewidth = 1) +
  annotate("text", x = Inf, y = Inf, label = lab_all,
           hjust = 1.05, vjust = 1.8, size = 4.5, colour = "grey20",
           parse = TRUE) +
  labs(
    title    = "All Countries",
    subtitle = "Full sample",
    x        = TeX("Country-average max wind speed, $\\bar{W}$ (m/s)"),
    y        = TeX("Log GDP per capita, 1970")
  ) +
  theme_nature()

# ── Combine -------------------------------------------------------------------

p_combined <- pA + pB +
  plot_annotation(
    title    = "Cyclone Exposure Is Orthogonal to Pre-Sample Income",
    subtitle = TeX("Country-average wind exposure $\\bar{W}$ vs. log GDP per capita in 1970.  OLS fit $\\pm$ 95% CI."),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15,
                                   colour = "grey10", hjust = 0.5,
                                   margin = margin(b = 4)),
      plot.subtitle = element_text(colour = "grey45", size = 12,
                                   hjust = 0.5, margin = margin(b = 8)),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

ggsave("figures/eda_wbar_gdp1970.png", p_combined,
       width = 12, height = 5.5, dpi = 300)
message("Saved figures/eda_wbar_gdp1970.png")

cat("Correlation (ever-struck):", round(r_struck, 3), "\n")
cat("Correlation (all):",         round(r_all,    3), "\n")
