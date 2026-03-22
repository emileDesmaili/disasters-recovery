# =============================================================================
# Direct Damages Investigation — Map of log(damage/GDP) by Country
# =============================================================================

library(tidyverse)
library(maps)

data <- readRDS("panel_damages.rds")

theme_nature_map <- function(base_size = 14) {
  theme_void(base_size = base_size) +
    theme(
      plot.title      = element_text(face = "bold", size = base_size + 1,
                                     colour = "grey10", margin = margin(b = 6),
                                     hjust = 0.5),
      legend.position = "bottom",
      legend.title    = element_text(colour = "grey20", size = base_size - 1,
                                     hjust = 0.5),
      legend.text     = element_text(colour = "grey30", size = base_size - 2),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin     = margin(10, 10, 10, 10)
    )
}

# ── Country-level: avg log(damage/GDP) in storm-years with positive damage --

country_damage <- data %>%
  filter(maxwind > 0, !is.na(ln_damage_lagGDP)) %>%
  group_by(countrycode) %>%
  summarise(
    avg_ln_damage = mean(ln_damage_lagGDP, na.rm = TRUE),
    n_storms      = n(),
    .groups = "drop"
  )

# ── Map data -----------------------------------------------------------------

world <- map_data("world") %>%
  mutate(iso3 = countrycode::countrycode(region, "country.name", "iso3c",
                                         warn = FALSE))

map_df <- world %>%
  left_join(country_damage, by = c("iso3" = "countrycode"))

p_map <- ggplot(map_df, aes(x = long, y = lat, group = group,
                             fill = avg_ln_damage)) +
  geom_polygon(colour = "white", linewidth = 0.08) +
  scale_fill_distiller(
    palette  = "Spectral",
    direction = -1,
    na.value = "grey88",
    name     = "Mean log(damage/GDP[t-1])\nin storm-years",
    guide    = guide_colorbar(barwidth = 12, barheight = 0.5,
                              title.position = "top", title.hjust = 0.5)
  ) +
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-60, 85)) +
  labs(title = "Mean log(Damage/GDP) by Country (EMDAT, 1970-2015)") +
  theme_nature_map()

ggsave("figures/map_direct_damage.png", p_map,
       width = 11, height = 6, dpi = 300)
message("Saved map_direct_damage.png")
