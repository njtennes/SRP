library(tidyverse)
library(ggplot2)

solar <- read.csv("/Users/nicktennes/Desktop/ALLSITES_solar.csv")

glimpse(solar)

# 1) Collapse to site-year means (this is your "yearly_means", but for every site)
site_year_means <- solar %>%
  group_by(site, year) %>%
  summarise(
    mean_DNI = mean(dni, na.rm = TRUE),
    mean_GHI = mean(ghi, na.rm = TRUE),
    mean_Tdry = mean(tdry, na.rm = TRUE),
    mean_Pres = mean(pres, na.rm = TRUE),
    mean_wspd = mean(wspd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    era = if_else(year < 1980, "Pre-1980", "Post-1980")
  )

site_labels <- c(
  CasaGrande = "Casa Grande, AZ",
  Deming_NM = "Deming, NM",
  Kingman_AZ = "Kingman, AZ",
  GrandCanyonJunction = "GC Junction, AZ",
  St.Johns_AZ = "St. Johns, AZ",
  Wilcox_AZ = "Willcox, AZ"
)

# 2) Panel figure: Mean DNI by year for all sites
ggplot(site_year_means, aes(x = year, y = mean_DNI)) +
  geom_point(size = 1.6, alpha = 0.8) +
  geom_smooth(
    method = "lm",
    formula = y ~ poly(x - 1980, 2, raw = TRUE),
    se = FALSE,
    linewidth = 1.0,
    color = "lightpink4"
  ) +
  geom_vline(xintercept = 1988, linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ site, labeller = labeller(site = site_labels), scales = "free_y") +
  labs(x = "Year", y = "Mean Annual DNI") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
n