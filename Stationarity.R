library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(patchwork)

# Read data
Xdf <- read.csv("/Users/nicktennes/Desktop/FullOutputMatrix.csv")

# Rename year column
Xdf <- Xdf %>%
  rename(Year = X)

# Site columns only
site_cols <- setdiff(names(Xdf), c("Year", "hour_of_year", "hour"))

# Long format: one row per site-hour
long <- Xdf %>%
  pivot_longer(
    cols = all_of(site_cols),
    names_to = "Site",
    values_to = "MW"
  )

# Annual mean output by site
annual_mean <- long %>%
  group_by(Site, Year) %>%
  summarize(
    mean_MW = mean(MW, na.rm = TRUE),
    .groups = "drop"
  )

trend_results <- annual_mean %>%
  group_by(Site) %>%
  do(tidy(lm(mean_MW ~ Year, data = .))) %>%
  ungroup() %>%
  filter(term == "Year") %>%
  select(Site, trend_MW_per_year = estimate, p_value = p.value)

trend_results <- trend_results %>%
mutate(
  Sig = case_when(
    p_value < 0.01 ~ "***",
    p_value < 0.05 ~ "**",
    p_value < 0.10 ~ "*",
    TRUE ~ ""
  )
)

t.test(results$MW_per_year)

ggplot(annual_mean, aes(x = Year, y = mean_MW)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ Site, scales = "free_y") +
  labs(
    x = "Year",
    y = "Annual Mean Output (MW)",
    title = "Annual Mean Output Trends by Site"
  )

Xdf <- Xdf %>%
  mutate(
    Decade = case_when(
      Year < 1990 ~ "1980s",
      Year < 2000 ~ "1990s",
      Year < 2010 ~ "2000s",
      Year < 2020 ~ "2010s",
      TRUE        ~ "2020s"
    )
  )

site_cols <- setdiff(
  names(Xdf),
  c("Year", "Decade", "hour", "hour_of_year")
)

long <- Xdf %>%
  pivot_longer(
    cols = all_of(site_cols),
    names_to = "Site",
    values_to = "MW"
  )

decade_mean <- long %>%
  group_by(Site, Decade) %>%
  summarize(
    Mean_MW = mean(MW, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(decade_mean,
       aes(x = Decade,
           y = Mean_MW,
           group = Site,
           color = Site)) +
  geom_line() +
  geom_point(size = 3)

decade_summary <- decade_mean %>%
  filter(Decade %in% c("1980s", "2010s")) %>%
  pivot_wider(
    names_from = Decade,
    values_from = Mean_MW
  ) %>%
  mutate(
    MW_Change = `2010s` - `1980s`,
    Pct_Change = 100 * MW_Change / `1980s`
  ) %>%
  arrange(desc(Pct_Change))

annual_mean <- annual_mean %>%
  mutate(
    Decade = case_when(
      Year < 1990 ~ "1980s",
      Year < 2000 ~ "1990s",
      Year < 2010 ~ "2000s",
      Year < 2020 ~ "2010s",
      TRUE        ~ "2020s"
    )
  )

results_anova <- data.frame()

for(site in unique(annual_mean$Site)) {
  
  df_site <- annual_mean %>%
    filter(Site == site)
  
  mod <- aov(mean_MW ~ Decade, data = df_site)
  
  pval <- summary(mod)[[1]][["Pr(>F)"]][1]
  
  results_anova <- rbind(
    results_anova,
    data.frame(
      Site = site,
      Pvalue = pval
    )
  )
}

results_anova
