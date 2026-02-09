library(dplyr)
library(gt) 
library(knitr)
library(kableExtra)

WIND <- read.csv('/Users/nicktennes/Documents/ERA5 Weather Files CLEAN/ALLSITES_wind.csv')

WIND <- WIND %>%
filter(year >= 1988)

wind_summary <- WIND %>%
  group_by(site) %>%
  summarise(
    n_obs = sum(!is.na(speed)),
    
    temperature_mean = mean(temperature, na.rm = TRUE),
    speed_mean = mean(speed, na.rm = TRUE),
    speed_sd  = sd(speed,  na.rm = TRUE),
    speed_min  = min(speed,  na.rm = TRUE),
    speed_max  = max(speed,  na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(site)

wind_yearly_means <- WIND %>%
  group_by(site, year) %>%
  summarise(
    temperature_mean = mean(temperature, na.rm = TRUE),
    pressure_mean    = mean(pressure,    na.rm = TRUE),
    speed_mean       = mean(speed,       na.rm = TRUE),
    .groups = "drop"
  )

wind_yearly_means <- wind_yearly_means %>%
  mutate(
    temperature_mean = as.numeric(temperature_mean),
    pressure_mean    = as.numeric(pressure_mean),
    speed_mean       = as.numeric(speed_mean)
  )

wind_mean_of_yearly_means <- wind_yearly_means %>%
  group_by(site) %>%
  summarise(
    n_years = n(),
    
    temperature_mean2 = mean(temperature_mean, na.rm = TRUE),
    temperature_var2  = var(temperature_mean),
    temperature_min2  = min(temperature_mean,  na.rm = TRUE),
    temperature_max2  = max(temperature_mean,  na.rm = TRUE),
    
    pressure_mean2 = mean(pressure_mean, na.rm = TRUE),
    pressure_var2  = var(pressure_mean,  na.rm = TRUE),
    
    speed_mean2 = mean(speed_mean, na.rm = TRUE),
    speed_var2  = var(speed_mean,  na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(site)

SOLAR <- read.csv('/Users/nicktennes/Documents/ERA5 Weather Files CLEAN/ALLSITES_solar.csv')

SOLAR <- SOLAR %>%
  filter(year >= 1988)

solar_daylight <- SOLAR %>%
  filter(!between(hour, 9, 5) )

solar_hour_summary <- SOLAR %>%
  group_by(hour) %>%
  summarise(
    n_obs = n(),
    
    ghi_mean = mean(ghi, na.rm = TRUE),
    dni_mean = mean(dni, na.rm = TRUE),
    dhi_mean = mean(dhi)
  
  ) %>%
  arrange(hour)

solar_site_summary <- SOLAR %>%
  group_by(site) %>%
  summarise(
    n_obs = n(),
    
    ghi_mean = mean(ghi, na.rm = TRUE),
    dni_mean = mean(dni, na.rm = TRUE),
    dhi_mean = mean(dhi)
  ) %>%
  arrange(site)

solar_site_summary <- SOLAR %>%
  group_by(site) %>%
  summarise(
    n_obs = n(),
    temp_avg = mean(tdry),
    ghi_mean = mean(ghi, na.rm = TRUE),
    ghi_sd = sd(ghi),
    ghi_max = max(ghi),
    ghi_min = min(ghi)
  ) %>%
  arrange(site)

solar_site_summary_daylight <- solar_daylight %>%
  group_by(site) %>%
  summarise(
    n_obs = n(),
    
    ghi_mean = mean(ghi, na.rm = TRUE),
    ghi_sd = sd(ghi),
    max_ghi = max(ghi),
    min_ghi = min(ghi),
  ) %>%
  arrange(site)
