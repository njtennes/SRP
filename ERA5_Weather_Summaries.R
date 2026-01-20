WIND <- read.csv('/Users/nicktennes/Documents/ERA5 Weather Files CLEAN/ALLSITES_wind.csv')

wind_annual_summary <- WIND %>%
  group_by(site) %>%
  summarise(
    n_obs = sum(!is.na(speed)),
    
    temperature_mean = mean(temperature, na.rm = TRUE),
    
    pressure_mean = mean(pressure, na.rm = TRUE),
    
    speed_mean = mean(speed, na.rm = TRUE),
    speed_var  = var(speed,  na.rm = TRUE),
    speed_min  = min(speed,  na.rm = TRUE),
    speed_max  = max(speed,  na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(site)
