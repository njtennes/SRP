library(tidyverse)

data_dir <- "/Users/nicktennes/Desktop/GCJ"
files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

read_weather_year <- function(f) {
  
  year <- stringr::str_extract(basename(f), "\\d{4}")
  
  # Row 1-2 metadata; row 3 is the header
  raw <- read.csv(
    f,
    skip = 2,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Trim whitespace in column names, fix blanks, then make names unique
  nm <- trimws(names(raw))
  nm[nm == "" | is.na(nm)] <- "X"
  names(raw) <- make.unique(nm, sep = "_")
  
  # Convert to numeric (safe if everything is numeric-like)
  raw <- raw %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
  
  raw %>%
    mutate(
      hour = row_number(),
      year = as.integer(year)
    ) %>%
    relocate(year, hour)
}

weather_long <- purrr::map_dfr(files, read_weather_year)

####
yearly_means <- weather_long %>%
  group_by(year) %>%
  summarise(
    mean_GHI  = mean(ghi,  na.rm = TRUE),
    mean_DNI  = mean(dni,  na.rm = TRUE),
    mean_Tdry = mean(tdry, na.rm = TRUE),
    mean_Pres = mean(pres, na.rm = TRUE),
    mean_wspd = mean(wspd, na.rm = TRUE)
  ) %>%
  ungroup()

yearly_means <- yearly_means %>%
  mutate(
    era = if_else(year < 1980, "Pre-1980", "Post-1980")
  )

t_GHI  <- t.test(mean_GHI  ~ era, data = yearly_means)
t_DNI  <- t.test(mean_DNI  ~ era, data = yearly_means)
t_Tdry <- t.test(mean_Tdry ~ era, data = yearly_means)
t_Pres <- t.test(mean_Pres ~ era, data = yearly_means)
t_wspd <- t.test(mean_wspd ~ era, data = yearly_means)

t_GHI
t_DNI
t_Tdry
t_Pres
t_wspd

library(broom)

T_table <- bind_rows(
  tidy(t_GHI)  %>% mutate(variable = "GHI"),
  tidy(t_DNI)  %>% mutate(variable = "DNI"),
  tidy(t_Tdry) %>% mutate(variable = "Temperature"),
  tidy(t_Pres) %>% mutate(variable = "Pressure"),
  tidy(t_wspd) %>% mutate(variable = "Wind Speed")
) %>%
  select(variable, estimate, statistic, p.value, conf.low, conf.high)

