library(tidyverse)

hrs <- 0:23

df <- tibble(
  hour = hrs,
  load = c(4200, 4050, 3950, 3900, 3950, 4100,
           4500, 5000, 5600, 6100, 6500, 6800,
           7000, 7100, 7300, 7600, 8000, 8350,
           8500, 8300, 7800, 7000, 6000, 5000),
  nuclear = rep(1200, 24),
  coal = c(1000, 1000, 1000, 1000, 1000, 1000,
           880, 880, 870, 860, 850, 850,
           840, 840, 830, 820, 800, 780,
           760, 860, 880, 920, 960, 1000),
  solar = c(0, 0, 0, 0, 0, 30,
            150, 450, 900, 1400, 1800, 2100,
            2400, 2500, 2150, 100, 2400, 1900,
            800, 200, 50, 0, 0, 0),
  wind = c(450, 750, 740, 410, 500, 380,
           350, 250, 300, 280, 260, 250,
           260, 170, 80, 320, 100, 500,
           500, 600, 500, 660, 800, 750)
) %>%
  mutate(
    residual_after_basics = load - (nuclear + coal + solar + wind),
    gas_cc = pmax(0, pmin(residual_after_basics, 3250)),
    gas_peaker = pmax(0, residual_after_basics - gas_cc)
  )

plot_df <- df %>%
  pivot_longer(
    cols = c(nuclear, coal, gas_cc, solar, wind, gas_peaker),
    names_to = "resource",
    values_to = "mw"
  ) %>%
  mutate(
    resource = factor(
      resource,
      levels = c("nuclear", "coal", "gas_cc", "solar", "wind", "gas_peaker"),
      labels = c("Nuclear", "Coal", "Combined Cycle", "Solar", "Wind", "Peakers")
    )
  ) %>%
  arrange(hour, resource)

ggplot(plot_df, aes(x = hour, y = mw, fill = resource, group = resource)) +
  geom_area(position = position_stack(reverse = TRUE), alpha = 0.95) +
  geom_line(
    data = df,
    aes(x = hour, y = load),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Nuclear" = "darkgreen",
      "Coal" = "gray50",
      "Combined Cycle" = "steelblue",
      "Solar" = "orange",
      "Wind" = "skyblue",
      "Peakers" = "slateblue3"
    ),
    breaks = c("Nuclear", "Coal", "Combined Cycle", "Solar", "Wind", "Peakers")
  ) +
  scale_x_continuous(breaks = seq(0, 23, by = 2)) +
  labs(
    x = "Hour of Day",
    y = "MW",
    fill = NULL,
    title = "Stylized Summer Load Shape and Resource Stack",
    subtitle = "*Not Real Data*"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

df2 <- df %>%
  mutate(
    solar_new = 3 * solar,
    wind_new  = 3 * wind
  )

# Calculate total added renewable energy
extra_energy <- sum(df2$solar_new + df2$wind_new - (df$solar + df$wind))

# Spread evenly across hours
extra_per_hour <- extra_energy / 24

df2 <- df2 %>%
  mutate(
    load_new = load + extra_per_hour
  )

df2 <- df2 %>%
  mutate(
    residual = load_new - (nuclear + coal + solar_new + wind_new),
    
    # Combined cycle capped
    gas_cc_new = pmax(0, pmin(residual, 3250)),
    
    residual_after_cc = residual - gas_cc_new,
    
    # Peakers capped
    gas_peaker_new = pmax(0, pmin(residual_after_cc, 2200)),
    
    # Anything left is unserved
    unserved = pmax(0, residual_after_cc - gas_peaker_new)
  )

df2 <- df2 %>%
  mutate(
    total_gen = nuclear + coal + solar_new + wind_new,
    curtailment = pmax(0, total_gen - load_new)
  )

plot_df2 <- df2 %>%
  select(hour, nuclear, coal, gas_cc_new, solar_new, wind_new, gas_peaker_new, unserved, curtailment) %>%
  pivot_longer(
    cols = -hour,
    names_to = "resource",
    values_to = "mw"
  ) %>%
  mutate(
    resource = factor(
      resource,
      levels = c("nuclear", "coal", "gas_cc_new", "solar_new", "wind_new", "gas_peaker_new", "curtailment", "unserved"),
      labels = c("Nuclear", "Coal", "Combined Cycle", "Solar", "Wind", "Peakers", "Curtailed Energy", "Unserved Energy")
    )
  )

ggplot(plot_df2, aes(x = hour, y = mw, fill = resource)) +
  geom_area(position = position_stack(reverse = TRUE), alpha = 0.95) +
  
  geom_line(
    data = df2,
    aes(x = hour, y = load_new),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  
  scale_fill_manual(values = c(
    "Nuclear" = "darkgreen",
    "Coal" = "gray50",
    "Combined Cycle" = "steelblue",
    "Solar" = "orange",
    "Wind" = "skyblue",
    "Peakers" = "slateblue3",
    "Curtailed Energy" = "red3",
    "Unserved Energy" = "black"
  )) +
  
  labs(
    x = "Hour of Day",
    y = "MW",
    subtitle = "Timing Mismatch Leads to Unserved Energy",
    title = "Scenario: Tripled Solar & Wind with Proportional Load Growth",
    fill = NULL
  ) +
  scale_y_continuous(
    breaks = seq(0, 12000, by = 2000),
    labels = scales::comma
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
  
