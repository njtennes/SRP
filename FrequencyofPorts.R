library(readxl)
library(tidyverse)
library(ggplot2)

Ports <- read_excel("/Users/nicktennes/Desktop/PortfolioValues.xlsx")
head(Ports)

Ports <- Ports %>%
  rename(`Medicine Bow Only` = MB) %>%
  rename('Kingman Only' = Kingman)

colMeans(Ports == 0)
colMeans(Ports < 15000)

frequency <- rbind(
  freq_zero = colMeans(Ports == 0)*100,
  freq_below_15k = colMeans(Ports < 15000)*100,
  freq_below_35k = colMeans(Ports < 35000)*100,
  freq_below_39.5k = colMeans(Ports < 39500)*100,
  freq_below_50k = colMeans(Ports < 50000)*100
)

Ports_long <- pivot_longer(Ports, everything(),
                           names_to = "portfolio",
                           values_to = "output")

Ports_long <- Ports_long %>%
  mutate('Portfolio Capacity Factor' = output/1000)

ggplot(Ports_long, aes(x = `Portfolio Capacity Factor`, fill = portfolio)) +
  geom_histogram(bins = 40) +
  scale_fill_manual(values = c(
    "Medicine Bow Only" = "lightblue2",
    "Kingman Only" = "coral3",
    "0.75" = "lightgreen",
    "1.5" = "darkolivegreen3",
    "3" = "darkolivegreen",
    "4.5" = "darkgreen"
  )) +
  facet_wrap(
    ~ portfolio,
    scales = "free",
    labeller = as_labeller(c(
      "0.75" = "gamma == 0.75",
      "1.5"  = "gamma == 1.5",
      "3"    = "gamma == 3",
      "4.5"  = "gamma == 4.5",
      "Medicine Bow Only" = "'Medicine Bow Only'",
      "Kingman Only" = "'Kingman Only'"
    ), label_parsed)
  ) +
  theme_minimal() +
  theme(legend.position = "none")
