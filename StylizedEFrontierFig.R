# A Fake Frontier from GPT

library(ggplot2)
library(dplyr)

set.seed(123)

# -------------------------------------------------
# 1. Create stylized efficient frontier (concave)
# -------------------------------------------------

risk <- seq(0, 10, length.out = 200)

# concave increasing function
frontier <- tibble(
  risk = risk,
  return = 4 * log(risk + 0.2)   # concave, increasing
)

# -------------------------------------------------
# 2. Create dominated single-site portfolios
# -------------------------------------------------

n_sites <- 10

single_sites <- tibble(
  risk = runif(n_sites, 0, 10)
) %>%
  mutate(
    # get frontier value at each risk
    frontier_val = 4 * log(risk + 0.5),
    # push points downward so they are dominated
    return = frontier_val - 3*runif(n_sites, 0.5, 2)
  )

# -------------------------------------------------
# 3. Plot
# -------------------------------------------------

ggplot() +
  geom_line(data = frontier, aes(risk, return),
            linewidth = 1.3, color = "blue3") +
  geom_point(data = single_sites, aes(risk, return),
             size = 2.5, alpha = 0.8, color = "red3") +
  labs(
    x = "Risk",
    y = "Return",
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )
