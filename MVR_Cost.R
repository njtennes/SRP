# ============================================
# Libraries
# ============================================
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(CVXR)
library(kableExtra)

#=================================

site_files <- c(
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/GC Junction Wind.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/Encino Wind.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/Medicine Bow Wind.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/Silver City Wind.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/Casa Grande Solar.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/Deming Solar.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/GC Junction Solar.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/Kingman Solar.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/St Johns Solar.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/Wilcox Solar.csv'
)
# ============================================
# gamma (1st one for single optimization)
# ============================================
gamma0     <- 1.5
gamma_grid <- 10^seq(-10, 3, length.out = 50)

site_names <- basename(site_files) |> tools::file_path_sans_ext()

site_meta <- tibble(
  file     = site_files,
  sitename = site_names
) %>%
  mutate(
    tech = case_when(
      str_detect(sitename, regex("wind",  ignore_case = TRUE)) ~ "Wind",
      str_detect(sitename, regex("solar", ignore_case = TRUE)) ~ "Solar",
      TRUE ~ "Unknown"
    ))

read_site_long_output <- function(path, sitename) {
  
  df <- readr::read_csv(path, show_col_types = FALSE)
  stopifnot("Hour" %in% names(df))
  df$Hour <- as.integer(df$Hour)
  
  df_long <- df |>
    tidyr::pivot_longer(
      -Hour,
      names_to  = "Year",
      values_to = "Output_kW"
    ) |>
    dplyr::mutate(
      Year      = as.integer(Year),
      Site      = sitename,
      Output_kW = as.numeric(Output_kW)
    ) |>
    dplyr::select(Year, Hour, Site, Output_kW)
  
  df_long
}

panel_out <- purrr::map2_dfr(
  site_files, site_names,
  ~ read_site_long_output(.x, .y)
)

# ============================================
# Build X matrix (hours x sites) amd the time data
# ============================================
build_X <- function(panel_df, value_col = "Output_kW") {
  W <- panel_df |>
    dplyr::arrange(Year, Hour, Site) |>
    tidyr::pivot_wider(names_from = Site, values_from = all_of(value_col)) |>
    dplyr::arrange(Year, Hour)
  
  site_cols <- setdiff(names(W), c("Year", "Hour"))
  
  list(
    X     = as.matrix(W[ , site_cols]),
    sites = site_cols,
    meta  = W[, c("Year", "Hour")]
  )
}

bx    <- build_X(panel_out, value_col = "Output_kW")
X     <- bx$X
sites <- bx$sites
meta  <- bx$meta

#write.csv(X, file = "/Users/nicktennes/Desktop/FullOutputMatrix.csv", row.names = FALSE)

cost_per_MWh <- c(
  "Medicine Bow Wind" = 60,
  "Encino Wind"       = 52.5,
  "Silver City Wind"  = 52.5,
  "GC Junction Wind"  = 51,
  "Kingman Solar"     = 29,
  "GC Junction Solar" = 29,
  "Casa Grande Solar" = 29,
  "Wilcox Solar"      = 29,
  "St Johns Solar"    = 29,
  "Deming Solar"      = 29
)
cost_per_MWh <- cost_per_MWh[sites]
# ============================================
# Derive HourOfDay, DayOfYear, Month
# ============================================
meta <- meta %>%
  mutate(
    HourOfDay = ((Hour - 1) %% 24),        # 0–23
    DayOfYear = ((Hour - 1) %/% 24) + 1    # 1–365
  ) %>%
  mutate(
    Month = case_when(
      DayOfYear <= 31                                      ~ 1,
      DayOfYear <= 31 + 28                                 ~ 2,
      DayOfYear <= 31 + 28 + 31                            ~ 3,
      DayOfYear <= 31 + 28 + 31 + 30                       ~ 4,
      DayOfYear <= 31 + 28 + 31 + 30 + 31                  ~ 5,
      DayOfYear <= 31 + 28 + 31 + 30 + 31 + 30             ~ 6,
      DayOfYear <= 31 + 28 + 31 + 30 + 31 + 30 + 31        ~ 7,
      DayOfYear <= 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31   ~ 8,
      DayOfYear <= 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 ~ 9,
      DayOfYear <= 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 ~ 10,
      DayOfYear <= 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30 ~ 11,
      TRUE ~ 12
    )
  )

# ============================================
# Mean & covariance of CF 
# ============================================
X_MW <- X / 1000

# Cost-adjusted output series: (MW) / ($/MWh)  ~  MWh per $
X_econ <- sweep(X_MW, 2, cost_per_MWh, FUN = "/")

Sigma <- cov(X_econ)
mu    <- colMeans(X_econ)
mu_output <- colMeans(X_MW)

# ============================================
#Site Summary Table
# ============================================
summary <- tibble(
  Site     = sites,
  Weighted_Output = mu,
  Weighted_Var = apply(X_econ, 2, var),
  Mean_output = mu_output
) |>
  dplyr::arrange(dplyr::desc(Weighted_Output))

kable(
  summary,
  format = "latex",
  booktabs = TRUE,
  caption = "Summary statistics for candidate renewable sites",
  digits = 2
) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down")
  )

con <- pipe("pbcopy", "w")
write.table(summary, con, sep = "\t", col.names = NA)
close(con)

summary_table <- as.data.frame(X) |>
  pivot_longer(
    cols = everything(),
    names_to = "Site",
    values_to = "Output"
  ) |>
  group_by(Site) |>
  summarise(
    Mean = mean(Output, na.rm = TRUE),
    StDev = sd(Output, na.rm = TRUE),
    Min = min(Output, na.rm = TRUE),
    Max = max(Output, na.rm = TRUE),
    .groups = "drop"
  )

summary_table

m <- meta$Month
y <- meta$Year

season_vec <- dplyr::case_when(
  m %in% c(12, 1, 2) ~ "Winter",
  m %in% c(3, 4, 5)  ~ "Spring",
  m %in% c(6, 7, 8)  ~ "Summer",
  m %in% c(9,10,11)  ~ "Fall",
  TRUE               ~ NA_character_
)

summary_monthly <- as_tibble(X) |>
  mutate(
    Month = m,
    Year  = y  
  ) |>
  pivot_longer(
    cols = -c(Month, Year),
    names_to  = "Site",
    values_to = "CF"
  ) |>
  group_by(Year, Month, Site) |>
  summarise(
    Yearly_Mean = (mean(CF, na.rm = TRUE) * 100),
    .groups = "drop"
  ) |>
  group_by(Month, Site) |>
  summarise(
    Mean_CF = mean(Yearly_Mean),
    SD_CF   = sd(Yearly_Mean),
    N_years = n(),
    SE_CF   = SD_CF / sqrt(N_years),
    .groups = "drop"
  )

summary_monthly <- summary_monthly |>
  mutate(
    tech = ifelse(
      grepl("solar", Site, ignore.case = TRUE),
      "Solar",
      "Wind"
    )
  )

site_order <- summary_monthly |>
  distinct(Site, tech) |>
  arrange(tech, Site) |>
  pull(Site)

# Apply factor order
summary_monthly <- summary_monthly |>
  mutate(
    Site = factor(Site, levels = site_order)
  )

summary_monthly <- summary_monthly |>
  mutate(
    tech = factor(tech, levels = c("Solar", "Wind"))
  )
# Get site lists by tech
solar_sites <- summary_monthly |>
  filter(tech == "Solar") |>
  distinct(Site) |>
  pull(Site)

wind_sites <- summary_monthly |>
  filter(tech == "Wind") |>
  distinct(Site) |>
  pull(Site)

# Generate shades (base R, no extra packages)
solar_colors <- colorRampPalette(c("#F4A6A6", "#8B0000"))(length(solar_sites))
wind_colors  <- colorRampPalette(c("#A6C8FF", "#08306B"))(length(wind_sites))

# Named vector for scale_color_manual()
site_colors <- c(
  setNames(solar_colors, solar_sites),
  setNames(wind_colors,  wind_sites)
)

ggplot(summary_monthly,
       aes(x = Month, y = Mean_CF, color = Site)) +
  geom_point(size = 2.6, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = Mean_CF - SE_CF,
        ymax = Mean_CF + SE_CF),
    width = 0.25,
    alpha = 0.6
  ) +
  scale_color_manual(values = site_colors) +
  scale_x_continuous(
    breaks = 1:12,
    labels = month.abb
  ) +
  labs(
    x = "Month",
    y = "Expected Capacity Factor (%)",
    color = "Site"
  ) +
  theme_minimal(base_size = 12)

h <- meta$HourOfDay

hourly_tech <- as_tibble(X) |>
  mutate(Hour = h) |>
  pivot_longer(
    cols = -Hour,
    names_to  = "Site",
    values_to = "CF"
  ) |>
  left_join(
    site_meta |> select(sitename, tech),
    by = c("Site" = "sitename")
  ) |>
  group_by(Hour, tech) |>
  summarise(
    Mean_CF = mean(CF, na.rm = TRUE),
    SD_CF   = sd(CF, na.rm = TRUE),
    N       = sum(!is.na(CF)),
    SE_CF   = SD_CF / sqrt(N),
    .groups = "drop"
  ) |>
  arrange(tech, Hour)

hourly_tech <- hourly_tech %>%
  mutate(Mean_CF = Mean_CF * 100)

ggplot(hourly_tech, aes(x = Hour, y = Mean_CF, color = tech)) +
  geom_line(size = 1, alpha = 1) +
  geom_errorbar(
    aes(ymin = Mean_CF - SE_CF,
        ymax = Mean_CF + SE_CF),
    width = 0.25,
    alpha = 0.6
  ) +
  scale_x_continuous(breaks = 0:23) +
  scale_color_manual(values = c("Solar" = "lightpink3", "Wind" = "#08306B")) +
  labs(x = "Hour of Day", y = "Expected Capacity Factor (%)", color = "Technology") +
  theme_minimal(base_size = 12)

# ============================================
# Correlation matrices
# ============================================
cor_X <- cor(X, use = "pairwise.complete.obs")
print(round(cor_X, 3))

con <- pipe("pbcopy", "w")
write.table(round(cor_X, 3), con, sep = "\t", col.names = NA)
close(con)

WindSolar <- X[, c("Medicine Bow Wind", "Encino Wind", "Silver City Wind", 
                   "GC Junction Wind", "Kingman Solar","GC Junction Solar",
                   "Casa Grande Solar", "Wilcox Solar", "St Johns Solar", 
                   "Deming Solar")]

cor_WindSolar <- cor(WindSolar, use = "pairwise.complete.obs")
round(cor_WindSolar, 3)

con <- pipe("pbcopy", "w")
write.table(round(cor_WindSolar, 3), con, sep = "\t", col.names = NA)
close(con)

#wind_only <- X[, c("Encino Wind", "Medicine Bow Wind", "Silver City Wind", "GC Junction Wind")]
#cor_wind <- cor(wind_only, use = "pairwise.complete.obs")
#round(cor_wind, 3)

#GCJ_only <- X[, c("GC Junction Solar", "GC Junction Wind")]
#cor_GCJ <- cor(GCJ_only, use = "pairwise.complete.obs")
#round(cor_GCJ, 3)

#WindSolar <- X[, c("Medicine Bow Wind", "Encino Wind", "Wilcox Solar", "Kingman Solar")]
#cor_WindSolar <- cor(WindSolar, use = "pairwise.complete.obs")
#round(cor_WindSolar, 3)

# ============================================
# Markowitz solver
# max  μᵀw - γ wᵀΣw  s.t. w >= 0, Σw = 1
# ============================================
solve_markowitz <- function(mu, Sigma, gamma = gamma0, solver = "OSQP") {
  n <- length(mu)
  w <- CVXR::Variable(n)
  obj <- t(mu) %*% w - gamma * CVXR::quad_form(w, Sigma)
  cons <- list(w >= 0, sum(w) == 1)
  
  prob <- CVXR::Problem(CVXR::Maximize(obj), cons)
  res  <- CVXR::solve(prob, solver = solver)
  
  wv <- as.numeric(res$getValue(w))
  list(
    status     = res$status,
    w          = wv,
    exp_output = sum(mu * wv),
    variance   = as.numeric(t(wv) %*% Sigma %*% wv),
    gamma      = gamma
  )
}

# ============================================
# Single-γ solution
# ============================================
sol0 <- solve_markowitz(mu, Sigma, gamma = gamma0)
cat("\n--- Single-γ solution ---\n")
print(
  tibble(Site = sites, Weight = round(sol0$w, 4)) |>
    dplyr::arrange(dplyr::desc(Weight))
)
cat(sprintf("γ=%g | E[MWh/$]=%.4f | Var=%.4f\n", sol0$gamma, sol0$exp_output, sol0$variance))

weights_table <- tibble(Site = sites, Weight = round(sol0$w, 4))

con <- pipe("pbcopy", "w")
write.table(weights_table, con, sep = "\t", col.names = NA)
close(con)

# ============================================
# Efficient frontier along gamma grid
# ============================================
front <- lapply(gamma_grid, function(g) solve_markowitz(mu, Sigma, gamma = g))

frontier <- tibble(
  gamma    = sapply(front, `[[`, "gamma"),
  exp_mwhdollar  = sapply(front, `[[`, "exp_output"),
  variance = sapply(front, `[[`, "variance")
) |>
  dplyr::arrange(variance)

cat("\n--- Efficient frontier (full year, head) ---\n")
print(head(frontier, 20))

# ggplot frontier + single-site points
frontier_labeled <- frontier %>%
  mutate(label_gamma = if_else(row_number() %% 5 == 0,
                               paste0("γ=", round(gamma, 2)),
                               NA_character_))

ggplot() +
  geom_line(data = frontier, aes(x = variance, y = exp_mwhdollar),
            linewidth = 0.7, color = "darkblue") +
  geom_point(data = frontier, aes(x = variance, y = exp_mwhdollar),
             size = 1.5, color = "darkblue", alpha = 0.7) +
  geom_point(data = summary, aes(x = Weighted_Var, y = Weighted_Output),
             color = "red", size = 2) +
  geom_text_repel(
    data = frontier_labeled,
    aes(x = variance, y = exp_mwhdollar, label = label_gamma),
    size = 3,
    color = "darkblue",
    na.rm = TRUE,
    max.overlaps = Inf
  ) +
  geom_text_repel(
    data = summary,
    aes(x = Weighted_Var, y = Weighted_Output, label = Site),
    size = 3,
    color = "red"
  ) +
  labs(
    x = "Variance MWHr^2/$ spent",
    y = "Expected MWHr^2/$ spent", 
    title = "Cost Weighted Efficient Frontier"
  ) +
  theme_minimal(base_size = 12)

###### ============================================
# Reliability-hour subset (custom windows)
# May–Sep 3–7pm, Dec–Mar 6–9am and 6–9pm
# ============================================

idx_rel <- which(
  # Summer late afternoon peak: May–Sep, 3–7pm
  (meta$Month %in% 6:8  & meta$HourOfDay %in% 16:20) #|
    # Winter morning peak: Dec–Mar, 6–9am
   # (meta$Month %in% c(12, 1, 2, 3) & meta$HourOfDay %in% 6:9) |
    # Winter evening peak: Dec–Mar, 6–9pm
    #(meta$Month %in% c(12, 1, 2, 3) & meta$HourOfDay %in% 18:21)
)

X_econ_rel <- X_econ[idx_rel, , drop = FALSE]

mu_rel    <- colMeans(X_econ_rel, na.rm = TRUE)
Sigma_rel <- stats::cov(X_econ_rel, use = "pairwise.complete.obs")

# Single-γ solution for reliability hours
sol_rel <- solve_markowitz(mu_rel, Sigma_rel, gamma = gamma0)

cat("\n=== Reliability hours: Single-γ solution ===\n")
print(
  tibble(Site = sites, Weight = round(sol_rel$w, 4)) |>
    dplyr::arrange(dplyr::desc(Weight))
)
cat(sprintf("status=%s | γ=%g | E[CF]_rel=%.4f | Var_rel=%.4f | n_hours=%d\n",
            sol_rel$status, sol_rel$gamma, sol_rel$exp_CF, sol_rel$variance, nrow(X_rel)))

fronttable_rel <- tibble(Site = sites, Weight = round(sol_rel$w, 4))

con <- pipe("pbcopy", "w")
write.table(fronttable_rel, con, sep = "\t", col.names = NA)
close(con)

# Frontier for reliability hours
front_rel <- lapply(gamma_grid, function(g) solve_markowitz(mu_rel, Sigma_rel, gamma = g))

frontier_rel <- tibble(
  gamma    = sapply(front_rel, `[[`, "gamma"),
  exp_CF   = sapply(front_rel, `[[`, "exp_CF"),
  variance = sapply(front_rel, `[[`, "variance"),
  status   = sapply(front_rel, `[[`, "status")
) |>
  dplyr::filter(status == "optimal") |>
  dplyr::arrange(variance)

cat("\n=== Reliability hours: Efficient frontier (head) ===\n")
print(head(frontier_rel, 20))

plot(frontier_rel$variance, frontier_rel$exp_CF, pch = 16,
     xlab = "Variance of CF (reliability hours)", ylab = "Expected CF (reliability hours)",
     main = "Efficient Frontier — Reliability Hours")
grid()

# Reliability-hour site stats + correlation (optional, mirrors your JA section)
site_stats_rel <- tibble(
  Site     = sites,
  Mean_CF  = mu_rel,
  Variance = apply(X_rel, 2, var, na.rm = TRUE)
) |>
  dplyr::arrange(dplyr::desc(Mean_CF))

cat("\n=== Reliability hours: Site stats ===\n")
print(site_stats_rel)

cor_X_rel <- cor(X_rel, use = "pairwise.complete.obs")
print(round(cor_X_rel, 3))
