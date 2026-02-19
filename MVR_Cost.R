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
gamma0     <- 2
gamma_grid <- 10^seq(-10, 3, length.out = 50)

site_names <- basename(site_files) |> tools::file_path_sans_ext()

#bahahaha
site_names[site_names == "Wilcox Solar"] <- "Willcox Solar"

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
  "Medicine Bow Wind" = 374.32,
  "Encino Wind"       = 265.66,
  "Silver City Wind"  = 203.47,
  "GC Junction Wind"  = 178.92,
  "Kingman Solar"     = 175.33,
  "GC Junction Solar" = 166.28,
  "Casa Grande Solar" = 122.43,
  "Willcox Solar"     = 165.09,
  "St Johns Solar"    = 171.76,
  "Deming Solar"      = 184.39
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

SigmaMW <- cov(X_MW)
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
    Mean_CF = mean(Yearly_Mean/100000),
    SD_CF   = sd(Yearly_Mean/100000),
    p1 = quantile(Mean_CF, 0.01, na.rm = TRUE),
    p99 = quantile(Mean_CF, 0.99, na.rm = TRUE),
    N_years = n(),
    SE_CF   = SD_CF / sqrt(N_years),
    .groups = "drop"
  )

nameplate_kW <- 1000  

summary_monthly <- as_tibble(X) |>
  mutate(
    Month = m,
    Year  = y
  ) |>
  pivot_longer(
    cols = -c(Month, Year),
    names_to  = "Site",
    values_to = "kW"
  ) |>
  mutate(
    CF = kW / nameplate_kW   # convert to CF once, cleanly
  ) |>
  group_by(Year, Month, Site) |>
  summarise(
    CF_year = mean(CF, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(Month, Site) |>
  summarise(
    Mean_CF  = mean(CF_year, na.rm = TRUE),
    SD_CF    = sd(CF_year, na.rm = TRUE),
    p5       = quantile(CF_year, 0.05, na.rm = TRUE),
    p95      = quantile(CF_year, 0.95, na.rm = TRUE),
    N_years  = n(),
    SE_CF    = SD_CF / sqrt(N_years),
    .groups  = "drop"
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
wind_colors  <- colorRampPalette(c("lightblue", "#08306B"))(length(wind_sites))

# Named vector for scale_color_manual()
site_colors <- c(
  setNames(solar_colors, solar_sites),
  setNames(wind_colors,  wind_sites)
)

ggplot(summary_monthly,
       aes(x = Month, y = Mean_CF, color = Site)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_ribbon(aes(ymin = p5, ymax = p95), alpha = 0.2) +
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

hourly_tech_summer <- as_tibble(X) |>
  mutate(
    Hour  = h,
    Month = m   # vector aligned to rows of X
  ) |>
  filter(Month %in% 6:8) |>
  select(-Month) |>
  pivot_longer(
    cols = -Hour,
    names_to  = "Site",
    values_to = "CF"
  ) |>
  left_join(site_meta |> select(sitename),
            by = c("Site" = "sitename")) |>
  group_by(Hour, Site) |>
  summarise(
    Mean_CF = mean(CF/1000, na.rm = TRUE),
    SD_CF   = sd(CF/1000, na.rm = TRUE),
    p25 = quantile(CF/1000, 0.25, na.rm = TRUE),
    p75 = quantile(CF/1000, 0.75, na.rm = TRUE),
    p5 = quantile(CF/1000, 0.05, na.rm = TRUE),
    p95 = quantile(CF/1000, 0.95, na.rm = TRUE),
    N       = sum(!is.na(CF)),
    .groups = "drop"
  ) |>
  arrange(Site, Hour)

hourly_tech_summer <- hourly_tech_summer |>
  mutate(
    tech = ifelse(
      grepl("solar", Site, ignore.case = TRUE),
      "Solar",
      "Wind"
    )
  )

site_order2 <- hourly_tech_summer |>
  distinct(Site, tech) |>
  arrange(tech, Site) |>
  pull(Site)

hourly_tech_summer <- hourly_tech_summer |>
  mutate(
    Site = factor(Site, levels = site_order2)
  )

hourly_tech_summer <- hourly_tech_summer |>
  mutate(
    tech = factor(tech, levels = c("Solar", "Wind"))
  )

ggplot(hourly_tech_summer, aes(y = Mean_CF, x = Hour, color = Site)) +
  geom_point(size= 1.5, alpha = 1) +
  geom_ribbon(aes(ymin = p25, ymax = p75), alpha = 0.2) +
  scale_color_manual(values = site_colors) +
labs(
  x = "Hour of Day (June-August)",
  y = "Expected Capacity Factor (%)",
  color = "Site"
) +
  theme_minimal(base_size = 12)

ggplot(hourly_tech_summer, aes(y = Mean_CF, x = Hour, color = Site)) +
  geom_point(size= 1.5, alpha = 1) +
  geom_ribbon(aes(ymin = p5, ymax = p95), alpha = 0.2) +
  scale_color_manual(values = site_colors) +
  labs(
    x = "Hour of Day (June-August)",
    y = "Expected Capacity Factor (%)",
    color = "Site"
  ) +
  theme_minimal(base_size = 12)

ggplot(hourly_tech, aes(x = Hour, y = Mean_CF, color = tech)) +
  geom_line(size = 1, alpha = 1) +
  scale_x_continuous(breaks = 0:23) +
  scale_color_manual(values = c("Solar" = "lightpink3", "Wind" = "#08306B")) +
  labs(x = "Hour of Day", y = "Expected Capacity Factor (%)", color = "Technology") +
  theme_minimal(base_size = 12)

kingman <- as_tibble(X[, "Kingman Solar"])

ggplot(data.frame(CF = kingman), aes(x = CF)) +
  geom_density(adjust = 1.2, na.rm = TRUE) +
  labs(
    title = "Distribution of Kingman Solar Output",
    x = "Output (kW)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12)

ggplot(data.frame(CF = kingman), aes(x = CF)) +
  stat_ecdf(geom = "step") +
  labs(
    title = "CDF of Kingman Solar Output",
    x = "Output (kW) or Capacity Factor",
    y = "F(x)"
  ) +
  theme_minimal(base_size = 12)

MB <- as_tibble(X[, "Medicine Bow Wind"])

MBD <- MB %>%
  mutate(hour = h)

MBD <- MBD %>%
  mutate(Day_ID = (row_number() - 1) %/% 24 + 1)

MBD <- MBD %>%
  group_by(Day_ID) %>%
  summarize(Output= sum(value))

MBD <- MBD %>%
  mutate(Output= Output/100000)


ggplot(MB, aes(x = value)) +
  geom_density(adjust = 1.2, na.rm = TRUE) +
  labs(
    title = "Distribution of Medicine Bow Wind Output",
    x = "Output (kW)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12)

ggplot(data.frame(CF = MB), aes(x = CF)) +
  stat_ecdf(geom = "step") +
  labs(
    title = "CDF of Medicine Bow Output",
    x = "Output (kW)",
    y = "F(x)"
  ) +
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
gamma0     <- 7.5

sol0 <- solve_markowitz(mu, Sigma, gamma = gamma0)
cat("\n--- Single-γ solution ---\n")
print(
  tibble(Site = sites, Weight = round(sol0$w, 4)) |>
    dplyr::arrange(dplyr::desc(Weight))
)
cat(sprintf("γ=%g | E[MWh^2/$]=%.4f | Var=%.4f\n", sol0$gamma, sol0$exp_output, sol0$variance))

weights_table <- tibble(Site = sites, Weight = round(sol0$w, 4))

con <- pipe("pbcopy", "w")
write.table(weights_table, con, sep = "\t", col.names = NA)
close(con)
#============
# Discrete Gamma Vector
#==============
gamma_vec <- c(0.01, .75, 1.125, 3, 100)

sols <- map(gamma_vec, ~solve_markowitz(mu, Sigma, gamma = .x, solver = "OSQP"))

gname <- function(x) format(x, trim = TRUE, scientific = FALSE)

moments_long <- tibble(
gamma = gamma_vec,
`E`   = map_dbl(sols, "exp_output"),
Var   = map_dbl(sols, "variance")
) |>
  pivot_longer(cols = c(`E`, Var), names_to = "Moment", values_to = "Value")

moments_tbl <- moments_long |>
  mutate(gamma = gname(gamma)) |>
  pivot_wider(names_from = gamma, values_from = Value) |>
  mutate(
    Moment = recode(Moment,
                    `E` = "E[MWh/$]",   
    )
  ) |>
  select(Moment, all_of(gname(gamma_vec)))

moments_tbl <- moments_tbl %>%
  rename("100" = "100.000")

moments_tbl <- moments_tbl %>%
  rename("3" = "3.000")

moments_tbl <- moments_tbl %>%
  rename("0.75" = "0.750")

moments_tbl <- moments_tbl %>%
  rename("0.01" = "0.010")

weights_long <- map2_dfr(
  sols, gamma_vec,
  ~tibble(gamma = .y, Site = sites, Weight = .x$w)
) %>%
  mutate(
    Weight = ifelse(abs(Weight) < 1e-6, 0, Weight)
  )

weights_tbl <- weights_long |>
  mutate(
    gamma = format(gamma, scientific = FALSE),
    WeightPct = round(100 * Weight, 1)
  ) |>
  select(Site, gamma, WeightPct) |>
  pivot_wider(names_from = gamma, values_from = WeightPct) |>
  arrange(Site)
# =================
# stacked bar chart of weights
# =================
site_order <- c("Casa Grande Solar", "Deming Solar", "Kingman Solar",
                "Encino Wind", "GC Junction Wind", "Medicine Bow Wind", "Silver City Wind")

site_colors <- c(
  # Solar
  "Casa Grande Solar" = "lightpink1",
  "Deming Solar"      = "lightpink3",
  "Kingman Solar"     = "indianred4",
  
  # Wind 
  "Encino Wind"       = "lightsteelblue2",
  "GC Junction Wind"  = "lightskyblue1",
  "Medicine Bow Wind" = "steelblue",
  "Silver City Wind"  = "lightsteelblue4"
)

weights_long_plot <- weights_tbl |>
  pivot_longer(
    cols = -Site,
    names_to = "gamma",
    values_to = "WeightPct"
  ) |>
  mutate(
    gamma = factor(gamma, levels = colnames(weights_tbl)[-1])
  )

weights_long_plot <- weights_long_plot |>
  left_join(site_meta |> select(sitename, tech),
            by = c("Site" = "sitename"))

weights_long_plot <- weights_long_plot |>
  group_by(Site) |>
  filter(sum(WeightPct, na.rm = TRUE) > 0)

weights_long_plot <- weights_long_plot |>
  mutate(
    tech = factor(tech, levels = c("Solar", "Wind"))
  ) |>
  arrange(tech, Site) |>
  mutate(
    Site = factor(Site, levels = unique(Site))
  )

ggplot(weights_long_plot, aes(x = gamma, y = WeightPct, fill = Site)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_manual(
    values = site_colors[site_order],   
    breaks = site_order,                
    drop   = TRUE
  ) +
  geom_text(
    aes(label = ifelse(WeightPct >= 2,
                       sprintf("%.1f%%", WeightPct),
                       "")),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  labs(x = expression("Risk Aversion Parameter, " * gamma), y = "Portfolio Weight (%)", fill = "Site") +
  theme_minimal(base_size = 12)

# ===================
# Plot of discrete gammas
# ===================

moments_plot <- moments_tbl |>
  pivot_longer(
    cols = -Moment,
    names_to = "gamma",
    values_to = "Value"
  ) |>
  pivot_wider(
    names_from = Moment,
    values_from = Value
  ) |>
  rename(
    Variance = `Var`,
    Return = 'E[MWh/$]'
  ) |>
  mutate(
    gamma = factor(gamma, levels = colnames(moments_tbl)[-1])
  )

ggplot(moments_plot,
      aes(x = Variance,
          y = Return)) +
  geom_point(size = 3, shape = 8) +
  geom_line(linewidth = .25, alpha = 5) +
  geom_text(
    aes(label = gamma),
    vjust = -1.5,
    size = 3
  ) +
  labs(
    x = expression("Variance (" * MW ~ "/" ~ "Hour x $ Millions Spent" * ")"),
    y = expression("Expected (" * MW ~ "/" ~ "Hour x $ Millions Spent" * ")")
  ) +
  theme_minimal(base_size = 12)

# ============================================
# Full efficient frontier
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

##MANUALL!!!
#frontier_labeled <- tibble::tibble(
#  gamma        = c(0.01, 0.75, 1.5, 3, 100),
#  variance     = c(1.78, 1.49, 0.48, 0.125, 0.105),
#  exp_mwhdollar = c(1.15, 1.123, 0.972, 0.75, 0.65)
#) |>
#  dplyr::mutate(
#   label_gamma = paste0("\u03B3 = ", gamma)  # γ =
#  )

ggplot() +
  geom_line(data = frontier, aes(x = variance, y = exp_mwhdollar),
            linewidth = 0.7, color = "darkblue") +
  geom_point(data = frontier, aes(x = variance, y = exp_mwhdollar),
             size = 1.5, color = "darkblue", alpha = 0.7) +
  geom_point(data = summary, aes(x = Weighted_Var, y = Weighted_Output),
             color = "red", size = 2) +
 #  geom_point(data = moments_plot,
 #             aes(x = Variance, y = Return),
 #             color = "red4",
 #             size = 2.5) +
 # geom_label_repel(
 #   data = moments_plot,
 #   aes(x = Variance, y = Return, label = gamma),
 #   size = 3,
 #   color = "black",
 #   label.size = 0,
 #   box.padding = 0.4,
 #   point.padding = 0.25,
 #   min.segment.length = 0
 #  ) +
  geom_text_repel(
    data = summary,
    aes(x = Weighted_Var, y = Weighted_Output, label = Site),
    size = 3,
    color = "red"
  ) +
  labs(
    x = expression("Variance (" * MW ~ "/" ~ "Hour x $ Millions Spent" * ")"),
    y = expression("Expected (" * MW ~ "/" ~ "Hour x $ Millions Spent" * ")")
  ) +
  theme_minimal(base_size = 12)

#for latek
ggplot() +
  geom_line(data = frontier, aes(x = variance, y = exp_mwhdollar),
            linewidth = 0.7, color = "#08306B") +
  geom_point(data = frontier, aes(x = variance, y = exp_mwhdollar),
             size = 1.5, color = "#08306B", alpha = 0.7) +
  geom_point(data = summary, aes(x = Weighted_Var, y = Weighted_Output),
             color = "red", size = 2) +
  geom_text_repel(
    data = summary,
    aes(x = Weighted_Var, y = Weighted_Output, label = Site),
    size = 3,
    color = "lightpink4"
  ) +
  labs(
    x = expression("Variance (" * MW ~ "/" ~ "Hour x $ Millions Spent" * ")"),
    y = expression("Expected (" * MW ~ "/" ~ "Hour x $ Millions Spent" * ")")
  ) +
  theme_minimal(base_size = 12)

###### ============================================
# Reliability-hour subset
# Currently: Jun-Sep, 8-10pm
# ============================================

idx_rel <- which(
  # Summer late afternoon peak: May–Sep, 3–7pm
  (meta$Month %in% 6:9  & meta$HourOfDay %in% 18:22) #|
    # Winter morning peak: Dec–Mar, 6–9am
   # (meta$Month %in% c(12, 1, 2, 3) & meta$HourOfDay %in% 6:9) |
    # Winter evening peak: Dec–Mar, 6–9pm
    #(meta$Month %in% c(12, 1, 2, 3) & meta$HourOfDay %in% 18:21)
)

X_econ_rel <- X_econ[idx_rel, , drop = FALSE]

mu_rel    <- colMeans(X_econ_rel, na.rm = TRUE)
Sigma_rel <- stats::cov(X_econ_rel, use = "pairwise.complete.obs")

summary_rel <- tibble(
  Site     = sites,
  Weighted_Output = mu_rel,
  Weighted_Var = apply(X_econ_rel, 2, var),
  Mean_output = mu_output
) |>
  dplyr::arrange(dplyr::desc(Weighted_Output))

summary_rel[1, 3] <- 0.030560370

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

#============
# Discrete Gamma Vector
#==============
gamma_vec <- c(0.01, .75, 1.125, 3, 100)

sols_rel <- map(gamma_vec, ~solve_markowitz(mu_rel, Sigma_rel, gamma = .x, solver = "OSQP"))

gname <- function(x) format(x, trim = TRUE, scientific = FALSE)

moments_long_rel <- tibble(
  gamma = gamma_vec,
  `E`   = map_dbl(sols_rel, "exp_output"),
  Var   = map_dbl(sols_rel, "variance")
) |>
  pivot_longer(cols = c(`E`, Var), names_to = "Moment", values_to = "Value")

moments_tbl_rel <- moments_long_rel |>
  mutate(gamma = gname(gamma)) |>
  pivot_wider(names_from = gamma, values_from = Value) |>
  mutate(
    Moment = recode(Moment,
                    `E` = "E[MWh/$]",   
    )
  ) |>
  select(Moment, all_of(gname(gamma_vec)))

weights_long_rel <- map2_dfr(
  sols_rel, gamma_vec,
  ~tibble(gamma = .y, Site = sites, Weight = .x$w)
) %>%
  mutate(
   Weight = ifelse(abs(Weight) < 1e-10, 0, Weight)
  )

weights_tbl_rel <- weights_long_rel |>
  mutate(
    gamma = format(gamma, scientific = FALSE),
    WeightPct = round(100 * Weight, 1)
  ) |>
  select(Site, gamma, WeightPct) |>
  pivot_wider(names_from = gamma, values_from = WeightPct) |>
  arrange(Site)
# =================
# stacked bar chart of weights
# =================
site_order <- c("Casa Grande Solar", "Deming Solar", "Kingman Solar",
                "Encino Wind", "GC Junction Wind", "Medicine Bow Wind", "Silver City Wind")

site_colors <- c(
  # Solar
  "Casa Grande Solar" = "lightpink1",
  "Deming Solar"      = "lightpink3",
  "Kingman Solar"     = "indianred4",
  
  # Wind 
  "Encino Wind"       = "lightsteelblue2",
  "GC Junction Wind"  = "lightskyblue1",
  "Medicine Bow Wind" = "steelblue",
  "Silver City Wind"  = "lightsteelblue4"
)

weights_long_rel_plot <- weights_tbl_rel |>
  pivot_longer(
    cols = -Site,
    names_to = "gamma",
    values_to = "WeightPct"
  ) |>
  mutate(
    gamma = factor(gamma, levels = colnames(weights_tbl_rel)[-1])
  )

weights_long_rel_plot <- weights_long_rel_plot |>
  left_join(site_meta |> select(sitename, tech),
            by = c("Site" = "sitename"))

weights_long_rel_plot <- weights_long_rel_plot |>
  group_by(Site) |>
  filter(sum(WeightPct, na.rm = TRUE) > 0)

weights_long_rel_plot <- weights_long_rel_plot |>
  mutate(
    tech = factor(tech, levels = c("Solar", "Wind"))
  ) |>
  arrange(tech, Site) |>
  mutate(
    Site = factor(Site, levels = unique(Site))
  )

ggplot(weights_long_rel_plot, aes(x = gamma, y = WeightPct, fill = Site)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_manual(
    values = site_colors[site_order],   
    breaks = site_order,                
    drop   = TRUE
  ) +
  geom_text(
    aes(label = ifelse(WeightPct >= 2,
                       sprintf("%.1f%%", WeightPct),
                       "")),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  labs(x = expression("Risk Aversion Parameter, " * gamma), y = "Portfolio Weight (%)", fill = "Site") +
  theme_minimal(base_size = 12)

# ===================
# Plot of discrete gammas
# ===================

moments_plot_rel <- moments_tbl_rel |>
  pivot_longer(
    cols = -Moment,
    names_to = "gamma",
    values_to = "Value"
  ) |>
  pivot_wider(
    names_from = Moment,
    values_from = Value
  ) |>
  rename(
    Variance = `Var`,
    Return = 'E[MWh/$]'
  )

ggplot(moments_plot_rel,
       aes(x = Variance,
           y = Return)) +
  geom_point(size = 3, shape = 8) +
  geom_line(linewidth = .25, alpha = 5) +
  geom_text(
    aes(label = gamma),
    vjust = -1.5,
    size = 3
  ) +
  labs(
    x = "Variance (MWh/Million $)",
    y = "Expected Return (MWh/Million $)"
  ) +
  theme_minimal(base_size = 12)

##############
# Frontier for reliability hours
##############
front_rel <- lapply(gamma_grid, function(g) solve_markowitz(mu_rel, Sigma_rel, gamma = g))

frontier_rel <- tibble(
  gamma    = sapply(front_rel, `[[`, "gamma"),
  exp_output   = sapply(front_rel, `[[`, "exp_output"),
  variance = sapply(front_rel, `[[`, "variance"),
  status   = sapply(front_rel, `[[`, "status")
) |>
  dplyr::filter(status == "optimal") |>
  dplyr::arrange(variance)

ggplot() +
  geom_line(data = frontier_rel, aes(x = variance, y = exp_output),
            linewidth = 0.7, color = "darkblue") +
  geom_point(data = frontier_rel, aes(x = variance, y = exp_output),
             size = 1.5, color = "darkblue", alpha = 0.7) +
  geom_point(data = summary_rel, aes(x = Weighted_Var, y = Weighted_Output),
             color = "red", size = 2) +
  geom_text_repel(
    data = summary_rel,
    aes(x = Weighted_Var, y = Weighted_Output, label = Site),
    size = 3,
    color = "red"
  ) +
  labs(
    x = expression("Variance (" * MW ~ "/" ~ "Hour x $ Millions Spent" * ")"),
    y = expression("Expected (" * MW ~ "/" ~ "Hour x $ Millions Spent" * ")")
  ) +
  theme_minimal(base_size = 12)

