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
gamma0     <- 2.5
gamma_grid <- 10^seq(-10, 3, length.out = 50)

# ============================================
# Nameplate capacities
# ============================================
nameplate_solar_kW <- 100000
nameplate_wind_kW  <- 100000

# ============================================
# Site metadata: name, solar or wind, nameplate capacity
# ============================================
site_names <- basename(site_files) |> tools::file_path_sans_ext()

site_names[site_names == "Wilcox Solar"] <- "Willcox Solar"

site_names

site_meta <- tibble(
  file     = site_files,
  sitename = site_names
) %>%
  mutate(
    tech = case_when(
      str_detect(sitename, regex("wind",  ignore_case = TRUE)) ~ "Wind",
      str_detect(sitename, regex("solar", ignore_case = TRUE)) ~ "Solar",
      TRUE ~ "Unknown"
    ),
    nameplate_kW = case_when(
      tech == "Solar" ~ nameplate_solar_kW,
      tech == "Wind"  ~ nameplate_wind_kW,
      TRUE ~ NA_real_
    )
  )

nameplate_vec <- setNames(site_meta$nameplate_kW, site_meta$sitename)

# ============================================
# CF Panel
# ============================================
read_site_long_cf <- function(path, sitename, nameplate_vec) {
  np <- nameplate_vec[[sitename]]
  if (is.null(np) || is.na(np)) {
    stop(sprintf("No nameplate_kW found for site '%s'", sitename))
  }
  
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
      Output_kW = Output_kW,
      CF        = Output_kW / np   
    ) |>
    dplyr::select(Year, Hour, Site, CF)
  
  df_long
}

# ============================================
# Build panel of CFs for all sites
# ============================================
panel_cf <- purrr::map2_dfr(
  site_files, site_names,
  ~ read_site_long_cf(.x, .y, nameplate_vec = nameplate_vec)
)

# ============================================
# Build X matrix (hours x sites) amd the time data
# ============================================
build_X <- function(panel_cf) {
  W <- panel_cf |>
    dplyr::arrange(Year, Hour, Site) |>
    tidyr::pivot_wider(names_from = Site, values_from = CF) |>
    dplyr::arrange(Year, Hour)
  
  site_cols <- setdiff(names(W), c("Year", "Hour"))
  
  list(
    X     = as.matrix(W[ , site_cols]),
    sites = site_cols,
    meta  = W[, c("Year", "Hour")]
  )
}

bx    <- build_X(panel_cf)
#Hour x Site Matrix
X     <- bx$X    
#A list of my sites (for labeling and joining metadata)
sites <- bx$sites
#time metadata
meta  <- bx$meta

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
mu    <- colMeans(X)
Sigma <- stats::cov(X)

# ============================================
#Site Summary Table
# ============================================
summary <- tibble(
  Site     = sites,
  Mean_CF  = mu,                    
  Variance = apply(X, 2, var)
) |>
  dplyr::arrange(dplyr::desc(Mean_CF))

summary <- summary %>%
  mutate(Mean_CF = Mean_CF * 100)
  

con <- pipe("pbcopy", "w")
write.table(summary, con, sep = "\t", col.names = NA)
close(con)

# ============================================
#Site x Seasonal Summary Table
# ============================================
m <- meta$Month

season_vec <- dplyr::case_when(
  m %in% c(12, 1, 2) ~ "Winter",
  m %in% c(3, 4, 5)  ~ "Spring",
  m %in% c(6, 7, 8)  ~ "Summer",
  m %in% c(9,10,11)  ~ "Fall",
  TRUE               ~ NA_character_
)

summary_seasonal <- as_tibble(X) |>
  mutate(Season = season_vec) |>
  pivot_longer(
    cols = -Season,
    names_to  = "Site",
    values_to = "CF"
  ) |>
  group_by(Season, Site) |>
  summarise(
    Mean_CF  = mean(CF, na.rm = TRUE),
    Variance = var(CF,  na.rm = TRUE),
    .groups  = "drop"
  ) |>
  arrange(Season, desc(Mean_CF))

summary_seasonal <- as_tibble(X) |>
  mutate(Season = season_vec) |>
  pivot_longer(
    cols = -Season,
    names_to  = "Site",
    values_to = "CF"
  ) |>
  group_by(Season, Site) |>
  summarise(
    Mean_CF  = mean(CF, na.rm = TRUE),
    Variance = var(CF,  na.rm = TRUE),
    .groups  = "drop"
  )

season_table <- summary_seasonal |>
  mutate(Cell = sprintf("%.3f", Mean_CF)) |>
  select(Season, Site, Cell) |>
  pivot_wider(names_from = Site, values_from = Cell) |>
  arrange(Season)

# ============================================
# Correlation matrices
# ============================================
cor_X <- cor(X, use = "pairwise.complete.obs")
print(round(cor_X, 3))

WindSolar <- X[, c("Medicine Bow Wind", "Encino Wind", "Silver City Wind", 
                   "GC Junction Wind", "Kingman Solar","GC Junction Solar",
                   "Casa Grande Solar", "Willcox Solar", "St Johns Solar", 
                   "Deming Solar")]

#cor_WindSolar <- cor(WindSolar, use = "pairwise.complete.obs")
#round(cor_WindSolar, 3)

#con <- pipe("pbcopy", "w")
#write.table(round(cor_X, 3), con, sep = "\t", col.names = NA)
#close(con)

#wind_only <- X[, c("Encino Wind", "Medicine Bow Wind", "Silver City Wind", "GC Junction Wind")]
#cor_wind <- cor(wind_only, use = "pairwise.complete.obs")
#round(cor_wind, 3)

#GCJ_only <- X[, c("GC Junction Solar", "GC Junction Wind")]
#cor_GCJ <- cor(GCJ_only, use = "pairwise.complete.obs")
#round(cor_GCJ, 3)

WindSolar <- X[, c("Medicine Bow Wind", "Encino Wind", "Willcox Solar", "Kingman Solar")]
cor_WindSolar <- cor(WindSolar, use = "pairwise.complete.obs")
round(cor_WindSolar, 3)

con <- pipe("pbcopy", "w")
write.table(round(cor_WindSolar, 3), con, sep = "\t", col.names = NA)
close(con)

# ============================================
# Markowitz solver
# max  μᵀw - γ wᵀΣw  s.t. w >= 0, Σw = 1
# ============================================
solve_markowitz <- function(mu, Sigma, gamma = gamma0, lb = NULL, ub = NULL, solver = "OSQP") {
  n <- length(mu)
  w <- CVXR::Variable(n)
  obj <- t(mu) %*% w - gamma * CVXR::quad_form(w, Sigma)
  cons <- list(w >= 0, sum(w) == 1)
  
  prob <- CVXR::Problem(CVXR::Maximize(obj), cons)
  res  <- CVXR::solve(prob, solver = solver)
  #im telling CVXR To maximize my objective function (obj) 
  #subject to consraints (cons)
  #res = solver executes optimization
  
  wv <- as.numeric(res$getValue(w))
  list(
    status   = res$status,
    w        = wv,
    exp_CF   = sum(mu * wv),
    variance = as.numeric(t(wv) %*% Sigma %*% wv),
    gamma    = gamma
  )
}

# ============================================
# Single-γ solution
# ============================================
gamma0     <- 3

sol0 <- solve_markowitz(mu, Sigma, gamma = gamma0)
cat("\n--- Single-γ solution ---\n")
print(
  tibble(Site = sites, Weight = round(sol0$w, 4)) |>
    dplyr::arrange(dplyr::desc(Weight))
)
 cat(sprintf("γ=%g | E[CF]=%.4f | Var=%.4f\n", sol0$gamma, sol0$exp_CF, sol0$variance))

fronttable <- tibble(Site = sites, Weight = round(sol0$w, 4))

con <- pipe("pbcopy", "w")
write.table(fronttable, con, sep = "\t", col.names = NA)
close(con)

#============
# Gamma Vector
#==============
gamma_vec <- c(0.01, .75, 1.5, 3, 4.5, 100)

sols <- map(gamma_vec, ~solve_markowitz(mu, Sigma, gamma = .x, solver = "OSQP"))

gname <- function(x) format(x, trim = TRUE, scientific = FALSE)

moments_long <- tibble(
  gamma = gamma_vec,
  `E`   = map_dbl(sols, "exp_CF"),
  Var   = map_dbl(sols, "variance")
) |>
  pivot_longer(cols = c(`E`, Var), names_to = "Moment", values_to = "Value")

moments_tbl <- moments_long |>
  mutate(gamma = gname(gamma)) |>
  pivot_wider(names_from = gamma, values_from = Value) |>
  mutate(
    Moment = recode(Moment,
                    `E` = "Expected Capacity Factor (%)",   
    )
  ) |>
  select(Moment, all_of(gname(gamma_vec)))

moments_tbl <- moments_tbl %>%
  rename("100" = "100.00")

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

site_order <- c("Casa Grande Solar", "Deming Solar", "Kingman Solar",
                "Encino Wind", "GC Junction Wind", "Medicine Bow Wind", "Silver City Wind")

site_colors <- c(
  # Solar
  "Casa Grande Solar" = "lightpink1",
  "Deming Solar"      = "lightpink3",
  "Kingman Solar"     = "indianred3",
  
  # Wind 
  "Encino Wind"       = "lightsteelblue2",
  "GC Junction Wind"  = "lightskyblue1",
  "Medicine Bow Wind" = "steelblue",
  "Silver City Wind"  = "lightsteelblue4"
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
    ECF = 'Expected Capacity Factor (%)'
  ) |>
  mutate(
    gamma = factor(gamma, levels = colnames(moments_tbl)[-1])
  )

ggplot(moments_plot,
       aes(x = Variance,
           y = ECF)) +
  geom_point(size = 3, shape = 8) +
  geom_line(linewidth = .25, alpha = 5) +
  geom_text(
    aes(label = gamma),
    vjust = -0.8,
    size = 3
  ) +
  labs(
    x = "Variance of Capacity Factor",
    y = "Expected Capacity Factor (%)"
  ) +
  theme_minimal(base_size = 12)

# ============================================
# Efficient frontier along gamma grid
# ============================================
front <- lapply(gamma_grid, function(g) solve_markowitz(mu, Sigma, gamma = g))

frontier <- tibble(
  gamma    = sapply(front, `[[`, "gamma"),
  exp_CF   = sapply(front, `[[`, "exp_CF"),
  variance = sapply(front, `[[`, "variance")
) |>
  dplyr::arrange(variance)

frontier <- frontier %>%
  mutate(exp_CF = exp_CF * 100)

cat("\n--- Efficient frontier (full year, head) ---\n")
print(head(frontier, 20))

#manual!!! change as needed
#gamma_points <- tibble::tibble(
 # gamma     = c(0.01, 0.75, 1.5, 3, 100),
 #  variance  = c(0.1368, 0.0844, 0.0474, 0.0373, 0.0243),
 #  exp_CF    = c(44.6, 43.8, 40.8, 38.9, 31.6))

# ggplot frontier + single-site points
ggplot() +
  geom_line(data = frontier, aes(x = variance, y = exp_CF),
            linewidth = 0.7, color = "#08306B") +
  geom_point(data = frontier, aes(x = variance, y = exp_CF),
             size = 1.5, color = "#08306B", alpha = 0.7) +
  geom_point(data = summary, aes(x = Variance, y = Mean_CF),
             color = "red3", size = 2) +
  geom_text_repel(nudge_y = 1.1, nudge_x = -0.003,
    data = moments_plot,
    aes(x = Variance, y = ECF*100, label = gamma),
    size = 3) +
  geom_text_repel(
    data = summary,
    aes(x = Variance, y = Mean_CF, label = Site),
    size = 3,
    color = "red"
  ) +
  labs(
    x = "Variance of Capacity Factor",
    y = "Expected Capacity Factor (%)",  # change to (%) if you multiply Mean_CF by 100
  ) +
  theme_minimal(base_size = 12)

# ============================================
# July–August subset
# ============================================
idx_JA <- which(meta$Month %in% c(7, 8))
X_JA   <- X[idx_JA, , drop = FALSE]

mu_JA    <- colMeans(X_JA)
Sigma_JA <- stats::cov(X_JA, use = "pairwise.complete.obs")

# Single-γ solution for Jul–Aug
sol_JA <- solve_markowitz(mu_JA, Sigma_JA, gamma = gamma0)

cat("\n=== July–August: Single-γ solution ===\n")
print(
  tibble(Site = sites, Weight = round(sol_JA$w, 4)) |>
    dplyr::arrange(dplyr::desc(Weight))
)
cat(sprintf("status=%s | γ=%g | E[CF]_JA=%.4f | Var_JA=%.4f\n",
            sol_JA$status, sol_JA$gamma, sol_JA$exp_CF, sol_JA$variance))

JA <- tibble(Site = sites, Weight = round(sol_JA$w, 4))

con <- pipe("pbcopy", "w")
write.table(JA, con, sep = "\t", col.names = NA)
close(con)

# Frontier for Jul–Aug
front_JA <- lapply(gamma_grid, function(g) solve_markowitz(mu_JA, Sigma_JA, gamma = g))

frontier_JA <- tibble(
  gamma    = sapply(front_JA, `[[`, "gamma"),
  exp_CF   = sapply(front_JA, `[[`, "exp_CF"),
  variance = sapply(front_JA, `[[`, "variance"),
  status   = sapply(front_JA, `[[`, "status")
) |>
  dplyr::filter(status == "optimal") |>
  dplyr::arrange(variance)

cat("\n=== July–August: Efficient frontier (head) ===\n")
print(head(frontier_JA, 8))

plot(frontier_JA$variance, frontier_JA$exp_CF, pch = 16,
     xlab = "Variance of CF (Jul–Aug)", ylab = "Expected CF (Jul–Aug)",
     main = "Efficient Frontier — July–August")
grid()

# Jul–Aug site stats
site_stats_JA <- tibble(
  Site     = sites,
  Mean_CF  = mu_JA,
  Variance = apply(X_JA, 2, var)
) |>
  dplyr::arrange(dplyr::desc(Mean_CF))

cat("\n=== July–August: Site stats ===\n")
print(site_stats_JA)

# Jul–Aug correlation matrix
cor_X_JA <- cor(X_JA, use = "pairwise.complete.obs")
print(round(cor_X_JA, 3))

con <- pipe("pbcopy", "w")
write.table(round(cor_X_JA, 3), con, sep = "\t", col.names = NA)
close(con)


###### ============================================
# Reliability-hour subset (custom windows)
# Comment in and out as you please
# ============================================

idx_rel <- which(
  # Summer late afternoon peak: May–Sep, 3–7pm
  (meta$Month %in% 6:9  & meta$HourOfDay %in% 16:21) #|
    # Winter morning peak: Dec–Mar, 6–9am
   # (meta$Month %in% c(12, 1, 2, 3) & meta$HourOfDay %in% 6:9) |
    # Winter evening peak: Dec–Mar, 6–9pm
   # (meta$Month %in% c(12, 1, 2, 3) & meta$HourOfDay %in% 18:21)
)

#generates the x matrix with the subsetted hours
X_rel <- X[idx_rel, , drop = FALSE]

# mean + variance
mu_rel    <- colMeans(X_rel, na.rm = TRUE)
Sigma_rel <- stats::cov(X_rel, use = "pairwise.complete.obs")

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

# Reliability-hour site stats + correlation
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


#######
#Daily Aggregates
#######

# ============================================
#Daily aggregation (daily mean CF)
# ============================================

# Day number within each Year (1..365) already exists as meta$DayOfYear.
# Create a unique day ID across the whole sample:
meta <- meta %>%
  mutate(DayID = paste0(Year, "-", sprintf("%03d", DayOfYear)))

# Aggregate hourly CF to daily mean for each site
X_day <- apply(X, 2, function(v) tapply(v, meta$DayID, mean, na.rm = TRUE))
colnames(X_day) <- sites
rownames(X_day) <- sort(unique(meta$DayID))  # matches tapply ordering (usually)

# Build daily metadata aligned to rows of X_day
meta_day <- meta %>%
  distinct(DayID, Year, DayOfYear, Month) %>%
  arrange(DayID)

# ============
# Mean & covariance
# ============
mu_day    <- colMeans(X_day, na.rm = TRUE)
Sigma_day <- stats::cov(X_day, use = "pairwise.complete.obs")

# Single-γ daily solution
sol0_day <- solve_markowitz(mu_day, Sigma_day, gamma = gamma0)

cat("\n--- Single-γ DAILY solution ---\n")
print(
  tibble(Site = sites, Weight = round(sol0_day$w, 4)) |>
    dplyr::arrange(dplyr::desc(Weight))
)
cat(sprintf("γ=%g | E[CF]_day=%.4f | Var_day=%.6f\n",
            sol0_day$gamma, sol0_day$exp_CF, sol0_day$variance))

# Efficient frontier (daily)
front_day <- lapply(gamma_grid, function(g) solve_markowitz(mu_day, Sigma_day, gamma = g))

frontier_day <- tibble(
  gamma    = sapply(front_day, `[[`, "gamma"),
  exp_CF   = sapply(front_day, `[[`, "exp_CF"),
  variance = sapply(front_day, `[[`, "variance"),
  status   = sapply(front_day, `[[`, "status")
) |>
  dplyr::filter(status == "optimal") |>
  dplyr::arrange(variance)

frontier_day <- frontier_day %>%
  mutate(exp_CF = exp_CF * 100)

cat("\n--- Efficient frontier (DAILY, head) ---\n")
print(head(frontier_day, 20))

# Add labels to your already-computed hourly frontier
frontier_hour <- frontier %>% mutate(Granularity = "Hourly")
frontier_day2  <- frontier_day %>% select(gamma, exp_CF, variance) %>% mutate(Granularity = "Daily")

front_compare <- bind_rows(frontier_hour, frontier_day2)

ggplot(front_compare, aes(x = variance, y = exp_CF, color = Granularity)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.6) +
  scale_color_manual(values = c("Daily" = "steelblue2", "Hourly" = "#08306B")) +
  labs(
    x = "Variance of Capacity Factor",
    y = "Expected Capacity Factor (%)"
  ) +
  theme_minimal(base_size = 12)

# Combine weights side-by-side
w_compare <- tibble(
  Site   = sites,
  Hourly = sol0$w,
  Daily  = sol0_day$w,
  Diff   = sol0_day$w - sol0$w
) |>
  mutate(across(c(Hourly, Daily, Diff), ~round(.x, 4))) |>
  arrange(desc(abs(Diff)))

cat("\n--- Weight comparison at gamma0 ---\n")
print(w_compare)

# Compare objective components at gamma0
metrics_compare <- tibble(
  freq     = c("Hourly", "Daily mean"),
  exp_CF   = c(sol0$exp_CF, sol0_day$exp_CF),
  variance = c(sol0$variance, sol0_day$variance)
)

metrics_wide <- metrics_compare %>%
  mutate(freq = recode(freq, "Daily mean" = "Daily_mean")) %>%  # optional: avoid spaces
  pivot_longer(cols = c(exp_CF, variance), names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = freq, values_from = value) %>%
  mutate(
    pct_change = (Hourly - Daily_mean) / Hourly * 100
  ) %>%
  mutate(
    metric = recode(metric,
                    exp_CF   = "Expected CF",
                    variance = "Variance")
  )
