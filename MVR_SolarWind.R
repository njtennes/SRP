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
gamma0     <- 2
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

con <- pipe("pbcopy", "w")
write.table(summary, con, sep = "\t", col.names = NA)
close(con)

# ============================================
# Correlation matrices
# ============================================
cor_X <- cor(X, use = "pairwise.complete.obs")
print(round(cor_X, 3))

WindSolar <- X[, c("Medicine Bow Wind", "Encino Wind", "Silver City Wind", 
                   "GC Junction Wind", "Kingman Solar","GC Junction Solar",
                   "Casa Grande Solar", "Wilcox Solar", "St Johns Solar", 
                   "Deming Solar")]

cor_WindSolar <- cor(WindSolar, use = "pairwise.complete.obs")
round(cor_WindSolar, 3)

con <- pipe("pbcopy", "w")
write.table(round(cor_X, 3), con, sep = "\t", col.names = NA)
close(con)

wind_only <- X[, c("Encino Wind", "Medicine Bow Wind", "Silver City Wind", "GC Junction Wind")]
cor_wind <- cor(wind_only, use = "pairwise.complete.obs")
round(cor_wind, 3)

GCJ_only <- X[, c("GC Junction Solar", "GC Junction Wind")]
cor_GCJ <- cor(GCJ_only, use = "pairwise.complete.obs")
round(cor_GCJ, 3)

WindSolar <- X[, c("Medicine Bow Wind", "Encino Wind", "Wilcox Solar", "Kingman Solar")]
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

cat("\n--- Efficient frontier (full year, head) ---\n")
print(head(frontier, 20))

# ggplot frontier + single-site points
ggplot() +
  geom_line(data = frontier, aes(x = variance, y = exp_CF),
            linewidth = 0.7, color = "darkblue") +
  geom_point(data = frontier, aes(x = variance, y = exp_CF),
             size = 1.5, color = "darkblue", alpha = 0.7) +
  geom_point(data = summary, aes(x = Variance, y = Mean_CF),
             color = "red", size = 2) +
  geom_text_repel(
    data = summary,
    aes(x = Variance, y = Mean_CF, label = Site),
    size = 3,
    color = "red"
  ) +
  labs(
    x = "Variance of Capacity Factor",
    y = "Expected Capacity Factor",  # change to (%) if you multiply Mean_CF by 100
    title = "Efficient Frontier with Single-Site Portfolios"
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
# May–Sep 3–7pm, Dec–Mar 6–9am and 6–9pm
# ============================================

idx_rel <- which(
  # Summer late afternoon peak: May–Sep, 3–7pm
  (meta$Month %in% 5:9  & meta$HourOfDay %in% 15:19) |
    # Winter morning peak: Dec–Mar, 6–9am
    (meta$Month %in% c(12, 1, 2, 3) & meta$HourOfDay %in% 6:9) |
    # Winter evening peak: Dec–Mar, 6–9pm
    (meta$Month %in% c(12, 1, 2, 3) & meta$HourOfDay %in% 18:21)
)

X_rel <- X[idx_rel, , drop = FALSE]

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
