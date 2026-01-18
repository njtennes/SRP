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
gamma0     <- .3
gamma_grid <- 10^seq(-10, 3, length.out = 50)

site_names <- basename(site_files) |> tools::file_path_sans_ext()

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

cost_per_MWh <- c(
  "Medicine Bow Wind" = 32,
  "Encino Wind"       = 40,
  "Silver City Wind"  = 40,
  "GC Junction Wind"  = 38,
  "Kingman Solar"     = 28,
  "GC Junction Solar" = 30,
  "Casa Grande Solar" = 29,
  "Wilcox Solar"      = 27,
  "St Johns Solar"    = 31,
  "Deming Solar"      = 33
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
  Mean_dollar = mu,
  Variance = apply(X_MW, 2, var),
  Mean_output = mu_output
) |>
  dplyr::arrange(dplyr::desc(Mean_dollar))

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
  exp_mwhdollar  = sapply(front, `[[`, "exp_output"),
  variance = sapply(front, `[[`, "variance")
) |>
  dplyr::arrange(variance)

cat("\n--- Efficient frontier (full year, head) ---\n")
print(head(frontier, 20))

# ggplot frontier + single-site points
ggplot() +
  geom_line(data = frontier, aes(x = variance, y = exp_mwhdollar),
            linewidth = 0.7, color = "darkblue") +
  geom_point(data = frontier, aes(x = variance, y = exp_mwhdollar),
             size = 1.5, color = "darkblue", alpha = 0.7) +
  geom_point(data = summary, aes(x = Variance, y = Mean_dollar),
             color = "red", size = 2) +
  geom_text_repel(
    data = summary,
    aes(x = Variance, y = Mean_dollar, label = Site),
    size = 3,
    color = "red"
  ) +
  labs(
    x = "Variance MWHr/$ spent",
    y = "Expected MWHr/$ spent", 
    title = "Efficient Frontier with Single-Site Portfolios"
  ) +
  theme_minimal(base_size = 12)

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
