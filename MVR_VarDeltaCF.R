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

# ============================================
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

# Risk parameter for optimization
gamma0     <- 1
gamma_grid <- 10^seq(-10, 3, length.out = 50)

# Nameplate capacities by technology
nameplate_solar_kW <- 75215.82
nameplate_wind_kW  <- 72000     

# ============================================
# Site metadata: sitename, tech, nameplate
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

# Sanity check: no missing nameplates
if (any(is.na(site_meta$nameplate_kW))) {
  print(site_meta)
  stop("Some sites have NA nameplate_kW. Fix 'site_meta' before continuing.")
}

# For fast lookup: nameplate_vec["Wilcox Solar"] -> 75215.82, etc.
nameplate_vec <- setNames(site_meta$nameplate_kW, site_meta$sitename)

# ============================================
# Function to read a single site's CF panel
# ============================================
read_site_long_cf <- function(path, sitename, nameplate_vec) {
  # look up this site's nameplate
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
      Output_kW = dplyr::if_else(Output_kW < 0, 0, Output_kW),
      CF        = Output_kW / np   # <-- site-specific nameplate
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
# Build X matrix (hours x sites) + meta
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
X     <- bx$X
sites <- bx$sites
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
# Helper: compute ramp-based covariance
#   - mean_CF = mean of CF (unchanged)
#   - Sigma_ramp = covariance of hour-to-hour ΔCF within each year
# ============================================
compute_ramp_cov <- function(X, years) {
  # split row indices by year so we don't diff across year boundaries
  idx_split <- split(seq_len(nrow(X)), years)
  
  X_diff_list <- lapply(idx_split, function(idx) {
    X_sub <- X[idx, , drop = FALSE]
    # diff down rows for each column (hour-to-hour CF changes)
    apply(X_sub, 2, diff)
  })
  
  X_diff <- do.call(rbind, X_diff_list)  # (total_hours - #years) x n_sites
  
  list(
    mu    = colMeans(X, na.rm = TRUE),  # same as usual expected CF
    Sigma = stats::cov(X_diff, use = "pairwise.complete.obs")
  )
}

# ============================================
# Full-year: ramp-based risk matrix
# ============================================
ramp_full <- compute_ramp_cov(X, meta$Year)
mu        <- ramp_full$mu
Sigma_ramp <- ramp_full$Sigma   # risk = variance of hourly CF changes

# ============================================
# Markowitz solver
#.max  μᵀw - γ wᵀΣw  s.t. w >= 0, Σw = 1
# ============================================
solve_markowitz <- function(mu, Sigma, gamma = gamma0, lb = NULL, ub = NULL, solver = "OSQP") {
  n <- length(mu)
  w <- CVXR::Variable(n)
  obj <- t(mu) %*% w - gamma * CVXR::quad_form(w, Sigma)
  cons <- list(w >= 0, sum(w) == 1)
  # if (!is.null(lb)) cons <- c(cons, list(w >= lb))
  # if (!is.null(ub)) cons <- c(cons, list(w <= ub))
  
  prob <- CVXR::Problem(CVXR::Maximize(obj), cons)
  res  <- CVXR::solve(prob, solver = solver)
  
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
# Single-γ solution (full year, ramp-based risk)
# ============================================
sol0 <- solve_markowitz(mu, Sigma_ramp, gamma = gamma0)
cat("\n--- Single-γ solution (full year, ramp risk) ---\n")
print(
  tibble(Site = sites, Weight = round(sol0$w, 4)) |>
    dplyr::arrange(dplyr::desc(Weight))
)
cat(sprintf("γ=%g | E[CF]=%.4f | Var(ΔCF)=%.6f\n", sol0$gamma, sol0$exp_CF, sol0$variance))

fronttable <- tibble(Site = sites, Weight = round(sol0$w, 4))

con <- pipe("pbcopy", "w")
write.table(fronttable, con, sep = "\t", col.names = NA)
close(con)

# ============================================
# Efficient frontier (full year, ramp-based risk)
# ============================================
front <- lapply(gamma_grid, function(g) solve_markowitz(mu, Sigma_ramp, gamma = g))

frontier <- tibble(
  gamma    = sapply(front, `[[`, "gamma"),
  exp_CF   = sapply(front, `[[`, "exp_CF"),
  variance = sapply(front, `[[`, "variance")
) |>
  dplyr::arrange(variance)

cat("\n--- Efficient frontier (full year, ramp risk, head) ---\n")
print(head(frontier, 10))

# Base R "baby plot"
plot(frontier$variance, frontier$exp_CF, pch = 20,
     xlab = "Variance of ΔCF (hour-to-hour)",
     ylab = "Expected CF",
     main = "Efficient Frontier (Full Year, Ramp-Based Risk)")
grid()

# Single-site stats for overlay (classical stats, just for context)
summary <- tibble(
  Site     = sites,
  Mean_CF  = mu,                      # in 0–1 units
  Var_CF   = apply(X, 2, var),        # classical variance around mean
  Tech     = site_meta$tech[match(sites, site_meta$sitename)]
) |>
  dplyr::arrange(dplyr::desc(Mean_CF))

con <- pipe("pbcopy", "w")
write.table(summary, con, sep = "\t", col.names = NA)
close(con)

# --- Single-site ramp variance ---
site_ramp_stats <- tibble(
  Site      = sites,
  Mean_CF   = mu,
  Var_Delta = diag(Sigma_ramp)   # diagonal of ΣΔ = variance of ΔCF for each site
)

# --- Plot efficient frontier + single-site points ---
ggplot() +
  # frontier curve
  geom_line(
    data = frontier,
    aes(x = variance, y = exp_CF),
    linewidth = 0.8,
    color = "darkblue"
  ) +
  geom_point(
    data = frontier,
    aes(x = variance, y = exp_CF),
    size = 1.5,
    color = "darkblue",
    alpha = 0.7
  ) +
  
  # single-site portfolios
  geom_point(
    data = site_ramp_stats,
    aes(x = Var_Delta, y = Mean_CF, color = "Single-Site"),
    size = 2.8
  ) +
  
  geom_text_repel(
    data = site_ramp_stats,
    aes(x = Var_Delta, y = Mean_CF, label = Site),
    size = 3.4,
    color = "black"
  ) +
  
  labs(
    x = "Variance of Hour-to-Hour ΔCF",
    y = "Expected Capacity Factor",
    title = "Efficient Frontier (Ramp-Based Risk)",
    color = ""
  ) +
  scale_color_manual(values = c("Single-Site" = "red")) +
  theme_minimal(base_size = 12)

# ============================================
# July–August subset (ramp-based risk)
# ============================================
idx_JA  <- which(meta$Month %in% c(7, 8))
X_JA    <- X[idx_JA, , drop = FALSE]
meta_JA <- meta[idx_JA, ]

ramp_JA <- compute_ramp_cov(X_JA, meta_JA$Year)
mu_JA   <- ramp_JA$mu
Sigma_JA_ramp <- ramp_JA$Sigma

# Single-γ solution for Jul–Aug (ramp)
sol_JA <- solve_markowitz(mu_JA, Sigma_JA_ramp, gamma = gamma0)

cat("\n=== July–August: Single-γ solution (ramp risk) ===\n")
print(
  tibble(Site = sites, Weight = round(sol_JA$w, 4)) |>
    dplyr::arrange(dplyr::desc(Weight))
)
cat(sprintf("status=%s | γ=%g | E[CF]_JA=%.4f | Var_JA(ΔCF)=%.6f\n",
            sol_JA$status, sol_JA$gamma, sol_JA$exp_CF, sol_JA$variance))

JA <- tibble(Site = sites, Weight = round(sol_JA$w, 4))

con <- pipe("pbcopy", "w")
write.table(JA, con, sep = "\t", col.names = NA)
close(con)

# Frontier for Jul–Aug (ramp risk)
front_JA <- lapply(gamma_grid, function(g) solve_markowitz(mu_JA, Sigma_JA_ramp, gamma = g))

frontier_JA <- tibble(
  gamma    = sapply(front_JA, `[[`, "gamma"),
  exp_CF   = sapply(front_JA, `[[`, "exp_CF"),
  variance = sapply(front_JA, `[[`, "variance"),
  status   = sapply(front_JA, `[[`, "status")
) |>
  dplyr::filter(status == "optimal") |>
  dplyr::arrange(variance)

cat("\n=== July–August: Efficient frontier (ramp risk, head) ===\n")
print(head(frontier_JA, 8))

plot(frontier_JA$variance, frontier_JA$exp_CF, pch = 16,
     xlab = "Variance of ΔCF (Jul–Aug)",
     ylab = "Expected CF (Jul–Aug)",
     main = "Efficient Frontier — July–August (Ramp-Based Risk)")
grid()

# Jul–Aug site stats (for context)
site_stats_JA <- tibble(
  Site     = sites,
  Mean_CF  = mu_JA,
  Var_CF   = apply(X_JA, 2, var)
) |>
  dplyr::arrange(dplyr::desc(Mean_CF))

cat("\n=== July–August: Site stats (classical variance) ===\n")
print(site_stats_JA)

# Optional: correlation of ramps (could be interesting later)
# (differences computed inside compute_ramp_cov; if you want them separately we can add that)
