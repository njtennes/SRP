library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(CVXR)

#NSRDB SITE FILES
site_files <- c(
  "/Users/nicktennes/Desktop/NSRDB Solar Output/CasaGrande.csv",
  "/Users/nicktennes/Desktop/NSRDB Solar Output/CGJunction.csv",
  "/Users/nicktennes/Desktop/NSRDB Solar Output/Deming.csv",
  "/Users/nicktennes/Desktop/NSRDB Solar Output/Kingman.csv",
  "/Users/nicktennes/Desktop/NSRDB Solar Output/StJohns.csv",
  "/Users/nicktennes/Desktop/NSRDB Solar Output/Wilcox.csv"
)

#ERA5 SITE FILES
site_files <- c(
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/Casa Grande Solar.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/Deming Solar.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/GC Junction Solar.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/Kingman Solar.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/St Johns Solar.csv',
  '/Users/nicktennes/Documents/SAM ERA5 Output 8760/Wilcox Solar.csv'
)


#setting the nameplate capacity 
nameplate_kW <- 75215.82

#choosing a risk parameter, gamma, to identify weight
gamma0       <- 1

#define a sequence of gammas (to find the efficient frontier)
gamma_grid   <- 10^seq(-10, 3, length.out = 50)

#cool tool that sites the sitename based on the name of my solar files
site_names <- basename(site_files) |> tools::file_path_sans_ext()

#create a masterfile of outputs, floor output to 0 (at night), turn to CF
read_site_long_cf <- function(path, sitename, nameplate_kW) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  stopifnot("Hour" %in% names(df))
  df$Hour <- as.integer(df$Hour)
  
  df_long <- df |>
    tidyr::pivot_longer(-Hour, names_to = "Year", values_to = "Output_kW") |>
    dplyr::mutate(
      Year = as.integer(Year),
      Site = sitename,
      Output_kW = dplyr::if_else(Output_kW < 0, 0, Output_kW),
      CF = Output_kW / nameplate_kW
    ) |>
    dplyr::select(Year, Hour, Site, CF)
  
  df_long
}

panel_cf <- purrr::map2_dfr(
  site_files, site_names,
  ~ read_site_long_cf(.x, .y, nameplate_kW)
)

#create the full capacity factor matrix
build_X <- function(panel_cf) {
  W <- panel_cf |>
    dplyr::arrange(Year, Hour, Site) |>
    tidyr::pivot_wider(names_from = Site, values_from = CF) |>
    dplyr::arrange(Year, Hour)
  site_cols <- setdiff(names(W), c("Year", "Hour"))
  list(
    X = as.matrix(W[ , site_cols]),
    sites = site_cols,
    meta = W[, c("Year", "Hour")]
  )
}

bx    <- build_X(panel_cf)
X     <- bx$X
sites <- bx$sites
meta  <- bx$meta

# Derive HourOfDay and DayOfYear from Hour
meta <- meta %>%
  mutate(
    HourOfDay = ((Hour - 1) %% 24),        # 0–23 convention
    DayOfYear = ((Hour - 1) %/% 24) + 1    # 1–365
  )

#estimate month
meta <- meta %>%
  mutate(
    Month = case_when(
      DayOfYear <= 31   ~ 1,
      DayOfYear <= 31 + 28   ~ 2,
      DayOfYear <= 31 + 28 + 31   ~ 3,
      DayOfYear <= 31 + 28 + 31 + 30   ~ 4,
      DayOfYear <= 31 + 28 + 31 + 30 + 31   ~ 5,
      DayOfYear <= 31 + 28 + 31 + 30 + 31 + 30   ~ 6,
      DayOfYear <= 31 + 28 + 31 + 30 + 31 + 30 + 31   ~ 7,
      DayOfYear <= 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31   ~ 8,
      DayOfYear <= 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30   ~ 9,
      DayOfYear <= 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31   ~ 10,
      DayOfYear <= 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30   ~ 11,
      TRUE ~ 12
    )
  )

mu    <- colMeans(X)
Sigma <- stats::cov(X, use = "pairwise.complete.obs")

summary <- tibble(
  Site     = sites,
  Mean_CF  = mu,                      
  Variance = apply(X, 2, var)
) |> arrange(desc(Mean_CF))

con <- pipe("pbcopy", "w")
write.table(summary, con, sep = "\t", col.names = NA)
close(con)

solve_markowitz <- function(mu, Sigma, gamma = gamma0, lb = NULL, ub = NULL, solver = "OSQP") {
  n <- length(mu)
  w <- CVXR::Variable(n)
  obj <- t(mu) %*% w - gamma * CVXR::quad_form(w, Sigma)
  cons <- list(w >= 0, sum(w) == 1)
  #if (!is.null(lb)) cons <- c(cons, list(w >= lb))
  #if (!is.null(ub)) cons <- c(cons, list(w <= ub))
  prob <- Problem(Maximize(obj), cons)
  res  <- solve(prob, solver = solver)
  wv   <- as.numeric(res$getValue(w))
  list(
    status   = res$status,
    w        = wv,
    exp_CF   = sum(mu * wv),
    variance = as.numeric(t(wv) %*% Sigma %*% wv),
    gamma    = gamma
  )
}

sol0 <- solve_markowitz(mu, Sigma, gamma = gamma0)
cat("\n--- Single-γ solution ---\n")
print(tibble(Site = sites, Weight = round(sol0$w, 4)) |> dplyr::arrange(dplyr::desc(Weight)))
cat(sprintf("γ=%g | E[CF]=%.4f | Var=%.4f\n", sol0$gamma, sol0$exp_CF, sol0$variance))

fronttable <- tibble(Site = sites, Weight = round(sol0$w, 4))

con <- pipe("pbcopy", "w")
write.table(fronttable, con, sep = "\t", col.names = NA)
close(con)

front <- lapply(gamma_grid, function(g) solve_markowitz(mu, Sigma, gamma = g))
frontier <- tibble(
  gamma    = sapply(front, `[[`, "gamma"),
  exp_CF   = sapply(front, `[[`, "exp_CF"),
  variance = sapply(front, `[[`, "variance")
) |> dplyr::arrange(variance)

cat("\n--- Efficient frontier (head) ---\n")
print(head(frontier, 10))


summary <- tibble(
  Site     = sites,
  Mean_CF  = mu,                      
  Variance = apply(X, 2, var)
) |> arrange(desc(Mean_CF))

con <- pipe("pbcopy", "w")
write.table(summary, con, sep = "\t", col.names = NA)
close(con)

#ggplot frontier
library(ggrepel)
library(ggplot2)
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
    y = "Expected Capacity Factor (%)",
    title = "Efficient Frontier with Single-Site Portfolios"
  ) +
  theme_minimal(base_size = 12)



cor_X <- cor(X, use = "pairwise.complete.obs")
print(round(cor_X, 3))

con <- pipe("pbcopy", "w")
write.table(round(cor_X, 3), con, sep = "\t", col.names = NA)
close(con)

#===========
# subset rows to July/Aug
idx_JA  <- which(meta$Month %in% c(7, 8))

X_JA    <- X[idx_JA, , drop = FALSE]

# July/Aug mean variance
mu_JA    <- colMeans(X_JA)
Sigma_JA <- stats::cov(X_JA, use = "pairwise.complete.obs")

# single gamma July/Aug
sol_JA <- solve_markowitz(mu_JA, Sigma_JA, gamma = gamma0)

cat("\n=== July–August: Single-γ solution ===\n")
print(tibble(Site = sites, Weight = round(sol_JA$w, 4)) |>
        dplyr::arrange(dplyr::desc(Weight)))
cat(sprintf("status=%s | γ=%g | E[CF]_JA=%.4f | Var_JA=%.4f\n",
            sol_JA$status, sol_JA$gamma, sol_JA$exp_CF, sol_JA$variance))

JA <- tibble(Site = sites, Weight = round(sol_JA$w, 4))

con <- pipe("pbcopy", "w")
write.table(JA, con, sep = "\t", col.names = NA)
close(con)

#frontier
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

frontier
plot(frontier_JA$variance, frontier_JA$exp_CF, pch = 16,
     xlab = "Variance of CF (Jul–Aug)", ylab = "Expected CF (Jul–Aug)",
     main = "Efficient Frontier — July–August")
grid()

# summary states
site_stats_JA <- tibble(
  Site     = sites,
  Mean_CF  = mu_JA,
  Variance = apply(X_JA, 2, var)
) |>
  dplyr::arrange(dplyr::desc(Mean_CF))

cat("\n=== July–August: Site stats ===\n")
print(site_stats_JA)

cor_X_JA <- cor(X_JA, use = "pairwise.complete.obs")
print(round(cor_X_JA, 3))

con <- pipe("pbcopy", "w")
write.table(round(cor_X_JA, 3), con, sep = "\t", col.names = NA)
close(con)
