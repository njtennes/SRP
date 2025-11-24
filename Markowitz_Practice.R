# =========================
# Libraries
# =========================
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(CVXR)

# =========================
# USER INPUTS
# =========================
site_files <- c(
  "/Users/nicktennes/Desktop/Solar/CasaGrande.csv",
  "/Users/nicktennes/Desktop/Solar/CGJunction.csv",
  "/Users/nicktennes/Desktop/Solar/Deming.csv",
  "/Users/nicktennes/Desktop/Solar/Kingman.csv",
  "/Users/nicktennes/Desktop/Solar/StJohns.csv",
  "/Users/nicktennes/Desktop/Solar/Wilcox.csv"
)

# Use a single nameplate for all solar sites (kW)
nameplate_kW <- 75000

weight_mode  <- "peak_emphasis"        # "equal", "daylight_only_risk", or "peak_emphasis"
peak_hours   <- 16:21          # used when weight_mode = "peak_emphasis"
peak_weight  <- 2.0
night_cf_threshold <- 0.01     # row max CF < threshold => night (for weighting only)

gamma_grid <- 10^seq(-3, 3, length.out = 25)  # for efficient frontier

gamma0 <- 1.0                                  # single point solution to print


read_site_long_cf <- function(path, sitename, nameplate_kW) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  stopifnot("Hour" %in% names(df))
  df$Hour <- as.integer(df$Hour)
  
  # Long format: Hour, Year, Output_kW
  df_long <- df %>%
    pivot_longer(-Hour, names_to = "Year", values_to = "Output_kW") %>%
    mutate(Year = as.integer(Year), Site = sitename)
  
  # === Take out negatives  ===
  df_long <- df_long %>%
    mutate(Output_kW = ifelse(Output_kW < 0, 0, Output_kW))
  
  # Convert to Capacity factor
  df_long %>%
    mutate(CF = Output_kW / nameplate_kW) %>%
    select(Year, Hour, Site, CF)
}

build_X <- function(panel_cf) {
  X_wide <- panel_cf %>%
    arrange(Year, Hour, Site) %>%
    pivot_wider(names_from = Site, values_from = CF) %>%
    arrange(Year, Hour)
  
  site_cols <- setdiff(names(X_wide), c("Year", "Hour"))
  
  # 
  list(
    X = as.matrix(X_wide[site_cols]),   # rows = hours×years, cols = sites
    meta = X_wide[, c("Year", "Hour")], # keep Year + Hour info for weighting
    site_names = site_cols
  )
}

compute_mu_Sigma <- function(X, hour_weights = NULL) {
  if (is.null(hour_weights)) {
    mu <- colMeans(X)
    Sigma <- stats::cov(X)
  } else {
    stopifnot(length(hour_weights) == nrow(X))
    w <- as.numeric(hour_weights); w <- w / sum(w)
    mu <- as.numeric(colSums(X * w))
    Xc <- sweep(X, 2, mu, "-")
    Sigma <- t(Xc) %*% (Xc * w)  # weighted second moment about mean
  }
  list(mu = mu, Sigma = Sigma)
}

make_hour_weights <- function(X, meta,
                              mode = "equal",
                              night_cf_threshold = 0.01,
                              peak_hours = 18:21,
                              peak_weight = 5.0) {
  w <- rep(1, nrow(X))
  if (mode == "equal") return(w)
  
  # identifying night hours
  row_max_cf <- apply(X, 1, max, na.rm = TRUE)
  is_night <- row_max_cf < night_cf_threshold
  
  # turn 8760 into the hour of the day 
  hour_of_day <- ((meta$Hour - 1) %% 24) + 1
  
  if (mode == "daylight_only_risk") {
    w <- as.numeric(!is_night)
    return(w)
  }
  
  if (mode == "peak_emphasis") {
    w[hour_of_day %in% peak_hours] <- peak_weight
    return(w)
  }
  
  stop("Unknown weight_mode.")
}

#MV WOOHOOOOOO finally 
solve_markowitz <- function(mu, Sigma, gamma = 1.0, lb = NULL, ub = NULL) {
  n <- length(mu)
  w <- Variable(n)
  obj <- t(mu) %*% w - gamma * quad_form(w, Sigma)
  cons <- list(w >= 0, sum(w) == 1)
  if (!is.null(lb)) cons <- c(cons, list(w >= lb))
  if (!is.null(ub)) cons <- c(cons, list(w <= ub))
  prob <- Problem(Maximize(obj), cons)
  res <- solve(prob)
  list(status = res$status,
       w = as.numeric(res$getValue(w)),
       obj = res$value,
       gamma = gamma)
}

solve_target_return <- function(mu, Sigma, r_target, lb = NULL, ub = NULL) {
  n <- length(mu)
  w <- Variable(n)
  obj <- quad_form(w, Sigma)
  cons <- list(w >= 0, sum(w) == 1, t(mu) %*% w >= r_target)
  if (!is.null(lb)) cons <- c(cons, list(w >= lb))
  if (!is.null(ub)) cons <- c(cons, list(w <= ub))
  prob <- Problem(Minimize(obj), cons)
  res <- solve(prob)
  list(status = res$status,
       w = as.numeric(res$getValue(w)),
       variance = as.numeric(res$value),
       r = as.numeric(t(mu) %*% res$getValue(w)))
}

# ----------------------------
# PIPELINE
# ----------------------------
site_names <- site_files %>% basename() %>% tools::file_path_sans_ext()

#Read all sites; floor negatives to 0; convert to CF
panel_cf <- purrr::map2_dfr(
  site_files, site_names,
  ~ read_site_long_cf(.x, .y, nameplate_kW = nameplate_kW)
)

# year x site matrix
bx <- build_X(panel_cf)
X    <- bx$X
meta <- bx$meta
sites <- bx$site_names

# 3) hourly weights
hour_weights <- make_hour_weights(
  X, meta,
  mode = weight_mode,
  night_cf_threshold = night_cf_threshold,
  peak_hours = peak_hours,
  peak_weight = peak_weight
)

# 4) Compute mu and Sigma in CF space
#    If weight_mode == "equal", pass NULL (equal weighting).
ms <- compute_mu_Sigma(X, hour_weights = if (weight_mode == "equal") NULL else hour_weights)
mu <- ms$mu
Sigma <- ms$Sigma

# 5) Solve one Markowitz for a specific gamme (risk tolerance)
sol0 <- solve_markowitz(mu, Sigma, gamma = gamma0)
weights0 <- tibble(Site = sites, Weight = sol0$w) %>% arrange(desc(Weight))
cat("\n--- Markowitz solution (gamma =", gamma0, ") ---\n")
print(weights0)
cat("Expected CF =", sum(mu * sol0$w), "  Variance =", as.numeric(t(sol0$w) %*% Sigma %*% sol0$w), "\n")

# 6) Efficient frontier 
front <- purrr::map(gamma_grid, ~ solve_markowitz(mu, Sigma, gamma = .x))
frontier <- tibble(
  gamma     = sapply(front, `[[`, "gamma"),
  exp_CF    = sapply(front, function(s) sum(mu * s$w)),
  variance  = sapply(front, function(s) as.numeric(t(s$w) %*% Sigma %*% s$w))
) %>% arrange(variance)

cat("\n--- Efficient frontier (first 6 rows) ---\n")
print(head(frontier, 6))

# 7) Target-return example (choose the 60th percentile expected CF on the frontier)
r_target <- quantile(frontier$exp_CF, 0.60)
tr <- solve_target_return(mu, Sigma, r_target = r_target)
tr_weights <- tibble(Site = sites, Weight = tr$w) %>% arrange(desc(Weight))
cat("\n--- Target-return solution (target exp_CF =", round(r_target, 4), ") ---\n")
print(tr_weights)
cat("Achieved exp_CF =", round(tr$r, 4), "  Variance =", round(tr$variance, 6), "\n")

plot(frontier$variance, frontier$exp_CF, pch=16,
     xlab="Variance of CF", ylab="Expected CF",
      main=paste0("Efficient Frontier (weight_mode=", weight_mode, ")"))
