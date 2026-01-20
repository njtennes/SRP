library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(janitor)

# --------------------------------
# USER PATHS
# --------------------------------
in_dir  <- '/Users/nicktennes/Documents/ERA5 Weather Files RAW/ERA5 Solar Weather'
out_dir <- '/Users/nicktennes/Documents/ERA5 Weather Files CLEAN'

# --------------------------------
# Parse filename: site + year
# Example:
# Solar_CasaGrande_1992_correctedDNI_30minuteshift_final.csv
# --------------------------------
parse_solar_name <- function(path) {
  fn <- basename(path)
  
  m <- str_match(fn, "^(?i)solar_(.+)_((?:19|20)\\d{2})_.*\\.csv$")
  if (any(is.na(m))) stop("Filename does not match expected pattern: ", fn)
  
  tibble(
    site = m[, 2],
    file_year = as.integer(m[, 3])
  )
}

# --------------------------------
# Read one solar file
# - skip first 2 metadata lines
# - line 3 becomes header
# --------------------------------
read_solar <- function(path) {
  meta <- parse_solar_name(path)
  
  df <- readr::read_csv(
    path,
    skip = 2,
    col_names = TRUE,
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    janitor::clean_names() %>%
    mutate(
      site = meta$site,
      file_year = meta$file_year
    ) %>%
    relocate(site, file_year)
  
  # Optional: hour_of_year within each file year (helpful later)
  # If you already trust (year,month,day,hour), you can remove this.
  df <- df %>%
    group_by(site, file_year) %>%
    mutate(hour_of_year = row_number()) %>%
    ungroup() %>%
    relocate(hour_of_year, .after = file_year)
  
  df
}

# --------------------------------
# Run on all solar files
# --------------------------------
files <- list.files(in_dir, pattern = "(?i)^solar_.*\\.csv$", full.names = TRUE)

solar <- purrr::map_dfr(files, read_solar)

# --------------------------------
# Write ONE long CSV per site
# --------------------------------
solar %>%
  group_by(site) %>%
  group_walk(~ {
    out_path <- file.path(out_dir, paste0(.y$site, "_solar.csv"))
    write_csv(.x, out_path)
  })

# --------------------------------
# Write master file (all sites)
# --------------------------------
write_csv(solar, file.path(out_dir, "ALLSITES_solar.csv"))
