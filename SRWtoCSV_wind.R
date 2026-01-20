library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(janitor)

# --------------------------------
# USER PATHS
# --------------------------------
in_dir  <- "/Users/nicktennes/Documents/ERA5 Weather Files RAW/ERA5 Wind Weather"
out_dir <- "/Users/nicktennes/Documents/ERA5 Weather Files CLEAN"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# --------------------------------
# Parse filename: site + year
# Ignores clean/corrected/etc
# --------------------------------
parse_name <- function(path) {
  fn <- basename(path)
  
  m <- str_match(fn, "^(?i)wind_(.+)_((?:19|20)\\d{2})_.*\\.srw$")
  if (any(is.na(m))) stop("Filename does not match expected pattern: ", fn)
  
  tibble(site = m[, 2], year = as.integer(m[, 3]))
}

# --------------------------------
# Read SRW: grab header row (line 3), then read data from line 6+
# --------------------------------
read_srw <- function(path) {
  meta  <- parse_name(path)
  lines <- read_lines(path, n_max = 5)
  
  # line 3 contains column names like: Temperature,Pressure,Direction,Speed
  header <- str_split(lines[3], ",", simplify = TRUE) |> as.character()
  header <- str_trim(header)
  
  # Data start at line 6 (skip 5 lines)
  df <- read_csv(
    path,
    skip = 5,
    col_names = header,
    col_types = cols(
      .default = col_double()
    ),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    clean_names() %>%
    mutate(
      site = meta$site,
      year = meta$year,
      hour = row_number()
    ) %>%
    relocate(site, year, hour)
  
  df
}

# --------------------------------
# Run on all SRW files
# --------------------------------
files <- list.files(in_dir, pattern = "\\.srw$", full.names = TRUE)

all_long <- map_dfr(files, read_srw)

# --------------------------------
# Write ONE long CSV per site
# --------------------------------
all_long %>%
  group_by(site) %>%
  group_walk(~ {
    out_path <- file.path(out_dir, paste0(.y$site, "_wind.csv"))
    write_csv(.x, out_path)
  })

# Optional: combined file
write_csv(all_long, file.path('/Users/nicktennes/Documents/ERA5 Weather Files CLEAN', "ALL_SITES_wind.csv"))
