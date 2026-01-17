library(dplyr)
library(readr)
library(stringr)
library(purrr)

#------------------------------------------------------------
# 1) Set folders
#------------------------------------------------------------
infolder  <- '/Users/nicktennes/Documents/SAM ERA5 Output 8760 dirty/Solar'   # dirty files
outfolder <- '/Users/nicktennes/Documents/SAM ERA5 Output 8760'   # cleaned files

dir.create(outfolder, showWarnings = FALSE)

#------------------------------------------------------------
# 2) List all CSVs (e.g., Deming.csv, StJohns.csv, etc.)
#------------------------------------------------------------
files <- list.files(infolder, pattern = "\\.csv$", full.names = TRUE)

#------------------------------------------------------------
# 3) Process each file
#------------------------------------------------------------
walk(files, function(f) {
  
  message("Cleaning: ", basename(f))
  
  # ---- Read dirty file WITH header row (run,1,2,3,...) ----
  raw <- read_csv(f, show_col_types = FALSE)
  
  # First column name (usually "run")
  first_col <- names(raw)[1]
  
  # ---- Row 1 (after header) has the long filenames & metadata ----
  header_row <- raw[1, ]
  
  # Build cleaned header:
  # - first column becomes "Hour"
  # - other columns: extract 4-digit year from filenames
  header_clean <- header_row %>%
    mutate(
      !!first_col := "Hour",
      across(-all_of(first_col), ~ {
        x_chr <- as.character(.x)
        yr <- str_extract(x_chr, "(19|20)\\d{2}")
        ifelse(is.na(yr), x_chr, yr)
      })
    )
  
  new_names <- as.character(header_clean[1, ])
  
  # ---- Drop the metadata row we just used as header source ----
  data <- raw[-1, ]
  names(data) <- new_names
  
  # ---- Identify year columns (everything except "Hour") ----
  other_cols <- setdiff(names(data), "Hour")
  
  year_nums <- suppressWarnings(as.integer(other_cols))
  valid     <- !is.na(year_nums)
  
  if (!any(valid)) {
    warning("No valid year columns found in file: ", basename(f),
            " — leaving as Hour + all other columns.")
    data_clean <- data %>%
      mutate(Hour = row_number())
  } else {
    year_cols <- other_cols[valid][order(year_nums[valid])]
    
    # Keep Hour + year columns; drop rows that are all NA in year cols
    data_clean <- data %>%
      select(Hour, all_of(year_cols)) %>%
      filter(!if_all(all_of(year_cols), ~ is.na(.))) %>%
      mutate(Hour = row_number())   # reindex 1..n (ideally 1..8760)
  }
  
  # ---- Write cleaned file to output folder ----
  outfile <- file.path(outfolder, basename(f))
  write_csv(data_clean, outfile)
})

