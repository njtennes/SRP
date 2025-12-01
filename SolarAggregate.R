library(dplyr)
library(readr)
library(purrr)
library(stringr)

# 1) Set your folder path
folder <- '/Users/nicktennes/Desktop/ERA5 Output Dirty/Silver City'  # <- change this
outfolder <- '/Users/nicktennes/Documents/SAM ERA5 Output 8760'

# 2) List all CSVs that look like Location_YYYY.csv
files <- list.files(folder, pattern = "SilverCity_\\d{4}\\.csv$", full.names = TRUE) # <- change this

# 3) Sort files by year so columns appear in ascending order
years <- str_extract(basename(files), "\\d{4}") %>% as.integer()
ord <- order(years)
files <- files[ord]
years <- years[ord]

read_and_label <- function(file) {
  yr <- str_extract(basename(file), "\\d{4}")
  read_csv(file, show_col_types = FALSE) %>%
    rename(
      Hour = `Time stamp`,
      !!yr := `System power generated | (kW)`
    )
}

dfs <- map(files, read_and_label)

# 5) Merge on Timestamp (full join to keep all stamps, if any small mismatches)
merged <- reduce(dfs, full_join, by = "Hour") %>%
  mutate(
    # Parse strings like "Jan 1, 12:00 am" as a datetime
    Hour_dt = as.POSIXct(
      paste("2001", Hour),            
      format = "%Y %b %e, %I:%M %p",
      tz = "UTC"
    )
  ) %>%
  arrange(Hour_dt) %>%                # sort in true time order
  mutate(Hour = row_number()) %>%     # 1,2,3,...,8760
  select(Hour, all_of(as.character(years)))

# 7) Write out the merged CSV
out_path <- file.path(outfolder, "Silver City.csv") # <- change this
write_csv(merged, out_path)
