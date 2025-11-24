library(dplyr)
library(readr)
library(purrr)
library(stringr)

# 1) Set your folder path
folder <- '/Users/nicktennes/Documents/SAM Solar Power 8760/Wilcox'  # <- change this
outfolder <- '/Users/nicktennes/Documents/SAM Solar Power 8760/Locations Aggregated'

# 2) List all CSVs that look like Location_YYYY.csv
files <- list.files(folder, pattern = "Wilcox_\\d{4}\\.csv$", full.names = TRUE)

# 3) Sort files by year so columns appear in ascending order
years <- str_extract(basename(files), "\\d{4}") %>% as.integer()
ord <- order(years)
files <- files[ord]
years <- years[ord]

# 4) Read each file and rename columns:
#    - "Time stamp" -> "Timestamp"
#    - "System power generated | (kW)" -> "<year>"
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
  arrange(Hour) %>%
  select(Hour, all_of(as.character(years)))

# 6) Optional: quick sanity check for duplicate timestamps
# any(duplicated(merged$Timestamp))

# 7) Write out the merged CSV
out_path <- file.path(outfolder, "Wilcox.csv")
write_csv(merged, out_path)
