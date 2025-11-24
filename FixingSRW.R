#Cleaning "Broken" SRW Files

library(stringr)

# Folder with your 44 Silver City SRW files
silver_dir <- "/Users/nicktennes/Desktop/ERA5 Wind Weather/Silver City"

# List all .srw files in that folder
silver_files <- list.files(
  path       = silver_dir,
  pattern    = "\\.srw$",
  full.names = TRUE
)

clean_silver_file <- function(path) {
  message("Processing: ", basename(path))
  
  # Read all lines from the file
  lines <- readLines(path)
  
  # Find the line that says "Data for Silver City, NM"
  data_idx <- grep("^Data for Silver City, NM", lines)[1]
  
  if (is.na(data_idx)) {
    warning("  Skipping: could not find 'Data for Silver City, NM' in ",
            basename(path))
    return(NULL)
  }
  
  # Extract the year from the filename (first 4-digit sequence)
  fname <- basename(path)
  year  <- str_extract(fname, "\\d{4}")
  
  if (is.na(year)) {
    warning("  No 4-digit year found in filename: ", fname,
            " (using blank year in header)")
    year <- ""
  }
  
  # Build the new first line (your desired header)
  # loc_id,SilverCity,USA,<year>,32.75,-108.2
  new_header <- sprintf("loc_id,SilverCity,NM,USA,%s,32.75,-108.2,1825", year)
  
  # Keep everything from "Data for Silver City, NM" down
  rest <- lines[data_idx:length(lines)]
  
  # Combine into cleaned file: new header + the rest
  cleaned_lines <- c(new_header, rest)
  
  # Write to a new file with _clean.srw suffix so originals stay untouched
  outpath <- file.path(
    silver_dir,
    paste0(tools::file_path_sans_ext(fname), "_clean.srw")
  )
  
  writeLines(cleaned_lines, outpath)
  
  message("  -> Wrote cleaned file: ", basename(outpath))
}

# Apply to all Silver City SRWs
invisible(lapply(silver_files, clean_silver_file))

cat("\nDone. Check the *_clean.srw files in:\n", silver_dir, "\n")
