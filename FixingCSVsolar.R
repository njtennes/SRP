############
########
#Deming
########
############

library(stringr)

deming_dir <- "/Users/nicktennes/Desktop/ERA5 Solar Weather/Deming"  # <-- change if needed

# List all CSVs 
deming_files <- list.files(
  path       = deming_dir,
  pattern    = "\\.csv$",
  full.names = TRUE
)

clean_deming_file <- function(path) {
  message("Processing: ", basename(path))
  
  lines <- readLines(path)
  
  # Find the line that starts the actual data header
  header_idx <- grep("^year,month,day,hour,ghi,dni,dhi,tdry,wspd,wdir,pres", lines)[1]
  
  if (is.na(header_idx)) {
    warning("  Skipping: couldn't find the data header in ", basename(path))
    return(NULL)
  }
  
  # Everything from the data header down is the real table
  data_block <- lines[header_idx:length(lines)]
  
  # Deming site metadata 
  header1 <- "Source,Latitude,Longitude,Time Zone,Elevation"
  header2 <- "ERA5,32.5,-108.8,-7,1320"
  
  # Combine into cleaned file
  cleaned_lines <- c(header1, header2, data_block)
  
  outpath <- file.path(
    deming_dir,
    paste0(tools::file_path_sans_ext(basename(path)), "_clean.csv")
  )
  
  writeLines(cleaned_lines, outpath)
  
  message("  -> Wrote cleaned file: ", basename(outpath))
}

# Apply to all Deming CSVs
invisible(lapply(deming_files, clean_deming_file))

cat("\nDone. Check the *_clean.csv files in:\n", deming_dir, "\n")

##############
#########
#Kingman
#########
##############

kingman_dir <- "/Users/nicktennes/Desktop/ERA5 Solar Weather/Kingman" 

# List allCSVs
kingman_files <- list.files(
  path       = kingman_dir,
  pattern    = "\\.csv$",
  full.names = TRUE
)

clean_kingman_file <- function(path) {
  message("Processing: ", basename(path))
  
  lines <- readLines(path)
  
  # Find the line that starts the actual data header
  header_idx <- grep("^year,month,day,hour,ghi,dni,dhi,tdry,wspd,wdir,pres", lines)[1]
  
  if (is.na(header_idx)) {
    warning("  Skipping: couldn't find the data header in ", basename(path))
    return(NULL)
  }
  
  # Everything from the data header down is the real table
  data_block <- lines[header_idx:length(lines)]
  
  # Kingman site metadata 
  header1 <- "Source,Latitude,Longitude,Time Zone,Elevation"
  header2 <- "ERA5,35.25,-114.2,-7,1016"
  
  # Combine into cleaned file
  cleaned_lines <- c(header1, header2, data_block)
  
  outpath <- file.path(
    kingman_dir,
    paste0(tools::file_path_sans_ext(basename(path)), "_clean.csv")
  )
  
  writeLines(cleaned_lines, outpath)
  
  message("  -> Wrote cleaned file: ", basename(outpath))
}

# Apply to all Deming CSVs
invisible(lapply(kingman_files, clean_kingman_file))

cat("\nDone. Check the *_clean.csv files in:\n", kingman_dir, "\n")

##############
#########
#Wilcox
#########
##############

wilcox_dir <- "/Users/nicktennes/Desktop/ERA5 Solar Weather/Wilcox" 

# List allCSVs
wilcox_files <- list.files(
  path       = wilcox_dir,
  pattern    = "\\.csv$",
  full.names = TRUE
)

clean_wilcox_file <- function(path) {
  message("Processing: ", basename(path))
  
  lines <- readLines(path)
  
  # Find the line that starts the actual data header
  header_idx <- grep("^year,month,day,hour,ghi,dni,dhi,tdry,wspd,wdir,pres", lines)[1]
  
  if (is.na(header_idx)) {
    warning("  Skipping: couldn't find the data header in ", basename(path))
    return(NULL)
  }
  
  # Everything from the data header down is the real table
  data_block <- lines[header_idx:length(lines)]
  
  # Wilcox site metadata 
  header1 <- "Source,Latitude,Longitude,Time Zone,Elevation"
  header2 <- "ERA5,32.25,-109.8,-7,1271"
  
  # Combine into cleaned file
  cleaned_lines <- c(header1, header2, data_block)
  
  outpath <- file.path(
    wilcox_dir,
    paste0(tools::file_path_sans_ext(basename(path)), "_clean.csv")
  )
  
  writeLines(cleaned_lines, outpath)
  
  message("  -> Wrote cleaned file: ", basename(outpath))
}

# Apply to all Wilcox CSVs
invisible(lapply(wilcox_files, clean_wilcox_file))

cat("\nDone. Check the *_clean.csv files in:\n", wilcox_dir, "\n")

##############
#########
#St. Johns
#########
##############

stjohns_dir <- "/Users/nicktennes/Desktop/ERA5 Solar Weather/St Johns" 

# List allCSVs
stjohns_files <- list.files(
  path       = stjohns_dir,
  pattern    = "\\.csv$",
  full.names = TRUE
)

clean_stjohns_file <- function(path) {
  message("Processing: ", basename(path))
  
  lines <- readLines(path)
  
  # Find the line that starts the actual data header
  header_idx <- grep("^year,month,day,hour,ghi,dni,dhi,tdry,wspd,wdir,pres", lines)[1]
  
  if (is.na(header_idx)) {
    warning("  Skipping: couldn't find the data header in ", basename(path))
    return(NULL)
  }
  
  # Everything from the data header down is the real table
  data_block <- lines[header_idx:length(lines)]
  
  # St. Johns site metadata 
  header1 <- "Source,Latitude,Longitude,Time Zone,Elevation"
  header2 <- "ERA5,34.5,-109.2,-7,1751"
  
  # Combine into cleaned file
  cleaned_lines <- c(header1, header2, data_block)
  
  outpath <- file.path(
    stjohns_dir,
    paste0(tools::file_path_sans_ext(basename(path)), "_clean.csv")
  )
  
  writeLines(cleaned_lines, outpath)
  
  message("  -> Wrote cleaned file: ", basename(outpath))
}

# Apply to all Wilcox CSVs
invisible(lapply(stjohns_files, clean_stjohns_file))

cat("\nDone. Check the *_clean.csv files in:\n", wilcox_dir, "\n")
