# Set your folder path
folder <- "/Users/nicktennes/SAM Downloaded Weather Files/Grand Canyon Junction, AZ"

# List all csv files in the folder
files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)

# Create new names by adding a prefix
new_names <- file.path(folder, paste0("GCJunction_", basename(files)))

# Add suffix before the .csv extension
new_names <- sub("\\.csv$", "_extra.csv", files)

# Rename the files
file.rename(files, new_names)
