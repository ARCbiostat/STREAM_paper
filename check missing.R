# Set your folder path
path <- "Simulation/models_results/low/streams"

# List all files
files <- list.files(path)

# Extract numbers from filenames
nums_in_folder <- as.numeric(gsub("[^0-9]", "", files))

# Expected numbers
expected <- 1:200

# Find missing
missing_nums <- setdiff(expected, nums_in_folder)


