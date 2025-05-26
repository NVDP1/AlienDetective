library(dplyr)

# Set the main folder path
main_folder <- "data/Output_M"

# Get list of all species folders (directories only)
species_folders <- list.dirs(path = main_folder, full.names = TRUE, recursive = FALSE)

# Initialize list to collect filtered data
filtered_data_list <- list()

# Loop through each species folder
for (folder in species_folders) {
  # List CSV files in the current species folder
  csv_files <- list.files(path = folder, pattern = "\\.csv$", full.names = TRUE)
  
  for (csv_file in csv_files) {
    # Read the CSV file
    data <- read.csv(csv_file, stringsAsFactors = FALSE )
    long_df <- data %>%
      pivot_longer(
        cols = ends_with("_seaway") | ends_with("_geodesic"),
        names_to = "location",
        values_to = "x"
      )%>%
      separate(location, into = c("location", "DistanceType"), sep = "_")
    long_sea <- long_df %>%
      filter(grepl("seaway", DistanceType))
    long_geo <- long_df %>%
      filter(grepl("geodesic", DistanceType))
    long_sea <- long_sea[!is.na(long_sea$x), ]
    # Check if 'Location' column exists
    if ("location" %in% names(long_sea)) {
      # Filter rows (e.g., Location == "New York")
      filtered_rows <- filter(long_sea, location == "RavennaH")
      
      # Optionally add a column with the species/folder name
      filtered_rows$Species <- basename(folder)
      
      # Add to list
      filtered_data_list[[length(filtered_data_list) + 1]] <- filtered_rows
    }
  }
}

# Combine all filtered data into one data frame
all_filtered_data <- bind_rows(filtered_data_list)

# View or export result
print(all_filtered_data)
write.csv(all_filtered_data, file.path(main_folder, "filtered_species_data.csv"), row.names = FALSE)

