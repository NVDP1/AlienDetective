library(dplyr)
library(tidyr)


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
      filtered_rows <- filter(long_sea, location == "TZS")
      
      # Optionally add a column with the species/folder name
      filtered_rows$Species <- basename(folder)
      
      # Add to list
      filtered_data_list[[length(filtered_data_list) + 1]] <- filtered_rows
    }
  }
}

# Combine all filtered data into one data frame
all_filtered_data <- bind_rows(filtered_data_list)

######     SAVE FILTERED CSVs AND PLOTS IN LOCATION FOLDERS     ############
library(ggplot2)

# Get all unique locations from the combined filtered data
locations <- unique(all_filtered_data$location)

location_graphs_folder <- file.path(main_folder, "Location_graphs")
if (!dir.exists(location_graphs_folder)) {
  dir.create(location_graphs_folder)
}

# Loop through unique locations
for (loc in locations) {
  # Clean the location name
  safe_loc <- gsub("[^a-zA-Z0-9_-]", "_", loc)
  
  # Create a folder for the location inside Location_graphs
  location_folder <- file.path(location_graphs_folder, safe_loc)
  if (!dir.exists(location_folder)) {
    dir.create(location_folder, recursive = TRUE)
  }
  
  # Subset data for this location
  loc_data <- subset(all_filtered_data, location == loc)
  loc_data <- loc_data[!is.na(loc_data$x), ]
  loc_data <- loc_data[!is.infinite(loc_data$x),]
  # Save filtered CSV
  csv_file_path <- file.path(location_folder, paste0("filtered_", safe_loc, ".csv"))
  write.csv(loc_data, csv_file_path, row.names = FALSE)
  
  # Create and save plot
  plot <- ggplot(loc_data, aes(x = x, fill = Species)) +
    geom_histogram(aes(y = after_stat(count / sum(count) * 100)),binwidth = 50, boundary = 0, position = "stack") +
    labs(title = paste0("Distribution of different species that are present in ", loc),
         x = "Sea distance (km)", y = "Percentage of total counts (%)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.key.size = unit(1, "lines"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_x_continuous(breaks = seq(0, 8000, by = 500), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 8000), ylim = c(0,20))
  
  plot_file_path <- file.path(location_folder, paste0("Distribution_of_species_in_", safe_loc, ".png"))
  ggsave(filename = plot_file_path, plot = plot, width = 2400, height = 1200, units = "px", dpi = 300)
}
max(loc_data$x)
