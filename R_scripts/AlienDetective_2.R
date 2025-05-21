# AlienDetective.R
# Main script


#############
### SETUP ###
#############

setwd("~/AlienDetective")
source("src/functions.R")

# Reset graphics settings
graphics.off()

.libPaths(c("/cfs/klemming/home/n/nielsvdp/Calculation_Distance/packages",
            "/cfs/klemming/pdc/software/dardel/23.12/eb/software/R/4.4.1-cpeGNU-23.12/lib64/R/library"))

# NB! No packages loaded here, only installed if missing. Better to use explicit namespaces instead [e.g. raster::extract() rather than just extract()].
# That way it's easier to maintain the code and see which packages are actually required as development progresses, and you also avoid clashes between
# package namespaces, making sure that the correct function is always used regardless of which other packages the user has installed and loaded.
message(">>> [INIT] Checking for required packages...")
packages <- c("rgbif", "sf", "sp", "gdistance", "geodist", "raster", "fasterize", "ggplot2", "rnaturalearth", "rnaturalearthdata","geosphere")
for (package in packages) {
  if(!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

setwd("/cfs/klemming/home/n/nielsvdp/Calculation_Distance")
# If installation fails for a package because it is "not available for your version of R", try installing from source
# install.packages(package, pkgType = "source")
# Install the rnaturalearthhires package for a higher-resolution vector map
# install.packages("devtools")
# remotes::install_github("ropensci/rnaturalearthhires")

message(">>> [INIT] Reading input data...")
args <- commandArgs(trailingOnly = TRUE)
species_location_path <- args[1]  # First argument: path to species file
location_coordinates_path <- args[2]  # Second argument: path to coordinates file
rasterized_path <- args[3]  # Third argument: path to rasterized sf-formatted world map (RDS file)
cost_matrix_path <- args[4] # Fourth argument: path to transition matrix (RDS file)
output_dir <- args[5]  # Fifth argument: output directory

# Set file paths
if (length(args) == 0) {
  species_location_path <- file.path("data", "Species_Location_NIS_Space.csv")
  location_coordinates_path <- file.path("data", "Coordinates_NIS.csv")
  land_polygons_path <- file.path("data", "land_polygons.shp")
  rasterized_path <- file.path("data", "rasterized_land_polygons.rds")
  cost_matrix_path <- file.path("data", "cost_matrix.rds")
  output_dir <- file.path("data", "output")
}

# if (!file.exists(land_polygons_path) && !file.exists(rasterized_path)) {
#   stop("Path to either an sf-formatted rasterized rds file or a polygon vector shape file must be provided.")
# }

# Read species-location presence/absence matrix
species_location <- read.csv(species_location_path, sep = ";")
# If there are more than one row per species, remove all but the first row for each species
species_location <- species_location[!duplicated(species_location$Specieslist),]
# Read table of coordinates for every location name (ObservatoryID)
location_coordinates <- read.csv(location_coordinates_path, sep = ";")

# INSERT LIST OF NATIVE SPECIES TO REMOVE NATIVE SPECIES FROM DF LIST

# Subselect species to run the script for (optional). Can also be used to exclude species, e.g. known natives, by negating the which function
species_subset <- c("Celleporaria brunnea")
species_location <- species_location[which(species_location$Specieslist %in% species_subset),]
#species_location <- species_location[c(2, 10, 57),] # Or subset a few species to try at random

required_columns <- c("decimalLatitude", "decimalLongitude", "year", "month", "country")


#########################
### MAP CONFIGURATION ###
#########################

# Load rasterized world map if it exists, otherwise load custom vector shapefile and rasterize it
message(">>> [MAP] Loading world map...")
land_shapefile <- file.path("data","land_polygons.shp")
if(file.exists(rasterized_path)) {
  r <- readRDS(rasterized_path)
} else {
  # Read vector map as sf object
  #land_polygons <- sf::st_read(land_polygons_path)
  land_polygons <- st_read(land_shapefile, quiet = TRUE)
  st_crs(land_polygons) <- 4326
  message(">>> [MAP] Rasterizing land polygons...")
  # Create raster
  r <- raster::raster(raster::extent(-180, 180, -90, 90), crs = sp::CRS("+init=EPSG:4326"), resolution = 0.1)
  # Rasterize vector map using fasterize
  r <- fasterize::fasterize(land_polygons, r, field = NULL, fun = "max")
  # Set sea cells to value 1 and land cells to NA (Opposite of what fasterize outputs)
  r <- raster::calc(r, function(x) ifelse(is.na(x), 1, NA))
  saveRDS(r, rasterized_path)
  rm(land_polygons)
  message(">>> [MAP] Rasterization done. Saved raster to \"", file.path(getwd(), rasterized_path), "\"")
}

if (file.exists(cost_matrix_path)) {
  message(">>> [MAP] Loading cost matrix...")
  cost_matrix <- readRDS(cost_matrix_path)
} else {
  message(">>> [MAP] Generating cost matrix...")
  # Create a transition object for adjacent cells
  cost_matrix <- gdistance::transition(r, transitionFunction = mean, directions = 16)
  # Set infinite costs to NA to prevent travel through these cells
  cost_matrix <- gdistance::geoCorrection(cost_matrix, type = "c", scl = FALSE)
  # Save transition matrix
  saveRDS(cost_matrix, file = cost_matrix_path)
  message(">>> [MAP] Saved cost matrix to \"", file.path(getwd(), cost_matrix_path), "\"")
}


#############################
### DISTANCES CALCULATION ###
#############################

for (species in species_location[,1]) {
  species_dir <- file.path(output_dir, gsub(" ", "_", species))
  gbif_occurrences_file <- file.path(species_dir, paste0(gsub(" ", "_", species), ".csv"))
  if (file.exists(gbif_occurrences_file)) {
    message(">>> [GBIF] Loading GBIF data for ", species)
    gbif_occurrences <- read.csv(gbif_occurrences_file, header = TRUE)
  } else {
    message(">>> [GBIF] Fetching GBIF data for ", species)
    gbif_occurrences <- fetch_gbif_data(species, fields = required_columns)
    if (!is.null(gbif_occurrences)) {
      if (!dir.exists(species_dir)) {
        dir.create(species_dir, recursive = TRUE)
      }
      write.csv(gbif_occurrences, file = gbif_occurrences_file, row.names = FALSE)
    } else {
      # Skip species that have no GBIF records
      next
    }
  }
  
  for (location in colnames(species_location[,-1])) {
    # Skip locations where the species hasn't been detected, determined by a read number cutoff (default 1 read).
    # Ideally, data from multiple marker genes should have been compiled into a single presence/absence table before, so there should only be 1 or 0.
    if (species_location[which(species_location[,1] == species), location] < 1) {
      next
    }
    
    # Get coordinates for the location of observation
    lat_raw <- location_coordinates[which(location_coordinates$Observatory.ID == location), "Latitude"]
    lon_raw <- location_coordinates[which(location_coordinates$Observatory.ID == location), "Longitude"]
    
    # Replace comma with dot and convert to numeric
    latitude <- as.numeric(gsub(",", ".", lat_raw))
    longitude <- as.numeric(gsub(",", ".", lon_raw))
    
    if (is.na(latitude) || is.na(longitude)) {
      warning("Invalid or missing coordinates for location: ", location)
      next
    }
    if (length(latitude) != 1 || length(longitude) != 1) {
      warning("Could not retrieve coordinates for location \"", location, "\"")
      next
    }
    
    # Run the distance calculations
    # Ensure query point is in the sea
    corrected_coords <- ensure_point_in_sea(latitude, longitude, r)
    if (any(is.na(corrected_coords))) {
      warning(">>> [DIST] Skipping location ", location, " for ", species, " – could not move point into sea")
      next
    }
    latitude <- corrected_coords[1]
    longitude <- corrected_coords[2]
    
    # Run the distance calculations
    message(">>> [DIST] Calculating distances to ", species, " occurrences from ", location)
    result <- calculate.distances(gbif_occurrences = gbif_occurrences,
                                  latitude = latitude,
                                  longitude = longitude,
                                  raster_map = r,
                                  cost_matrix = cost_matrix)
    
    if (!(is.null(result$seaway) & is.null(result$geodesic))) {
      gbif_occurrences[,paste0(location, "_seaway")] <- result$sea_distances
      gbif_occurrences[,paste0(location, "_geodesic")] <- result$geodesic_distances
    } else {
      gbif_occurrences[,paste0(location, "_seaway")] <- NA
      gbif_occurrences[,paste0(location, "_geodesic")] <- NA
    }
  }
  
  # Save to csv file
  write.csv(gbif_occurrences, file = gbif_occurrences_file, row.names = FALSE)
  message("")
}


################
### PLOTTING ###
################

library("ggplot2")
library(tidyr)
library(dplyr)
required_columns <- c("latitude", "longitude", "year", "month", "country", "basisOfRecord")

# Iterate over species for which a csv file exists
for (species in species_location[,1]) {
  species_ <- gsub(" ", "_", species)
  species_dir <- file.path(output_dir, species_)
  if (file.exists(file.path(species_dir, paste0(species_, ".csv")))) {
    distance_df <- read.csv(file.path(species_dir, paste0(species_, ".csv")))
    long_df <- distance_df %>%
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
    
  } else {
    warning("No output directory found for species \"", species, "\". Skipping plotting.")
    next
  }
  
    # Assign year categories
    year_categories <- c("1965-1985", "1985-1990", "1990-1995",
                         "1995-2000", "2000-2005", "2005-2010",
                         "2010-2015", "2015-2020", "2020-2025")
    long_sea$year_category <- sapply(long_sea$year, assign_year_category)
    long_geo$year_category <- sapply(long_sea$year, assign_year_category)
    # clean dataframe from rows with Inf in them
    long_sea <- long_sea[is.finite(long_sea$x), ]
    # select only distances below 40000km
    long_sea <- subset(long_sea, x < 40000)
    
    for (loc in unique(long_sea$location)) {
      sea_loc_data <- long_sea[long_sea$location == loc, ]
      geo_loc_data <- long_geo[long_geo$location == loc, ]
      
      combined_distances <- rbind(sea_loc_data, geo_loc_data)
      
      # Plot functions by location
      plot_dist_sea <- plot.dist.sea(
        species = species,
        location = loc,
        distances = sea_loc_data$x,
        output_dir = species_dir
      )
      
      plot_both <- plot.dist.both(
        species = species,
        location = loc,
        distances = combined_distances$x,
        output_dir = species_dir
      )
      
      plot_country <- plot.dist.by.country(
        species = species,
        location = loc,
        distances = sea_loc_data$x,
        output_dir = species_dir
      )
      
      plot_year <- plot.dist.by.year(
        species = species,
        location = loc,
        distances = sea_loc_data$x,
        output_dir = species_dir
      )
    }
}
plot_dist_sea
plot_both
plot_country
plot_year


library(ggplot2)
library(gganimate)

ggplot(long_sea, aes(x = longitude, y = latitude)) +
  geom_point(aes(color = location, shape = location), size = 3, alpha = 0.8) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = paste0("Spatiotemporal spread of ", species, " Year: {frame_time}"), x = "Longitude", y = "Latitude") +
  transition_time(year) +
  shadow_mark(past = TRUE, alpha = 0.3)
animate(last_plot(), renderer = gifski_renderer(), width = 800, height = 600)
anim_save("species_spread.gif")


library(dplyr)
library(geosphere)  # for distHaversine
library(ggplot2)

long_sea$species_name <- "Balanus trigonus"
# Step 1: Compute annual centroids
centroids <- long_sea %>%
  group_by(species_name, year) %>%
  summarise(
    mean_lat = mean(latitude, na.rm = TRUE),
    mean_lon = mean(longitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(species, year)

# Step 2: Calculate distance between consecutive centroids per species
centroids <- centroids %>%
  group_by(species_name) %>%
  mutate(
    prev_lat = lag(mean_lat),
    prev_lon = lag(mean_lon),
    dist_km = distHaversine(
      cbind(prev_lon, prev_lat),
      cbind(mean_lon, mean_lat)
    ) / 1000  # meters to km
  ) %>%
  ungroup()
ggplot(centroids, aes(x = year, y = dist_km, group = species, color = species_name)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Centroid Movement Over Time",
    x = "Year",
    y = "Distance moved (km)"
  ) +
  theme_minimal() +
  scale_color_viridis_d(option = "plasma") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )


