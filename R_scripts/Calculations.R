# ------------------------------------------------------------------------------
# Calculate sea and fly distances from ARMS locations to GBIF occurrence points
# ------------------------------------------------------------------------------

# Clear environment
rm(list = ls())
graphics.off()

# ------------------------------------------------------------------------------
# Load Packages
# ------------------------------------------------------------------------------
cat(">>> [INIT] Loading packages...\n")
.libPaths(c("/cfs/klemming/home/n/nielsvdp/Calculation_Distance/packages",
            "/cfs/klemming/pdc/software/dardel/23.12/eb/software/R/4.4.1-cpeGNU-23.12/lib64/R/library"))
suppressMessages({
  library(gdistance)
  library(dplyr)
  library(geosphere)
  library(rgbif)
  library(ggplot2)
  library(tidyr)
  library(worrms)
  library(sf)
  library(sp)
  library(raster)
  library(rnaturalearth)
  library(doParallel)
  library(foreach)
})

# ------------------------------------------------------------------------------
# Setup parallelisation
# ------------------------------------------------------------------------------
# registerDoParallel(cores=10)

# ------------------------------------------------------------------------------
# Set paths and parameters
# ------------------------------------------------------------------------------
setwd("/cfs/klemming/home/n/nielsvdp/Calculation_Distance")
input_species_file <- "Species_Location_NIS_Space.csv"
input_coordinates_file <- "Coordinates_NIS.csv"
output_directory <- "Output_calculations"
land_shapefile <- "land_polygons.shp"

# Create required output dirs
dir.create("OccurrenceData", showWarnings = FALSE)
dir.create(file.path(output_directory, "sea_distances"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_directory, "fly_distances"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_directory, "errors"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
cat(">>> [DATA] Reading species and coordinate files...\n")
df <- read.csv(input_species_file, sep = ";")
Coordinates <- read.csv(input_coordinates_file, sep = ";")

# Convert species-location to long format
long <- df %>%
  pivot_longer(cols = -Specieslist, names_to = "Location", values_to = "Presence") %>%
  filter(Presence > 0)

# Optionally subset species for testing
target_species <- c("Fenestrulina delicia")
long <- long %>% filter(Specieslist %in% target_species)

# ------------------------------------------------------------------------------
# Load land polygons for sea distance masking
# ------------------------------------------------------------------------------
cat(">>> [MAP] Loading land polygons...\n")
land_polygons <- st_read(land_shapefile, quiet = TRUE)
st_crs(land_polygons) <- 4326

# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------

ensure_columns <- function(df, required_cols) {
  for (col in required_cols) {
    if (!col %in% names(df)) {
      df[[col]] <- NA
    }
  }
  return(df)
}

fetch_gbif_occurrences <- function(species_name, required_columns) {
  cat(">>> [GBIF] Fetching occurrences for: ", species_name, "\n")
  basis_types <- c("OBSERVATION", "MACHINE_OBSERVATION", "HUMAN_OBSERVATION", 
                   "MATERIAL_SAMPLE", "LIVING_SPECIMEN", "OCCURRENCE")
  all_data <- list()
  
  for (basis in basis_types) {
    Sys.sleep(2)  # be gentle with GBIF servers
    batch <- tryCatch({
      occ_data(scientificName = species_name, hasCoordinate = TRUE, 
               basisOfRecord = basis, continent = "europe", limit = 10000)$data
    }, error = function(e) {
      cat(">>> [WARN] Could not fetch for basis: ", basis, "\n")
      NULL
    })
    if (!is.null(batch) && nrow(batch) > 0) {
      batch <- ensure_columns(batch, required_columns)
      all_data[[length(all_data)+1]] <- batch[, required_columns]
    }
  }
  
  if (length(all_data) > 0) {
    result <- do.call(rbind, all_data)
    result <- result[!is.na(result$decimalLatitude) & !is.na(result$decimalLongitude), ]
    colnames(result) <- c("Longitude", "Latitude", "Year", "Month", "Country")
    return(result)
  } else {
    return(NULL)
  }
}

# ------------------------------------------------------------------------------
# Main Occurrence Download Loop
# ------------------------------------------------------------------------------

required_columns <- c("decimalLongitude", "decimalLatitude", "year", "month", "country")

for (species in unique(long$Specieslist)) {
# foreach(species = unique(long$Specieslist),
#         .packages = c("rgbif"),
#         .export = c("ensure_columns", "fetch_gbif_occurrences")) %dopar% {
  file_path <- file.path("OccurrenceData", paste0(species, ".csv"))
  if (file.exists(file_path)) {
    cat(">>> [SKIP] File already exists for ", species, "\n")
    next
  }
  
  occ_data <- fetch_gbif_occurrences(species, required_columns)
  
  if (is.null(occ_data) || nrow(occ_data) == 0) {
    error_path <- file.path(output_directory, "errors", paste0("error_", gsub(" ", "_", species), ".txt"))
    writeLines(paste("No usable data for", species), error_path)
    next
  }
  
  write.csv(occ_data, file_path, row.names = FALSE)
  cat(">>> [OK] Wrote GBIF data for ", species, "\n")
}

# ------------------------------------------------------------------------------
# Checking if input Coordinates are on land, if so, move to sea
# ------------------------------------------------------------------------------
cat(">>> [DATA] Checking if input Coordinates are in sea ...")

world <- ne_countries(scale = "medium", returnclass = "sf")  # Load medium or large scale natural earth countries as an sf (simple features) object
r <- raster(extent(-180, 180, -90, 90), res = 0.1)           # Create a raster object with a global extent and resolution of 0.1 degrees
r <- rasterize(world, r, field = 1, fun = max, na.rm = TRUE) # Rasterize the 'world' sf object, assigning a value of 1 to cells with country presence
costs <- reclassify(r, cbind(1, Inf))                        # Reclassify the raster: convert all values of 1 to Inf (infinity)
costs[is.na(costs)] <- 1    # Replace NA values in the 'costs' raster with 1

transition_matrix <- "transitMatrix.rds"
if (!file.exists(transition_matrix)) {
  # Create a transition object for adjacent cells
  transitMatrix <- transition(costs, transitionFunction = function(x) 1/mean(x), directions = 16)
  # Set infinite costs to NA to prevent travel through these cells
  transitMatrix <- geoCorrection(transitMatrix, scl = TRUE)
  # Save/Load transition matrix
  saveRDS(transitMatrix, file = "transitMatrix.rds")
  
} else {
  transitMatrix <- readRDS(file = "transitMatrix.rds")
}

cat(">>> [DATA] Checking if input Coordinates are in sea ...")
for (i in 1:nrow(Coordinates)) {
  cat("\nChecking coordinates for", Coordinates$Observatory.ID[i], ":", Coordinates$Longitude[i], Coordinates$Latitude[i], "\n")
  # Extract coordinates for the current row
  longitude <- as.numeric(gsub(",", ".", Coordinates$Longitude[i]))
  latitude <- as.numeric(gsub(",", ".", Coordinates$Latitude[i]))
  
  # Define the point
  point <- SpatialPoints(cbind(longitude, latitude), proj4string = CRS(proj4string(r)))
  
  # Check if the point is on land or sea
  if (!is.na(raster::extract(r, point))) {
    cat("> Point on land detected, searching for nearest connected sea point...\n")
    
    # Get the transition matrix (sparse)
    trans_matrix <- transitionMatrix(transitMatrix)
    
    # Get all cells that are connected (non-zero transitions)
    connected_cells <- which(rowSums(trans_matrix != 0) > 0)
    connected_coords <- xyFromCell(r, connected_cells)
    
    # Keep only those that fall on sea 
    is_sea <- is.na(raster::extract(r, connected_coords))
    sea_coords <- connected_coords[is_sea, , drop = FALSE]
    
    if (nrow(sea_coords) == 0) {
      warning("No valid connected sea cells found.")
    } else {
      # Compute geodesic distance to all valid sea coords
      dists <- geosphere::distVincentySphere(coordinates(point), sea_coords)
      
      # Find the index of the closest sea coordinate
      nearest_idx <- which.min(dists)
      
      # Update point to the closest connected sea coordinate
      new_coords <- sea_coords[nearest_idx, , drop = FALSE]
      point <- SpatialPoints(new_coords, proj4string = CRS(proj4string(r)))
      
      cat("> Moved point to nearest connected sea at: ", new_coords[1], new_coords[2], "\n")
    }
  } else {
    cat("> Point is already in sea.\n")
  }
  
  # Update the dataframe with the new coordinates (if moved)
  Coordinates$Longitude[i] <- coordinates(point)[1]
  Coordinates$Latitude[i] <- coordinates(point)[2]
}
cat(">>> [DONE] All coordinates updated to nearest sea point")

# ------------------------------------------------------------------------------
# Distance Calculation (Put your existing function here)
# ------------------------------------------------------------------------------

# Iterate over species_name and location_name
Calculation_seadistance <- function(species_name, species_location){
  # Initialize an empty list to store error messages
  error_messages <- list()

  # Define a helper function to add error messages
  add_error_message <- function(message) {
    error_messages <<- c(error_messages, message)
  }
  # Print species name and location
  print(paste0("species_name: ", species_name))
  print(paste0("species_location: ", species_location))

  # Try to read occurrence data
  OccurrenceData <- tryCatch({
    occurrence_coord <- paste0("OccurrenceData/", species_name, ".csv")
    read.csv(occurrence_coord, header = TRUE)
  }, error = function(e) {
    add_error_message(paste("Error reading occurrence data for", species_name, ":", e$message))
    return(NULL)
  })

  if (is.null(OccurrenceData)) return(list(result = NA, error_messages = error_messages))

  # Try to get coordinates for ARMS location
  # use grep to get the row from the Coordinates df where location_name is present
  location_row_index <- grep(species_location, Coordinates$Observatory.ID)
  if (length(location_row_index) == 0) {
    add_error_message(paste("Location not found in Coordinates for", species_location))
    return(list(result = NA, error_messages = error_messages))
  }
  # save longitude and latitude for the row that you selected with grep
  longitude <- as.numeric(Coordinates[location_row_index, "Longitude"])
  latitude <- as.numeric(Coordinates[location_row_index, "Latitude"])
  # make a dataframe out of the longitude and latitude called samplelocation
  samplelocation <- data.frame(Latitude = latitude, Longitude = longitude)

  # check if OccurrenceData has coordinates
  if (nrow(OccurrenceData) < 1) {
    add_error_message("OccurrenceData has no coordinates")
    return(list(result = NA, error_messages = error_messages))
  }

  ##########################################################################
  # DISTANCES CALCULATION
  ##########################################################################

  # Load world data and prepare the raster
  #world <- ne_countries(scale = "large", returnclass = "sf")  # Load medium or large scale natural earth countries as an sf (simple features) object
  #r <- raster(extent(-180, 180, -90, 90), res = 0.1)           # Create a raster object with a global extent and resolution of 0.1 degrees
  #r <- rasterize(world, r, field = 1, fun = max, na.rm = TRUE) # Rasterize the 'world' sf object, assigning a value of 1 to cells with country presence
  #costs <- reclassify(r, cbind(1, Inf))                        # Reclassify the raster: convert all values of 1 to Inf (infinity)
  #costs[is.na(costs)] <- 1    # Replace NA values in the 'costs' raster with 1


  # Initialize lists to store distances
  sea_distances <- c()
  flying_distances <- c()

  # for loop to iterate over OccurrenceData
  for (row in 1:nrow(OccurrenceData)) {
    print(paste0("Calculating latitude: ", OccurrenceData[row, 2], " and longitude: ", OccurrenceData[row, 1]))
    print(paste0("for ", species_location, " latitude, longitude: ", samplelocation[,2], " ", samplelocation[,1]))

    #################
    ## SEA DISTANCE##
    #################

  #  transition_matrix <- "transitMatrix.rds"
   # if (!file.exists(transition_matrix)) {
    #  # Create a transition object for adjacent cells
     # transitMatrix <- transition(costs, transitionFunction = function(x) 1/mean(x), directions = 16)
      # Set infinite costs to NA to prevent travel through these cells
      #transitMatrix <- geoCorrection(transitMatrix, scl = TRUE)
      # Save/Load transition matrix
      #saveRDS(transitMatrix, file = "transitMatrix.rds")

    #} else {
     # transitMatrix <- readRDS(file = "transitMatrix.rds")
    #}


    # Define points using correct projection
    point1 <- SpatialPoints(cbind(samplelocation$Longitude, samplelocation$Latitude), proj4string = CRS(proj4string(r)))
    point2 <- SpatialPoints(cbind(OccurrenceData[row, 1], OccurrenceData[row, 2]), proj4string = CRS(proj4string(r)))
    
    # Check if the OccurrenceData point is on land, if so, move to closest sea
    if (!is.na(raster::extract(r, point2))) {
      cat("Point on land detected, searching for nearest connected sea point...\n")
      
      # Get the transition matrix (sparse)
      trans_matrix <- transitionMatrix(transitMatrix)
      
      # Get all cells that are connected (have non-zero transitions)
      connected_cells <- which(rowSums(trans_matrix != 0) > 0)
      connected_coords <- xyFromCell(r, connected_cells)
      
      # Keep only those that fall on sea 
      is_sea <- is.na(raster::extract(r, connected_coords))
      sea_coords <- connected_coords[is_sea, , drop = FALSE]
      
      if (nrow(sea_coords) == 0) {
        warning("No valid connected sea cells found.")
      } else {
        # Compute geodesic distance to all valid sea coords
        dists <- geosphere::distVincentySphere(coordinates(point2), sea_coords)
        
        # Find the index of the closest sea coordinate
        nearest_idx <- which.min(dists)
        
        # Update point2 to the closest connected sea coordinate
        new_coords <- sea_coords[nearest_idx, , drop = FALSE]
        point2 <- SpatialPoints(new_coords, proj4string = CRS(proj4string(r)))
        
        cat("Moved point to nearest connected sea at: ", new_coords[1], new_coords[2], "\n")
      }
    }
    
    sea_distance <- tryCatch({
      # Coerce points to SpatialPointsDataFrame for compatibility with gdistance
      point1_df <- SpatialPointsDataFrame(coords = point1, data = data.frame(id = 1))
      point2_df <- SpatialPointsDataFrame(coords = point2, data = data.frame(id = 2))
      # Compute cost distance
      cost_distance <- costDistance(transitMatrix, point1_df, point2_df)
      # Calculate the shortest path
      shortest_path <- shortestPath(transitMatrix, point1_df, point2_df, output = "SpatialLines")

      # Plotting the shortest path and the world map
      #plot(r, main = "Shortest Water Path")
      #plot(world, add = TRUE, col = "grey")
      #plot(shortest_path, add = TRUE, col = "blue", lwd = 2)
      #points(point1_df, col = "red", pch = 20)
      #points(point2_df, col = "green", pch = 20)

      # Assuming 'shortest_path' is your SpatialLines object from the shortestPath function
      # First, ensure the CRS is set on the original SpatialLines object
      crs_info <- proj4string(shortest_path)  # or use crs(shortest_path) if using `sp`

      # If it's not set, set it here, assuming the original data was in WGS 84 (EPSG:4326)
      if (is.na(crs_info)) {
        proj4string(shortest_path) <- CRS("+init=epsg:4326")
      }

      # Convert SpatialLines to sf object
      shortest_path_sf <- st_as_sf(shortest_path)

      # Confirm CRS is set for sf object, if not, set it:
      if (is.na(st_crs(shortest_path_sf))) {
        st_crs(shortest_path_sf) <- 4326  # EPSG code for WGS 84
      }

      # Transform to a suitable projected CRS for distance calculation (e.g., UTM zone 33N)
      shortest_path_utm <- st_transform(shortest_path_sf, 32633)  # UTM zone 33N

      # Calculate the length in meters
      path_length <- st_length(shortest_path_utm)
      path_length <- as.numeric(path_length)

      # Print the length
      print(paste0("distance through sea in m: ", path_length))
      sea_distances <- append(sea_distances, path_length)

    }, error = function(e) {
      add_error_message(paste("An error occurred during sea distance calculation for", species_name, "in", species_location, ":", e$message))
      sea_distances <<- append(sea_distances, NA)

    }) # trycatch() closed
  } # iteration over OccurrenceData stopped

  # for loop to iterate again over OccurrenceData for fly distances
  for (row in 1:nrow(OccurrenceData)) {

    # for this calculation, longitude comes first and then latitude!!!
    # Define points using correct projection
    point1 <- SpatialPoints(cbind(samplelocation$Longitude, samplelocation$Latitude), proj4string = CRS(proj4string(r)))
    point2 <- SpatialPoints(cbind(OccurrenceData[row, 1], OccurrenceData[row, 2]), proj4string = CRS(proj4string(r)))

    # Check if the OccurrenceData point is on land
    if (!is.na(raster::extract(r, point2))) {
      add_error_message(paste("Point on land detected for ", species_name, " at ", OccurrenceData[row, 1], "", OccurrenceData[row, 2]))
      flying_distances <- append(flying_distances, Inf)
      next
    }

    #####################
    ## FLYING DISTANCE ##
    #####################

    # Calculate the straight-line distance (accounting for the Earth's curvature)
    straight_line_distance <- distHaversine(coordinates(point1), coordinates(point2))
    print(paste("Straight line distance:", straight_line_distance, "meters"))
    flying_distances <- append(flying_distances, straight_line_distance)

  } # iteration over OccurrenceData stopped

  # Create data frame if lengths match
  create_data_frame <- function(distances, year, month, country) {
    if (length(distances) == length(year) && length(year) == length(month) && length(month) == length(country)) {
      return(data.frame(
        x = distances,
        year = year,
        month = month,
        country = country
      ))
    } else {
      stop("Lengths of vectors do not match. Please ensure all vectors have the same length.")
    }
  }
  ### CHECKS IF NECESSARY ###
   cat("Length of flying_distances: ", length(flying_distances), "\n")
   cat("content of flying_distances: ", flying_distances, "\n")
   cat("Length of sea_distances: ", length(sea_distances), "\n")
   cat("content of sea_distances: ", sea_distances, "\n")
   cat("Length of OccurrenceData$Year: ", length(OccurrenceData$Year), "\n")
   cat("Length of OccurrenceData$Month: ", length(OccurrenceData$Month), "\n")
   cat("Length of OccurrenceData$Country: ", length(OccurrenceData$Country), "\n")

  # Create sea data frame
  sea_data <- create_data_frame(sea_distances, OccurrenceData$Year, OccurrenceData$Month, OccurrenceData$Country)

  # Create fly data frame
  fly_data <- create_data_frame(flying_distances, OccurrenceData$Year, OccurrenceData$Month, OccurrenceData$Country)

  # Define file paths
  sea_distance_file <- paste0("Output_calculations/sea_distances/", species_name, "_distancesTo_", species_location, ".csv")
  fly_distance_file <- paste0("Output_calculations/fly_distances/", species_name, "_distancesTo_", species_location, ".csv")

  # Create directories if they do not exist
  if (!dir.exists("Output_calculations/sea_distances")) {
    dir.create("Output_calculations/sea_distances", recursive = TRUE)
  }

  if (!dir.exists("Output_calculations/fly_distances")) {
    dir.create("Output_calculations/fly_distances", recursive = TRUE)
  }

  # Write data frames to CSV files
  write.csv(sea_data, file = sea_distance_file, row.names = FALSE)
  write.csv(fly_data, file = fly_distance_file, row.names = FALSE)
  # Return the result and the list of error messages
  list(result = list(sea_distances = sea_distances, flying_distances = flying_distances), error_messages = error_messages)
  cat("CSV files are written to Output_calculations")
}

results <- lapply(seq_len(nrow(long)), function(i) Calculation_seadistance(long$Specieslist[i], long$Location[i]))

# results <- foreach(i = seq_len(nrow(long)),
#                    .packages = c("rnaturalearth", "raster", "gdistance", "sp", "geosphere", "sf"),
#                    .export = c("Calculation_seadistance", "Coordinates")) %dopar% {
#                      Calculation_seadistance(long$Specieslist[i], long$Location[i])
#                    }

# stopImplicitCluster()

# ------------------------------------------------------------------------------
# Done
# ------------------------------------------------------------------------------
cat(">>> [DONE] All GBIF occurrence downloads and setup complete.\n")
