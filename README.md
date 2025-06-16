# AlienDetective
Species of maritime fauna all over the world are known to travel great distances in the oceans and seas. The goal of this workflow is to be able to detect them by calculating sea distances (going around land) from a sample location to the occurrence data from this species (from GBIF). This workflow focuses on occurrence data in Europe.
![logo3](https://github.com/IrisVP/AlienDetective/assets/151626670/21dd7508-bd81-448a-a096-db07bace2515)
### Input

The “Input” folder has all input data to execute the R scripts with. The files are in .csv format:
-	Coordinates_NIS.csv: Contains the names of the samplelocations in the column ‘Observatory.ID’, along with their latitudes and longitudes. <br />
- Metadata.csv: Contains metadata per sample, including sample region, country, latitudes, longitudes, sample dates, etc. Note: Metadata should be prepared for personal use in this workflow. The 'Coordinates.csv' file is an example of this preparation. <br />
-	Species_Locations_NIS_Space.csv: Contains a ‘Specieslist’ column with all species names. The other columns are sample location names, with the data in these columns representing absence or presence of a sample taken from a species for a specific location.<br />

### Short description of R scripts

All R scripts can be found in "R_scripts" folder. For detailed documentation: see the 'Documentation' directory for .Rmd and html files

#### Functions.R
/!\ Run this on a server. This is computationally expensive. <br />

Has all functions involved to calculate distances: <br />
•	fetch_gbif_data -> fetches data from gbif <br />
•	ensure_point_in_sea -> checks if point from gbif is in sea or not <br />
•	move_occurrences_to_sea -> if point is detected on land move it to the nearest location in sea <br />
•	calculate.distances -> calculate the sea and geodesic distance <br />
•	Plotting functions

#### AlienDetective_2_M.R
/!\ Run this on a server. This is computationally expensive. <br />
Uses Functions.R script with data from input folder. Processes the input data and put output in correct folders.

##### Output
The output is in the directory 'Output_M' directory. In this directory all alien species folder are made with a csv file of the data and the output graphs. Also a folder "Location_graphs" is found in this directory what shows the graphs and csv files for every ARMS location. 

1. The first graph shows the distribution of one alien species that is present in certain ARMS locations.
2. The second graph shows the sea distances per species and sample location.
3. The third graph is a combined graph of flying distances and sea distances.
4. The fourth graph is similar to the first graph but colored by country.
5. The fifth graph is also similar to the first graph but colored by year category.

Patterns can be observed in these histograms. The x-axis represents distances in km, and the y-axis represents a percentage of the total counts frequency. An alien species can be detected when there is a spike in frequency at a large distance. Using the country plot, the origin region of the species can be identified. Using the year plot, the year in which the species migrated to that location can be determined.

/!\ WARNING /!\ 
Land_polygons file must be downloaded to make the raster map. This is too big to put on the GitHub and also add the cost_matrix.