# Master-thesis-Assessment of structural and topographic influences on greenhouse gas flux observations using a footprint model and LIDAR measurements in a peatland ecosystem-
Assessment of structural and topographic influences on greenhouse gas flux observations using a footprint model and LIDAR measurements in a peatland ecosystem is the title of the Master thesis
R scripts for all the performed analysis and visualization for the Master thesis.

Descriptions of the scripts:

Master_thesis.rmd- contains all the codes for the analysis and visualizations performed during the Master thes

Master_thesis.html- contains the knitted html version of the rmd file.

spatial_FFP- contains the modified FFP footprint model which was integrated to accept a spatial roughness parameter to compute the footprints within the rmd file

user_defined_functions- contains codes for other user defined functions that were useful in the analysis and visualizations in the rmd script.



To run the original model by Kljun et al. (2015), the script can be obtained from https://footprint.kljun.net/ which was used to compute the footprints within the rmd file



Descriptions of the data:

Ec_tower- point shapefile of the eddy covariance tower

field_poly- polygon shapefile of the Amtsvenn site 

rf_training- training data for the random forest model

wind_direction- shapefile of the wind sectors of directions with a buffer distance

flux_data- Eddy covariance data used in the study
