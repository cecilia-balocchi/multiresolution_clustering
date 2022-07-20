## Aggregate crime data from the csv downloaded for the Philadelphia Police Department

`get_spatial_westphilly.R` creates the shapefiles and spatial objects for the West Philadelphia area, by identifying the region and restricting the original shapefiles (which are saved in `data/shapefile/`), and saves it in `data/westphilly.rdata`.

`get_density_westphilly.R` filters crime data (using the files saved in `data/crime/`,) for the neighborhoods in the West Philadelphia area (using `data/westphilly.rdata`) and calculates crime density as the number of violent crimes per squared meters - note that it will be transformed into the number of violent crimes per squared kilometers before running the analysis. It outputs the csv files used for running the analysis (`data/Ydensity_westphilly.csv`, `data/Y_LRdensity_westphilly.csv`, `data/Mapping_westphilly.csv` and `data/InvMapping_westphilly.csv`).

For convenience the compiled datasets `crime_tracts2019.rdata` and `crime_agr2019.rdata` have been saved in the repository, in `data/crime/`. To create these files, the data from the Philadelphia Police Department was accessed in February 2019.
However, they can be recreated using `get_crime.R` (which calls some helper functions in `subsetup/clean_crime.R`). This requires downloading the raw data from the Philadelphia Police Department and saving it in `data/csv/` (see below for details). 
Similarly, the shapefiles `phillytracts` and `phillyblockgroup` are saved as RData files in `data/shapefile/`, but can be recreated using the instructions below.

#### Instructions to create shapefiles, crime_agr2019 and crime_tracts2019 from the raw data

These steps are not strictly necessary, as I have saved the needed output in `data/`.

1. Download the crime CSV datasets from https://www.opendataphilly.org/dataset/crime-incidents for each year from 2006 to 2018. Combine the CSV files for all the years into one CSV file (make sure you don't copy the header multiple times) named `incidents_part1_part2.csv` and save it in `data/csv`. Note: I was able to download all the incidents from 2006 to 2018 (included) as one csv file from <a href="https://phl.carto.com/api/v2/sql?filename=incidents_part1_part2&format=csv&skipfields=cartodb_id,the_geom,the_geom_webmercator&q=SELECT%20*%20,%20ST_Y(the_geom)%20AS%20lat,%20ST_X(the_geom)%20AS%20lng%20FROM%20incidents_part1_part2%20WHERE%20dispatch_date_time%20%3E=%20%272006-01-01%27%20AND%20dispatch_date_time%20%3C%20%272019-01-01%27">this link</a>.

2. Download the shapefile for Philadelphia's census tracts from Open Data Philly for [census tracts](https://www.opendataphilly.org/dataset/census-tracts) and [block groups](https://www.opendataphilly.org/dataset/census-block-groups) and save the compressed folder in `data/shapefiles` and unzip it there, so that the shapefiles are in `data/shapefiles/tracts_pc2010/` (for example).  Note that the folder will contain a .shp file and other files with the same name but different extensions. Do not delete or change the name to those files. Also, note that I used the shapefiles from 2010.

3. Import the shapefiles in R using the R function `rgdal::readOGR` and save as Rdata files.

4. Run the script `get_data/get_crime.R`. This will create `crime_tract2019.rdata` and `crime_agr2019.rdata`, using the functions contained in `subsetup/clean_crime.R`.

Some of this code was created with the help of [Colman Humphrey](https://github.com/ColmanHumphrey/). The scripts here have been adapted from another of my reposotories [Urban-project](https://github.com/cecilia-balocchi/Urban-project).
