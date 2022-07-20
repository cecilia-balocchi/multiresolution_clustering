### This file creates the object "crime_agr2019.rdata" and "crime_tract2019.rdata" from scratch
### For convenience those files have already been saved in the folder data/crime/
### So this script can be skipped.

# Note the crime dataset (CSV) downloaded from the Philly Police Dept 
# (https://www.opendataphilly.org/dataset/crime-incidents)
# should be saved in `data/csv/` as `incidents_part1_part2.csv`
# The sub-directory `data/shapefile` contains the geographic files 
# for aggregating over the neighborhoods (these shapefiles where downloaded from 
# opendataphilly.org, read into R using 'readOGR' and saved as Rdata objects)

source("get_data/subsetup/clean_crime.R")
# crimes.new is a dataset where each row represents a crime and we save time/location and crimetype info
crimes.new <- clean_crime.from_csv("csv/incidents_part1_part2.csv", data_path = "data/")
save(crimes.new, file = 'data/crime/crimes.new2019.rdata')

#### Let's start from census tracts ####

# Now we restrict to using only the crimes up to (including) year 2018 (so only using 13 years, i.e. columns).
# in crime_agr each line gives the number of crimes for a region and a year.
crime_tract2019 <- aggregate_crime.from_rdata("data/crime/crimes.new2019.rdata", ncols = 13, variable = "tract")
crime_tract2019$tr.violent <- ihst(crime_tract2019$violent)    # Inverse Hyperbolic Sine Transformation (not used)
crime_tract2019$log.violent <- log(crime_tract2019$violent+1)  # (not used)
crime_tract2019$year19 <- crime_tract2019$year-2005
save(crime_tract2019, file = "data/crime/crime_tract2019.rdata")

#### Now we do the same for census block groups ####
crime_agr2019 <- aggregate_crime.from_rdata("data/crime/crimes.new2019.rdata", ncols = 13, variable = "blockgroup")
crime_agr2019$tr.violent <- ihst(crime_agr2019$violent)         # Inverse Hyperbolic Sine Transformation (not used)
crime_agr2019$log.violent <- log(crime_agr2019$violent+1)       # (not used)
crime_agr2019$year19 <- crime_agr2019$year-2005
save(crime_agr2019, file = "data/crime/crime_agr2019.rdata")