load('data/crime/crime_agr2019.rdata') # for block groups
load('data/crime/crime_tract2019.rdata') # for tracts
load("data/shapefile/phillyblockgroup")
load("data/shapefile/phillytracts")
tmp <- load("data/westphilly.rdata")

library(rgdal)
library(dplyr)

tmp <- rgdal::readOGR("data/shapefile/Census_Block_Groups_2010/Census_Block_Groups_2010.shp")
# I checked and ALAND10 from tmp, if aggregated, gives the ALAND10 from tracts

crime_agr2019 <- merge(crime_agr2019, tmp@data[,c("GEOID10","ALAND10")], by = "GEOID10")
crime_tract2019 <- merge(crime_tract2019, tracts@data[,c("GEOID10","ALAND10")], by = "GEOID10")

crime_agr2019$crime_density <- crime_agr2019$violent/crime_agr2019$ALAND10          # NOTE: ALAND10 is in square meters!
crime_tract2019$crime_density <- crime_tract2019$violent/crime_tract2019$ALAND10    # NOTE: ALAND10 is in square meters!

crime_tract2019_restricted <- crime_tract2019[!is.na(crime_tract2019$crime_density),] # Nothing is NA
crime_tract2019agr <- crime_tract2019_restricted %>% group_by(X) %>% summarize(GEOID10 = unique(GEOID10),
                                                      crime_density = mean(crime_density), 
                                                      violent = mean(violent))

crime_agr2019_restricted <- crime_agr2019[!is.na(crime_agr2019$crime_density),]
crime_agr2019agr <- crime_agr2019_restricted %>% group_by(X) %>% summarize(GEOID10 = unique(GEOID10),
                                                                           crime_density = mean(crime_density),
                                                                           violent = mean(violent))

## block groups
index <- crime_agr2019agr$GEOID10 %in% wp_GEOID_bgroups
crime_agr2019wp <- crime_agr2019agr[index,]

tmp <- match(bgroup_wp_index,as.numeric(crime_agr2019wp$X))
# all.equal(bgroup_wp_index,as.numeric(crime_agr2019wp$X[tmp])) # true
yHR_density <- as.numeric(crime_agr2019wp[tmp,"crime_density"][[1]])
yHR_violent <- as.numeric(crime_agr2019wp[tmp,"violent"][[1]])

## tracts
index <- crime_tract2019agr$GEOID10 %in% wp_GEOID_tracts
crime_tract2019wp <- crime_tract2019agr[index,]

tmp <- match(tracts_wp_index,as.numeric(crime_tract2019wp$X))
# all.equal(tracts_wp_index,as.numeric(crime_tract2019wp$X[tmp])) # true
yLR_density <- as.numeric(crime_tract2019wp[tmp,"crime_density"][[1]])
yLR_violent <- as.numeric(crime_tract2019wp[tmp,"violent"][[1]])

length(yLR_density) #56
length(yHR_density) #202

str <- "density_westphilly.csv"
write.table(x = yHR_density, row.names = FALSE, col.names = FALSE, sep=",", file = paste0("data/Y",str))
write.table(x = yLR_density, row.names = FALSE, col.names = FALSE, sep=",", file = paste0("data/Y_LR",str))

### Save Mapping and InvMapping

# ordered by census tract, contains the indices (name, mainID) of the block groups contained in each census tract
mapping <- matrix(NA, nrow = length(bgroup_wp_index), ncol = 2)
index <- 1
for(i in 1:length(tracts_wp_index)){
  rest <- as.numeric(bg_in_tr_list[[i]])
  for(j in 1:length(rest)){
    localID = which(bgroup_wp_index == rest[j])
    mapping[index,1] <- localID-1 
    mapping[index,2] <- i-1 
    index <- index+1
  }
}
# ordered by the names (mainID) contains the index of Mapping
inv_mapping <- matrix(NA, nrow = length(bgroup_wp_index), ncol = 2)
for(i in 1:nrow(mapping)){
  inv_mapping[i,1] <- which(mapping[,1]==i-1)-1
  inv_mapping[i,2] <- mapping[which(mapping[,1]==i-1),2] 
}

str <- "_westphilly.csv"
write.table(mapping, row.names = FALSE, col.names = FALSE, sep=",", file = paste0("data/Mapping",str))
write.table(inv_mapping, row.names = FALSE, col.names = FALSE, sep=",", file = paste0("data/InvMapping",str))

