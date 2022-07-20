library(dplyr)
library(rgeos)
load("data/shapefile/phillyblockgroup")
load("data/shapefile/phillytracts")
load("data/shapefile/phillyblock")

n_b <- length(phillyblock)
n_bg <- length(phillyblockgroup)
n_tr <- length(tracts)

xy <- array(NA, dim = c(n_bg,2))
for(i in 1:n_bg)
  xy[i,] <- phillyblockgroup@polygons[[i]]@labpt

xytr <- array(NA, dim = c(n_tr,2))
for(i in 1:n_tr)
  xytr[i,] <- tracts@polygons[[i]]@labpt

xyb <- array(NA, dim = c(n_b,2))
for(i in 1:n_b)
  xyb[i,] <- phillyblock@polygons[[i]]@labpt

west_philly_tract <- which((xytr[,1]<=-75.19) & (xytr[,2]>39.941)  & (xytr[,2]<40.01))

wp_GEOID_tracts <- tracts$GEOID10[west_philly_tract]
wp_GEOID_tracts <- levels(wp_GEOID_tracts)[wp_GEOID_tracts]

l <- list()
index <- 1
for(x in west_philly_tract){
  l[[index]] <- tracts@polygons[[x]]
  index <- index + 1
}
SP <- SpatialPolygons(l)
tracts_wp_index <- west_philly_tract
tracts_wp_name <- tracts_wp_index-1
tracts_wp_poly <- SpatialPolygonsDataFrame(Sr = SP, data = data.frame(index = tracts_wp_index, 
                                                                      row.names = tracts_wp_name,
                                                                      GEOID10 = tracts$GEOID10[west_philly_tract]))

sp <- SpatialPoints(xy)
bg_tr_cc <- gContains(spgeom1 = tracts_wp_poly, sp, byid = TRUE)
tmp <- which(bg_tr_cc, arr.ind = T) 
# everything is 1-indexed
west_philly_bgroup <- as.numeric(tmp[,1])
# all.equal(unique(tmp[,1]), as.numeric(tmp[,1])) # TRUE
bgroup_wp_index <- west_philly_bgroup

wp_GEOID_bgroups <- phillyblockgroup$GEOID10[bgroup_wp_index]
wp_GEOID_bgroups <- levels(wp_GEOID_bgroups)[wp_GEOID_bgroups]

l <- list()
index <- 1
for(x in bgroup_wp_index){
  l[[index]] <- phillyblockgroup@polygons[[x]]
  index <- index + 1
}
SP <- SpatialPolygons(l)
bgroup_wp_poly <- SpatialPolygonsDataFrame(Sr = SP, data = data.frame(index = bgroup_wp_index, 
                                                                      row.names = bgroup_wp_index,
                                                                      GEOID10 = phillyblockgroup$GEOID10[bgroup_wp_index]))

sp <- SpatialPoints(xyb)
b_tr_cc <- gContains(spgeom1 = tracts_wp_poly, sp, byid = TRUE)
b_bg_cc <- gContains(spgeom1 = bgroup_wp_poly, sp, byid = TRUE)
tmp1 <- which(b_tr_cc, arr.ind = T); tmp <- tmp1
# tmp2 <- which(b_bg_cc, arr.ind = T) 
# all.equal(sort(tmp1[,1]),sort(tmp2[,1])) #TRUE
## everything is 1-indexed
west_philly_block <- as.numeric(tmp[,1])
# all.equal(unique(tmp[,1]), as.numeric(tmp[,1])) # TRUE
block_wp_index <- west_philly_block
block_wp_name <- block_wp_index-1

wp_GEOID_blocks <- phillyblock$GEOID10[block_wp_index]
wp_GEOID_blocks <- levels(wp_GEOID_blocks)[wp_GEOID_blocks]

l <- list()
index <- 1
for(x in block_wp_index){
  l[[index]] <- phillyblock@polygons[[x]]
  index <- index + 1
}
SP <- SpatialPolygons(l)
block_wp_poly <- SpatialPolygonsDataFrame(Sr = SP, data = data.frame(index = block_wp_index, 
                                                                      row.names = block_wp_name,
                                                                      GEOID10 = phillyblock$GEOID10[block_wp_index]))


pbgwp <- fortify(bgroup_wp_poly)
pctwp <- fortify(tracts_wp_poly)
pbwp <- fortify(block_wp_poly)


bg_in_tr_list <- apply(bg_tr_cc, MARGIN = 2, which)
b_in_bg_list <- apply(b_bg_cc, MARGIN = 2, which)
b_in_tr_list <- apply(b_tr_cc, MARGIN = 2, which)

varlist <- c("tracts_wp_index","tracts_wp_name","tracts_wp_poly",
             "bgroup_wp_index","bgroup_wp_poly",
             "block_wp_index", "block_wp_name", "block_wp_poly",
             "wp_GEOID_tracts", "wp_GEOID_bgroups","wp_GEOID_blocks",
             "bg_in_tr_list", "b_in_tr_list", "b_in_bg_list")
save(list = varlist, file = "data/westphilly.rdata")


# library(ggplot2)
# p1 <- ggplot() +  geom_polygon(data=pbgwp, aes(x=long,y=lat,group=group),
#                                color=alpha("red",0.2), alpha =0.1) + coord_map() 
# p1 <- p1+  geom_polygon(data=pctwp, aes(x=long,y=lat,group=group),
#                color=alpha("black",0.3), alpha =0)
# p1 <- p1+  geom_polygon(data=pbwp, aes(x=long,y=lat,group=group),
#                         color=alpha("blue",0.1), alpha =0)
  
# p1
