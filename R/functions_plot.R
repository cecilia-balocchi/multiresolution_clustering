require(rgeos)
require(reshape2)
require(ggplot2)

get_ids_HR <- function(idsLR, tmp_array){
  zHR <- numeric(dim(tmp_array)[1])
  for(cl in unique(idsLR)){
    indexLR <- which(idsLR == cl)
    indexHR <- which(tmp_array[,4] %in% indexLR)
    zHR[indexHR] <- cl
  }
  return(zHR)
}
get_ids_HR_col <- function(idsLR, tmp_array_col){
  zHR <- numeric(length(tmp_array_col))
  for(cl in unique(idsLR)){
    indexLR <- which(idsLR == cl)
    indexHR <- which(tmp_array_col %in% indexLR)
    zHR[indexHR] <- cl
  }
  return(zHR)
}

get_ids_from_pairwise <- function(pairwise){
  N <- dim(pairwise)[1]
  visited <- c()
  clusters <- list()
  for(i in 1:N){
    if(!(i %in% visited)){
      samecl <- which(pairwise[i,] == 1)
      samecl <- c(i, samecl)
      visited <- c(visited, samecl)
      clusters <- c(clusters, list(samecl))
    }
  }
  clusters_id <- rep(NA, N)
  id <- 1
  for(xx in clusters){
    clusters_id[xx] <- id
    id <- id + 1
  }
  return(clusters_id)
}

# functions for plotting

# since we are working with a subset of the data, I need to separate my subset (index) 
# with the outside (index_outside)
border_list <- function(clusters_id, index, index_outside, w){
  # within my subset, I find which pairs of indices belong to different clusters
  # if i,j belong to different clusters, a will contain an element (i,j)
  # moreover, if i shares a border with someone outside, I will add that pair as well
  a <- list()
  windex <- w[index, index]
  tmp <- 1
  for(i in 1:length(index)){
    js <- which(windex[i,]==1) 
    for(j in js){
      if(clusters_id[i] != clusters_id[j]){
        if(j > i){
          a[[tmp]] <- c(index[i],index[j])
          tmp <- tmp + 1
        }
      }
    }
    if(!is.null(index_outside)){
      ind_outside <- which(w[index[i], index_outside]==1)
      if(length(ind_outside)>0){
        for(x in ind_outside){
          a[[tmp]] <- c(index[i],index_outside[x])
          tmp <- tmp + 1
        }
      }
    }
  }
  a
}

border_SL_list <- function(border_list, phillyblockgroup){
  # using the border_list from before, I need to get the intersection of that pair
  # if the intersection is a Spatial Line, I create a list of the lines
  l <- list()
  iter <- 1
  if(length(border_list) == 0){
    return(list())
  }
  for(n in 1:length(border_list)){
    i <- border_list[[n]][1]
    j <- border_list[[n]][2]
    # cat(n,":",i,j)
    A <- SpatialPolygons(list(Polygons(list(phillyblockgroup@polygons[[i]]@Polygons[[1]]),ID=1)))
    B <- SpatialPolygons(list(Polygons(list(phillyblockgroup@polygons[[j]]@Polygons[[1]]),ID=1)))
    tmp <- gIntersection(A,B)
    if(class(tmp)[1] == "SpatialLines"){
      tmp1 <- tmp@lines[[1]]
      tmp1@ID <- as.character(iter)
      l[[iter]] <- tmp1
      iter <- iter + 1
      # cat("\n")
    } else {
      if(class(tmp)[1] == "SpatialCollections"){
        if(.hasSlot(tmp, "lineobj")){
          for(ni in 1:length(tmp@lineobj@lines)){
            tmp1 <- tmp@lineobj@lines[[ni]]
            tmp1@ID <- as.character(iter)
            l[[iter]] <- tmp1
            iter <- iter + 1
          }
          # cat("\n")
        }
      } else {
        # cat(" error!\n")
      }
    }
  }
  l
}
# this is used to draw the external border for regions that are at the external part of the city
external_border_list <- function(l, array){ 
  # I will extend the list l by adding some spatial lines constructed from an array of points
  #id is the starting id I can use
  id <- length(l) + 1
  n <- dim(array)[1]-1
  l1 <- list()
  for(i in 1:n){
    l1[[i]] <- Line(array[c(i,i+1),])
  }
  l[[id]]<- Lines(l1, ID = as.character(id))
  l
}

get_SLDF <- function(z, index, index_outside, w.sym, shp, external = FALSE, ext_line = NULL){
  bl <- border_list(clusters_id = z, index = index, index_outside = index_outside, w = w.sym)
  SL1 <- border_SL_list(bl, shp)
  if(external){
    if(is.null(ext_line)){
      warning("ext_line needed")
    } else {
      SL1 <- external_border_list(SL1, ext_line)
    }
  }
  SL_complete <- SpatialLines(SL1)
  mySLDF <- SpatialLinesDataFrame(SL_complete, data = data.frame(ID = rep(1,length(SL1)),row.names = 1:length(SL1)))
  return(mySLDF)
}

plot_pairwise <- function(pw, z, index = NULL){
  n <- length(z)
  if(is.null(index)){
    index2 <- sort(z, index.return=T)$ix
    ticks <- c(0,cumsum(table(z))) + 0.5
  } else {
    index2 <- index 
    table2 <- tabulate(z)[unique(z[index2])]
    ticks <- c(0,cumsum(table2)) + 0.5
  }
  mm <- melt(pw[index2,index2])
  p <- ggplot(mm) + geom_tile(aes(x=Var1, y=Var2,fill = value),
                              colour = "white") + scale_fill_gradient(low = "white", high = "steelblue")
  for(x in ticks){
     p <- p + geom_vline(xintercept = x) + geom_hline(yintercept = x)
  }
  
  p
}



