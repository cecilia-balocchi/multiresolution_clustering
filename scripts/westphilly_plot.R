library(mcclust.ext)
library(ggplot2)
library(spdep)
library(rgeos)
library(dplyr)
library(gridExtra)
library(mapproj)
source("R/functions_plot.R")
source("R/functions_analysis.R")
load("data/westphilly.rdata")

# Prepare spatial data
forttract <- fortify(tracts_wp_poly)
fortbgroup <- fortify(bgroup_wp_poly)
nLR <- length(tracts_wp_poly)
nHR <- length(bgroup_wp_poly)


list.poly <- poly2nb(bgroup_wp_poly)
w <- matrix(data = 0, nrow = nHR, ncol = nHR)
for(i in 1:nHR){
  w[i, list.poly[[i]]] <- 1
}
w.sym <- w + t(w) - w*t(w)

list.poly <- poly2nb(tracts_wp_poly)
wLR <- matrix(data = 0, nrow = nLR, ncol = nLR)
for(i in 1:nLR){
  wLR[i, list.poly[[i]]] <- 1
}
wLR.sym <- wLR + t(wLR) - wLR*t(wLR)


################## Plot the raw data ##################

# Import Y data
str <- "density_westphilly.csv"
Y <- read.csv(paste0("data/Y",str), header = F)
Y <- as.vector(Y[[1]])
YLR <- read.csv(paste0("data/Y_LR",str), header = F)
YLR <- as.vector(YLR[[1]])
str <- "_westphilly.csv"
Mapping_mat <- read.csv(paste0("data/Mapping",str), header = F)
Mapping <- as.vector(as.matrix(Mapping_mat))
InvMapping_mat <- read.csv(paste0("data/InvMapping",str), header = F)
InvMapping <- as.vector(as.matrix(InvMapping_mat))
Mapping_mat = Mapping_mat + 1
n_hr <- length(Y)
n_lr <- length(YLR)

Y <- Y * 1e+6
YLR <- YLR * 1e+6

Var <- c(Y, YLR); var <- c(Y); varLR <- c(YLR);# varbl <- c(Ybl)
limits <- c(min(Var, na.rm=T), max(Var, na.rm=T))
MI <- limits[1]; MA <- limits[2]; 
miHR <- min(var, na.rm=T); maHR <- max(var, na.rm=T)
zHR <- (0-miHR)/(maHR-miHR); zHR_adj <-(maHR-miHR)/(MA-MI) * zHR + (miHR-MI)/(MA-MI)
mi <- min(varLR, na.rm=T); ma <- max(varLR, na.rm=T)
zLR <- (0-mi)/(ma-mi); zLR_adj <-(ma-mi)/(MA-MI) * zHR + (mi-MI)/(MA-MI)


## Plot of real data, with tract borders and block groups borders
my.data1 <- data.frame(id = bgroup_wp_index, Y = as.numeric(Y))
my.data1$id <- as.character(my.data1$id)
plotData_real <- left_join(fortbgroup,my.data1, by = "id")
p3 <- ggplot() + geom_polygon(data=plotData_real, aes(x=long,y=lat,group=group, fill = Y),
                              color=alpha("red",0.3),size = 0.2, alpha =0.8) + coord_map() +
  scale_fill_distiller(palette = "Spectral", limits = c(MI,MA), values = c(0,zHR_adj,1))+ labs(fill = "Crime\ndensity")#labs(fill = "Log\ncrime")
p3 <- p3 + geom_polygon(data=forttract, aes(x=long,y=lat,group=group),
                        color=alpha("blue",0.8), size = 0.3, alpha =0) + 
  theme(panel.background = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5)) + ggtitle("High-resolution data")
p3

my.data1 <- data.frame(id = tracts_wp_name, Y = as.numeric(YLR))
my.data1$id <- as.character(my.data1$id)
plotData_real <- left_join(forttract,my.data1, by = "id")
p3b <- ggplot() + geom_polygon(data=plotData_real, aes(x=long,y=lat,group=group, fill = Y),
                               color=alpha("blue",0.8),size = 0.3, alpha =0.8) + coord_map() +
  scale_fill_distiller(palette = "Spectral", limits = c(MI,MA), values = c(0,zLR_adj,1))+ labs(fill = "Crime\ndensity")+#labs(fill = "Log\ncrime")
  theme(panel.background = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5)) + ggtitle("Low-resolution data")
p3b

ps <- list(p3b,p3)
filename_string <- "figures/realdata_density_wp_scaled.png"
ggsave(filename = filename_string, plot = arrangeGrob(grobs = ps, nrow = 1),
       device = "png", width = 5, height = 2, units = "cm", scale = 5) #path = "partplots",



################## Plot the estimated partitions ##################

## fixed eta
## m = 0.5^2, v = 0.1
temp <- load(paste0("results/westphilly/wp_density_new_burn5000adj0save50000gibbs3k0.1_m0.25_v0.1_eta2_0.5_2_scaled_nchains20_task1.rdata")) 
out1 <- out_nHDP
temp <- load(paste0("results/westphilly/wp_density_new_burn5000adj0save50000gibbs3k0.1_m0.25_v0.1_eta2_0.5_2_scaled_nchains20_task2.rdata")) 
out2 <- out_nHDP

burnin <- 10000
skip <- 50
keep <- seq(burnin, save, by = skip)

out <- list()
out$highres <- rbind(out1$highres[keep,], out2$highres[keep,])
out$lowres <- rbind(out1$lowres[keep,], out2$lowres[keep,])
out$highres_mean <- rbind(out1$highres_mean[keep,], out2$highres_mean[keep,])
postmean <- colMeans(out$highres_mean)

### Now we find and plot the low resolution clusters

tmp_hier <- get_partition(out$lowres, hier = T, Kmode = T)
tmp_wade <- get_partition(out$lowres, hier = F)
pw_lr <- tmp_hier$pw; z_lr_wade <- tmp_wade$z; z_lr_hier <- tmp_hier$z
k_wade <- length(unique(z_lr_wade)); z_lr_hier2 <- get_partition(out$lowres, hier = T, K = k_wade)$z

vi_hier <- VI(matrix(z_lr_hier, nrow = 1), out$lowres)
vi_hier2 <- VI(matrix(z_lr_hier2, nrow = 1), out$lowres)
vi_wade <- VI(matrix(z_lr_wade, nrow = 1), out$lowres)

vi_min <- which.min(c(vi_hier, vi_hier2, vi_wade)); vi_min
if(vi_min == 1){
  cat("hier has smaller VI\n")
  z_lr <- z_lr_hier
} else if(vi_min == 2){
  cat("hier2 has smaller VI\n")
  z_lr <- z_lr_hier2
} else {
  z_lr <- z_lr_wade
}

mySLDF_lr <- get_SLDF(z = z_lr, index = 1:n_lr, index_outside = NULL,
                      w.sym = wLR.sym, shp = tracts_wp_poly)

Ycl_lr <- as.numeric(YLR)
for(zk in unique(z_lr)){
  ind <- which(z_lr == zk)
  Ycl_lr[ind] <-mean(Ycl_lr[ind])
}


my.data2 <- data.frame(id =tracts_wp_name, cluster = as.factor(z_lr), value = Ycl_lr)
my.data2$id <- as.character(my.data2$id)
plotData_part_lr <- left_join(forttract,my.data2, by = "id")

p1 <- ggplot() +  geom_polygon(data=plotData_part_lr, aes(x=long,y=lat,group=group, fill = value),
                               color=alpha("blue", 0.5),size = 0.3, alpha =0.8) + coord_map() +
  scale_fill_distiller(palette = "Spectral", name = "Cluster\naverage",limits = c(MI,MA), values = c(0,zLR_adj,1))+
  theme(panel.background = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste0("Low-resolution clusters (K = ",length(unique(z_lr)),")"))
p1 <- p1 + geom_line(data = mySLDF_lr, aes(x=long,y=lat,group=group),size = 0.8, color=alpha("black",1))
p1

### Now we find and plot the high resolution clusters

tmp_hier <- get_partition(out$highres, hier = T,  Kmode = T)
tmp_wade <- get_partition(out$highres, hier = F)
pw_hr <- tmp_wade$pw; z_hr_wade <- tmp_wade$z; z_hr_hier <- tmp_hier$z
k_wade <- length(unique(z_hr_wade)); z_hr_hier2 <- get_partition(out$highres, hier = T, K = k_wade)$z

vi_hier <- VI(matrix(z_hr_hier, nrow = 1), out$highres)
vi_hier2 <- VI(matrix(z_hr_hier2, nrow = 1), out$highres)
vi_wade <- VI(matrix(z_hr_wade, nrow = 1), out$highres)

vi_min <- which.min(c(vi_hier, vi_hier2, vi_wade)); vi_min
if(vi_min == 1){
  cat("hier has smaller VI\n")
  z_hr <- z_hr_hier
} else if(vi_min == 2){
  cat("hier2 has smaller VI\n")
  z_hr <- z_hr_hier2
} else {
  z_hr <- z_hr_wade
}

mySLDF_hr <- get_SLDF(z = z_hr, index = 1:n_hr, index_outside = NULL,
                      w.sym = w.sym, shp = bgroup_wp_poly)

Ycl_hr <- as.numeric(Y)
for(zk in unique(z_hr)){
  ind <- which(z_hr == zk)
  Ycl_hr[ind] <-mean(Ycl_hr[ind])
}

my.data <- data.frame(id = bgroup_wp_index, value =Y, value_cl = Ycl_hr,
                      cluster = as.factor(z_hr), pmean = postmean)
my.data$id <- as.character(my.data$id)
plotData <- left_join(fortbgroup,my.data, by = "id")
p2 <- ggplot() +  geom_polygon(data=plotData, aes(x=long,y=lat,group=group, fill = value_cl),
                               color=alpha("red", 0.5),size = 0.2, alpha =0.8) +
  scale_fill_distiller(type = "div", palette = 9, limits = c(MI,MA), values = c(0,zHR_adj,1),
                       name = "Cluster\naverage") + coord_map()+
  geom_polygon(data=forttract, aes(x=long,y=lat,group=group),color=alpha("blue", 0.5),size = 0.3, alpha =0) +
  theme(panel.background = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste0("High-resolution clusters (K = ",length(unique(z_hr)),")"))
p2 <- p2 + geom_line(data = mySLDF_hr, aes(x=long,y=lat,group=group),size = 0.8, color=alpha("black",1))
p2

str_input <- "Density_nHDP_m025_V01_fixed_etas_2_05_2"
ps <- list(p1,p2)
# do.call("grid.arrange", c(ps, nrow = 1))          # to display in R
filename_string <- paste0("figures/",str_input,".png")
ggsave(filename = filename_string, plot = arrangeGrob(grobs = ps, nrow = 1),
       device = "png", width = 5, height = 2, units = "cm", scale = 5) #path = "partplots",

