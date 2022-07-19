source("R/functions_analysis.R")
source("R/paralleltemp.R")
str <- "density_westphilly.csv"

## this can be imported from command line script, or set manually to be equal to 1 or 2
## TASK = 3,4 would run the algorithm with the k-means partition as starting point.
TASK <- strtoi(Sys.getenv("SGE_TASK_ID"))

set.seed(943)
seed1 = sample(1000000,1)
seed2 = sample(1000000,1)
seed3 = sample(1000000,1)
seed4 = sample(1000000,1)

# Import data and standardize
Y <- read.csv(paste0("data/Y",str), header = F)
Y <- as.vector(Y[[1]])
YLR <- read.csv(paste0("data/Y_LR",str), header = F)
YLR <- as.vector(YLR[[1]])
str <- "_westphilly.csv"
Mapping_mat <- read.csv(paste0("data/Mapping",str), header = F)
Mapping <- as.vector(as.matrix(Mapping_mat))
InvMapping_mat <- read.csv(paste0("data/InvMapping",str), header = F)
InvMapping <- as.vector(as.matrix(InvMapping_mat))

Y <- scale(Y)
meanHR <- attributes(Y)$`scaled:center`
sdHR <- attributes(Y)$`scaled:scale`
Y <- Y[,1]
YLR <- (YLR - meanHR)/sdHR  # actually not used

n <- length(Y)
n2 <- length(YLR)

## setting parameters
burnin <- 0
adjust <- 0 
skip <- 10
save <- 50000
save_nDP <- save/skip
gibbs <- 3

# within cluster sd is 0.25
m <- 0.5^2
v <- 0.1
alpha <- 2 + m^2/v
beta <- m*(1 + m^2/v)

kH = 1/10
alphaH = alpha
betaH = beta

kL = 1/10 
alphaL = alpha
betaL = beta/3 # sort of average number of bg in ct

eta_LR <- 2
eta_CT <- 0.5
eta_TD <- 2

alpha2 <- 0.0

nchains <- 20
temps <- seq(0,1,length = 1+nchains)[-1]

## find K-means partition for initialization in TASKs 3,4
d <- dist(Y)
Kmax <- 10
nstarts <- 25
km_silhouette <- numeric(Kmax)
for(k in 2:min(Kmax, n-1)){
  km = kmeans(Y, centers = k, nstart = nstarts)
  tmp_sil <- silhouette(x = km$cluster, dist = d) # not defined for k = 1
  km_silhouette[k] <- mean(tmp_sil[, 3])
}

k_hr <- which.max(km_silhouette)
km = kmeans(Y, centers = k_hr, nstart = nstarts)
z_hr <- km$cluster 

km2_silhouette <- numeric(Kmax)
# find proportions of HR clusters in each LR unit, use as data for second km
prop_hr <- array(0, dim = c(n2, k_hr))
for(i in 1:n2){
  # index <- tmp_array[tmp_array[,2] == i,1] # equal to: Mapping_mat[(Mapping_mat[,2] == i - 1),1]+1
  index <- Mapping_mat[(Mapping_mat[,2] == i - 1),1]+1
  cl_index <- as.numeric(names(table(z_hr[index])))
  prop_hr[i,cl_index] <- as.numeric(table(z_hr[index]))
  prop_hr[i,] <- prop_hr[i,]/sum(prop_hr[i,])
}
d2 = dist(prop_hr)

for(k in 2:min(Kmax, n2-1,nrow(unique(prop_hr)))){
  km2 = kmeans(prop_hr, centers = k, nstart = nstarts)
  tmp_sil <- silhouette(x = km2$cluster, dist = d2)
  km2_silhouette[k] <- mean(tmp_sil[, 3])
}

k_lr <- which.max(km2_silhouette)
km2 = kmeans(prop_hr, centers = k_lr, nstart = nstarts)
z_lr <- km2$cluster

## These are the ones we used in the paper:
## these two start from everyone alone or together
if(TASK == 1){
  out_nHDP <- paralleltemp(Y, YLR = NULL, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = FALSE, one_vs_manyHR_bool = FALSE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = kH, alphaH = alphaH, betaH = betaH, 
                            seed = seed1, alpha1 = 1.0, alpha2 = 0.0)
}
if(TASK == 2){
  out_nHDP <- paralleltemp(Y, YLR = NULL, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = FALSE, one_vs_manyHR_bool = TRUE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = kH, alphaH = alphaH, betaH = betaH, 
                            seed = seed2, alpha1 = 1.0, alpha2 = 0.0)
}

## Optional, not used in paper:
## these two start from kmeans
if(TASK == 3){
  out_nHDP <- paralleltemp(Y, YLR = NULL, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            gammaLR_init = z_lr, gammaHR_init = z_hr,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = kH, alphaH = alphaH, betaH = betaH, 
                            seed = seed3, alpha1 = 1.0, alpha2 = 0.0)
}
if(TASK == 4){
  out_nHDP <- paralleltemp(Y, YLR = NULL, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            gammaLR_init = z_lr, gammaHR_init = z_hr,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = kH, alphaH = alphaH, betaH = betaH, 
                            seed = seed4, alpha1 = 1.0, alpha2 = 0.0)
}

filename <- paste0("results/westphilly/wp_density_new_burn",burnin,"adj",adjust,"save",save,"gibbs",gibbs,
                   "k",kH,"_m",m,"_v",v,"_eta",eta_LR,"_",eta_CT,"_",eta_TD,
                   "_scaled_nchains",nchains,"_task",TASK,".rdata")
varlist <- c("out_nHDP",
             "kH","alphaH","betaH",
             "kL","alphaL","betaL",
             "eta_LR","eta_CT","eta_TD", "alpha2",
             "burnin", "adjust", "save", "skip", "save_nDP", "gibbs",
             "nchains", "temps")
save(list = varlist, file = filename)                      
