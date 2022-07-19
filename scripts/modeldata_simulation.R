library(mcclust.ext)
library(mcclust)
library(spdep)

source("R/functions_plot.R")
source("R/functions_analysis.R")

source("R/paralleltemp.R")
source("R/paralleltemp_noLR.R")
source("R/nDP.R")

### Parameters to run the script
# This script was run on a high performance cluster, 
# with the script being called from a command line script, including the line
# `R --no-save --args arg1 arg2 arg3 arg4 < scripts/modeldata_all.R`
# where the vector (arg1, arg2, arg3, arg4) is passed in the following as `new_args`.
# The first argument determines the simulation number and method
# the second determines the nHDP hyperparameters (config),
# the third fixes the number of HR units within LR ones, 
# the fourth determines the distance between clusters. 
# If running the simulation from RStudio, you can fix the values of the second, third and fourth arguments
# and have the first component range for the desired simulations and methods.

# In the reported simulations, the second argument (`config`) was fixed to either 1 or 2,
# the third argument (`opt_m`) was fixed to either 2 or 4,
# while the third argument (`opt_dist`) was fixed to 2, which is considered a level of moderate cluster separation.

new_args <- commandArgs(TRUE)
tmp <- as.numeric(new_args[1])
sim_number <- (tmp-1) %% 50 + 1
method <- (tmp-1) %/% 50 + 1
sim <- as.character(sim_number)
config <- as.character(new_args[2])
if(config == 1){
  etas_str <- "353"
} else if(config == 2){
  etas_str <- "111"
}
load(paste0("data/newsim/nCT10nBGinCT100_etas",etas_str,"_sim",sim,".rdata"))
opt_m <- as.numeric(new_args[3])
opt_dist <- as.numeric(new_args[4])

cat(sim_number, "\n")
cat(config, "\n")
cat(opt_m, "\n")
cat(opt_dist, "\n")
cat(method, "\n")

dists <- c(2, 2.5, 3, 4)
dists_str <- c("2", "25", "3", "4")

dist <- dists[opt_dist]

### Generation of hierarchical structure of HR (bg = block groups) and LR (ct = census tracts) units### Generation of hierarchical structure of HR (bg = block groups) and LR (ct = census tracts) units

n_ct = 10                        # number of low-res units
if(opt_m == 1){                  # number of high-res units in each low-res units
  n_bg_in_ct = 5 
} else if(opt_m == 2){
  n_bg_in_ct = 10
} else if(opt_m ==3){
  n_bg_in_ct = 25
} else if(opt_m ==4){
  n_bg_in_ct = 50
} else if(opt_m ==5){
  n_bg_in_ct = 100
}
n_bg <- n_ct * n_bg_in_ct      # total number of high-res units
bg_in_ct <- t(array(1:n_bg, dim = c(n_bg_in_ct, n_ct)))

tmp_array <- array(NA, dim = c(n_bg,2))
for(i in 1:n_ct){
  for(x in bg_in_ct[i,]){
    tmp_array[x, ] <- c(x, i)
  }
}

Mapping <- tmp_array-1
InvMapping <- array(NA, dim = c(n_bg,2))
for(i in 1:nrow(Mapping)){
  index <- which(Mapping[,1]==i-1)
  InvMapping[i,1] <- index-1
  InvMapping[i,2] <- Mapping[index,2]
}

nHR <- length(unique(tmp_array[,1]))
nLR <- length(unique(tmp_array[,2]))

### Create dataset that will store the results

# Possible methods considered
configs <- c("km_nHDP", "km", "nDP", 
             "nHDP1_0","nHDP1_0_sampleeta", 
             "nHDPheur", "HDPoracle", "HDPone", "HDPmany")
# Various measures to evaluate methods
colnames_results <- c("R2_est","R2_bayes","R2_orac",
                  "R2_bayes_mean", "R2_bayes_sd", 
                  "R2_OS_est", "R2_OS_bayes", "R2_OS_orac",
                  "MSE_est","MSE_bayes","MSE_orac",
                  "MSE_bayes_mean", "MSE_bayes_sd", 
                  "MSE_Gjs_est", "MSE_Gjs_bayes", "MSE_Gjs_orac", 
                  "VarY", 
                  "B","B_bayes","B_bayes_sd","B_bayes2","B_bayes_sd2",
                  "R", "Radj", "R_bayes","R_bayes_sd", "Radj_bayes", "Radj_bayes_sd", 
                  "VI_est", "VI_bayes", "VI_bayes_sd", 
                  "K_est", "K_bayes","K_bayes_sd","K_true",
                  "lprior_est","lprior_bayes","lprior_bayes_iqr")
result <- array(NA, dim = c(length(configs), length(colnames_results)))
colnames(result) <- colnames_results
rownames(result) <- configs
resultHR <- result
resultLR <- result
# To store running times
timing <- numeric(length(configs))
names(timing) <- configs
# To store partition estimates
zHR_est <- array(NA, dim = c(length(configs), n_bg))
zLR_est <- array(NA, dim = c(length(configs), n_ct))
rownames(zHR_est) <- configs
rownames(zLR_est) <- configs
# To store acceptance rates of MH algorithms 
accept <- array(NA, dim = c(length(configs),3))
rownames(accept) <- configs

### Data generating process: sample the partitions

set.seed(609 + sim_number)
seed1 <- sample(1000000, 1) 
seed2 <- sample(1000000, 1)

sigma <- 0.5

## this was computed by the generate_newsim_data script
zLR_HR <- get_ids_HR_col(zLR, tmp_array[,2])

rests <- list()
n_rests <- numeric(kLR)
zHR <- numeric(nHR)
for(h in 1:kLR){
  rests[[h]] <- which(zLR_HR == h)
  n_rests[h] <- length(rests[[h]])
  clLR <- which(zLR == h)
  for(j in clLR){
    ind <- tmp_array[tmp_array[,2] == j,1]
    if(kLR == 1){
      x <- rmultinom(length(ind), size = 1, prob = Gjs)
    } else {
      x <- rmultinom(length(ind), size = 1, prob = Gjs[h,])
    }
    z_j <- apply(x, 2, function(x) which(x == 1))
    zHR[ind] <- z_j
  }
}
zHR_ordered <- reorder_part(zHR)
Gjs_ordered <- Gjs[,unique(zHR)]

Gjs_unordered <- Gjs
zHR_unordered <- zHR
Gjs <- Gjs_ordered
zHR <- zHR_ordered

nD <- length(unique(zHR))

# This is the true proportion of HR units for each LR *cluster*
prop <- array(0, dim = c(kLR, nD))
for(h in 1:kLR){
  index <- as.numeric(names(table(zHR[rests[[h]]])))
  prop[h,index] <- as.numeric(table(zHR[rests[[h]]]))/n_rests[h]
}

# This is the true proportion of HR units for each LR *unit*
propHR <- array(0, dim = c(nLR, nD))
for(i in 1:nLR){
  index <- tmp_array[tmp_array[,2] == i,1] # equal to: Mapping_mat[(Mapping_mat[,2] == i - 1),1]+1
  cl_index <- as.numeric(names(table(zHR[index])))
  propHR[i,cl_index] <- as.numeric(table(zHR[index]))
  propHR[i,] <- propHR[i,]/sum(propHR[i,])
}
new_rest <- find_unique_dist(propHR, eps = 0.2)
zLR_emp <- new_rest
zLR_HR_emp <- get_ids_HR_col(zLR, tmp_array[,2])

### Data generating process: sample the HR data
yHR <- numeric(nHR)
yHR_outsample <- numeric(nHR)
mHR <- numeric(nHR)
limit <- dist*(nD-1)/2
means <- seq(-limit,limit,length = nD)
for(d in 1:nD){
  cl <- which(zHR == d)
  yHR[cl] <- means[d] + rnorm(length(cl), mean = 0, sd = sigma)
  mHR[cl] <- means[d]
}
yHR_outsample <- mHR + rnorm(nHR, mean = 0, sd = sigma)

Gjs_means <- rowSums(t(t(Gjs) * means))
Gjs_means <- Gjs_means[zLR]

yLR <- numeric(nLR)
yLR_outsample <- numeric(nLR)
mLR <- numeric(nLR)

for(i in 1:nLR){
  index <- tmp_array[tmp_array[,2] == i,1]
  yLR[i] <- mean(yHR[index])
  yLR_outsample[i] <- mean(yHR_outsample[index])
  mLR[i] <- mean(mHR[index])
}

yHRcl_true <- as.numeric(yHR)
for(zk in unique(zHR)){
  ind <- which(zHR == zk)
  yHRcl_true[ind] <-mean(yHRcl_true[ind])
}
yLRcl_true <- as.numeric(yLR)
for(zk in unique(zLR)){
  ind <- which(zLR == zk)
  yLRcl_true[ind] <-mean(yLRcl_true[ind])
}
Y <- yHR
n <- nHR
n2 <- nLR

dataset<-matrix(c(tmp_array[,2],Y),ncol=2)        # final data set

### Parameters for running the methods

burnin <- 2000 
adjust <- 0
skip <- 10
save <- 10000
save_nDP <- save/skip
gibbs <- 3

k0 = 0.01
alpha0 = 5
beta0 = 1

alpha0_LR = alpha0
beta0_LR = beta0/n_bg_in_ct

eta_LR <- alpha_LR
eta_CT <- alpha0_HR 
eta_TD <- gamma_HR  

nchains <- 10
temps <- seq(0,1,length = 1+nchains)[-1]

grid = seq(0.5, 20, by= 0.5)
eta_m <- 3; eta_sd <- 2
prior = dnorm(grid, mean = eta_m, sd = eta_sd)
prior = prior/sum(prior)

varlist <- c("burnin", "adjust", "skip", "save", "save_nDP", "gibbs",
  "k0", "alpha0", "beta0", "alpha0_LR", "beta0_LR", "eta_LR", "eta_CT", "eta_TD", "nchains", "temps",
  "grid", "eta_m", "eta_sd", "prior")

if(method == 1){
  library(cluster)

  ptm <- proc.time()

  d <- dist(Y)
  Kmax <- 20
  nstarts <- 25
  km_silhouette <- numeric(Kmax)
  km2_silhouette <- numeric(Kmax)
  for(k in 2:min(Kmax, nHR-1)){
    km = kmeans(Y, centers = k, nstart = nstarts)
    tmp_sil <- silhouette(x = km$cluster, dist = d) # not defined for k = 1
    km_silhouette[k] <- mean(tmp_sil[, 3])
  }
  
  k_hr <- which.max(km_silhouette)
  km = kmeans(Y, centers = k_hr, nstart = nstarts)
  z_hr <- km$cluster 

  # find proportions of HR clusters in each LR unit, use as data for second km
  prop_hr <- array(0, dim = c(nLR, k_hr))
  for(i in 1:nLR){
    index <- tmp_array[tmp_array[,2] == i,1] # equal to: Mapping_mat[(Mapping_mat[,2] == i - 1),1]+1
    cl_index <- as.numeric(names(table(z_hr[index])))
    prop_hr[i,cl_index] <- as.numeric(table(z_hr[index]))
    prop_hr[i,] <- prop_hr[i,]/sum(prop_hr[i,])
  }
  d2 = dist(prop_hr)

  for(k in 2:min(Kmax, nLR-1,nrow(unique(prop_hr)))){
    km2 = kmeans(prop_hr, centers = k, nstart = nstarts)
    tmp_sil <- silhouette(x = km2$cluster, dist = d2)
    km2_silhouette[k] <- mean(tmp_sil[, 3])
  }

  k_lr <- which.max(km2_silhouette)
  km2 = kmeans(prop_hr, centers = k_lr, nstart = nstarts)
  z_lr <- km2$cluster

  out_nHDP1 <- paralleltemp(Y, YLR = NULL, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = FALSE, one_vs_manyHR_bool = FALSE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = k0, alphaH = alpha0, betaH = beta0, 
                            kL = k0, alphaL = alpha0_LR, betaL = beta0_LR, 
                            gammaLR_init = z_lr, gammaHR_init = z_hr,
                            seed = seed1, alpha1 = 1.0, alpha2 = 0.0)
  out_nHDP2 <- paralleltemp(Y, YLR = NULL, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = FALSE, one_vs_manyHR_bool = FALSE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = k0, alphaH = alpha0, betaH = beta0, 
                            kL = k0, alphaL = alpha0_LR, betaL = beta0_LR, 
                            gammaLR_init = z_lr, gammaHR_init = z_hr,
                            seed = seed2, alpha1 = 1.0, alpha2 = 0.0)
  
  ptm <- proc.time() - ptm

  timing["km_nHDP"] <- ptm[3]
  index <- seq(1,save+1, by = skip)
  out_nHDP <- list()
  out_nHDP$highres <- rbind(out_nHDP1$highres[index,], out_nHDP2$highres[index,])
  out_nHDP$lowres <- rbind(out_nHDP1$lowres[index,], out_nHDP2$lowres[index,])
  out_nHDP$highres_mean <- rbind(out_nHDP1$highres_mean[index,], out_nHDP2$highres_mean[index,])
  out_nHDP$lowres_mean <- rbind(out_nHDP1$lowres_mean[index,], out_nHDP2$lowres_mean[index,])
  index <- 1:nrow(out_nHDP$lowres)

  tmp <- get_partition(out_nHDP$highres) # sara's method
  pw_hr <- tmp$pw; z_hr <- tmp$z

  yHRcl_bayes <- colMeans(out_nHDP$highres_mean)
  yHRcl <- as.numeric(yHR)
  for(zk in unique(z_hr)){
    ind <- which(z_hr == zk)
    yHRcl[ind] <-mean(yHRcl[ind])
  }

  tmp <- get_partition(out_nHDP$lowres)
  pw_lr <- tmp$pw; z_lr <- tmp$z
  ## this is because yLR = NULL
  lowres_mean1 <- get_lowres_mean(out_nHDP$lowres, yLR)
  yLRcl_bayes <- colMeans(lowres_mean1)
  ##
  yLRcl <- as.numeric(yLR)
  for(zk in unique(z_lr)){
    ind <- which(z_lr == zk)
    yLRcl[ind] <-mean(yLRcl[ind])
  }

  resultHR["km_nHDP",] <- write_resultsHR_var(yHR, yHRcl, yHRcl_bayes, yHRcl_true, mHR,
                                       z_hr, zHR, out_nHDP$highres, out_nHDP$highres_mean, yHR_outsample)
  resultLR["km_nHDP",] <- write_resultsLR_var(yLR, yLRcl, yLRcl_bayes, yLRcl_true, mLR, Gjs_means,
                                       z_lr, zLR, out_nHDP$lowres, lowres_mean1, yLR_outsample)
  accept["km_nHDP",] <- (out_nHDP1$accept[nchains,] + out_nHDP2$accept[nchains,])/2
  zHR_est["km_nHDP",] <- z_hr
  zLR_est["km_nHDP",] <- z_lr
}

if(method == 2){
  library(cluster)

  ptm <- proc.time()

  d <- dist(Y)
  Kmax <- 20
  nstarts <- 25
  km_silhouette <- numeric(Kmax)
  km2_silhouette <- numeric(Kmax)
  for(k in 2:min(Kmax, nHR-1)){
    km = kmeans(Y, centers = k, nstart = nstarts)
    tmp_sil <- silhouette(x = km$cluster, dist = d) # not defined for k = 1
    km_silhouette[k] <- mean(tmp_sil[, 3])
  }
  
  k_hr <- which.max(km_silhouette)
  km = kmeans(Y, centers = k_hr, nstart = nstarts)
  z_hr <- km$cluster 

  # find proportions of HR clusters in each LR unit, use as data for second km
  prop_hr <- array(0, dim = c(nLR, k_hr))
  for(i in 1:nLR){
    index <- tmp_array[tmp_array[,2] == i,1] # equal to: Mapping_mat[(Mapping_mat[,2] == i - 1),1]+1
    cl_index <- as.numeric(names(table(z_hr[index])))
    prop_hr[i,cl_index] <- as.numeric(table(z_hr[index]))
    prop_hr[i,] <- prop_hr[i,]/sum(prop_hr[i,])
  }
  d2 = dist(prop_hr)

  for(k in 2:min(Kmax, nLR-1,nrow(unique(prop_hr)))){
    km2 = kmeans(prop_hr, centers = k, nstart = nstarts)
    tmp_sil <- silhouette(x = km2$cluster, dist = d2)
    km2_silhouette[k] <- mean(tmp_sil[, 3])
  }

  k_lr <- which.max(km2_silhouette)
  km2 = kmeans(prop_hr, centers = k_lr, nstart = nstarts)
  z_lr <- km2$cluster

  ptm <- proc.time() - ptm

  timing["km"] <- ptm[3]
  out_km <- list()
  out_km$highres <- z_hr
  out_km$lowres <- z_lr
  out_km$silhouette <- km_silhouette
  out_km$silhouette2 <- km2_silhouette
  
  varlist <- c(varlist, "Kmax", "out_km", "nstarts","prop_hr")

  yHRcl <- as.numeric(yHR)
  for(zk in unique(z_hr)){
    ind <- which(z_hr == zk)
    yHRcl[ind] <-mean(yHRcl[ind])
  }

  yLRcl <- as.numeric(yLR)
  for(zk in unique(z_lr)){
    ind <- which(z_lr == zk)
    yLRcl[ind] <-mean(yLRcl[ind])
  }

  resultHR["km",] <- write_resultsHR_var(yHR, yHRcl, NA, yHRcl_true, mHR,
                                       z_hr, zHR, array(NA, dim=c(2,nHR)), array(NA, dim=c(2,nHR)), yHR_outsample)
  resultLR["km",] <- write_resultsLR_var(yLR, yLRcl, NULL, yLRcl_true, mLR, Gjs_means,
                                       z_lr, zLR, array(NA, dim=c(2,nLR)), array(NA, dim=c(2,nLR)), yLR_outsample)
  zHR_est["km",] <- z_hr
  zLR_est["km",] <- z_lr
}

if(method == 3){
  ptm <- proc.time()
  out_nDP1 <- summarize_nDP(dataset, 
    burnin  = burnin, skip = skip, save = save_nDP, 
    eta_LR = eta_LR, eta_HR = eta_CT*eta_TD,
    seed = seed1)
  out_nDP2 <- summarize_nDP(dataset, 
    burnin  = burnin, skip = skip, save = save_nDP, 
    eta_LR = eta_LR, eta_HR = eta_CT*eta_TD,
    seed = seed2)
  ptm <- proc.time() - ptm

  timing["nDP"] <- ptm[3]
  out_nDP <- list()
  out_nDP$highres <- rbind(out_nDP1$highres, out_nDP2$highres)
  out_nDP$lowres <- rbind(out_nDP1$lowres, out_nDP2$lowres)
  out_nDP$highres_mean <- rbind(out_nDP1$highres_mean, out_nDP2$highres_mean)

  tmp <- get_partition(out_nDP$highres) # sara's method
  pw_hr <- tmp$pw; z_hr <- tmp$z

  yHRcl_bayes <- colMeans(out_nDP$highres_mean)
  yHRcl <- as.numeric(yHR)
  for(zk in unique(z_hr)){
    ind <- which(z_hr == zk)
    yHRcl[ind] <-mean(yHRcl[ind])
  }

  tmp <- get_partition(out_nDP$lowres)
  pw_lr <- tmp$pw; z_lr <- tmp$z
  ## this is because yLR = NULL
  lowres_mean1 <- get_lowres_mean(out_nDP$lowres, yLR)
  yLRcl_bayes <- colMeans(lowres_mean1)
  ##
  yLRcl <- as.numeric(yLR)
  for(zk in unique(z_lr)){
    ind <- which(z_lr == zk)
    yLRcl[ind] <-mean(yLRcl[ind])
  }

  resultHR["nDP",] <- write_resultsHR_var(yHR, yHRcl, yHRcl_bayes, yHRcl_true, mHR,
                                       z_hr, zHR, out_nDP$highres, out_nDP$highres_mean, yHR_outsample)
  resultLR["nDP",] <- write_resultsLR_var(yLR, yLRcl, yLRcl_bayes, yLRcl_true, mLR, Gjs_means,
                                       z_lr, zLR, out_nDP$lowres, lowres_mean1, yLR_outsample)
  zHR_est["nDP",] <- z_hr
  zLR_est["nDP",] <- z_lr
}

### these are examples, not really considered in the final simulations reported in the paper

# nHDP1_0
if(method == 4){
  ## fixed eta
  ptm <- proc.time()
  out_nHDP1 <- paralleltemp(Y, YLR = NULL, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = FALSE, one_vs_manyHR_bool = FALSE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = k0, alphaH = alpha0, betaH = beta0, 
                            kL = k0, alphaL = alpha0_LR, betaL = beta0_LR, 
                            seed = seed1, alpha1 = 1.0, alpha2 = 0.0)
  out_nHDP2 <- paralleltemp(Y, YLR = NULL, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = FALSE, one_vs_manyHR_bool = FALSE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = k0, alphaH = alpha0, betaH = beta0, 
                            kL = k0, alphaL = alpha0_LR, betaL = beta0_LR, 
                            seed = seed2, alpha1 = 1.0, alpha2 = 0.0)
  ptm <- proc.time() - ptm
  timing["nHDP1_0"] <- ptm[3]
  index <- seq(1,save+1, by = skip)
  out_nHDP <- list()
  out_nHDP$highres <- rbind(out_nHDP1$highres[index,], out_nHDP2$highres[index,])
  out_nHDP$lowres <- rbind(out_nHDP1$lowres[index,], out_nHDP2$lowres[index,])
  out_nHDP$highres_mean <- rbind(out_nHDP1$highres_mean[index,], out_nHDP2$highres_mean[index,])
  out_nHDP$lowres_mean <- rbind(out_nHDP1$lowres_mean[index,], out_nHDP2$lowres_mean[index,])
  index <- 1:nrow(out_nHDP$lowres)

  tmp <- get_partition(out_nHDP$highres) # sara's method
  pw_hr <- tmp$pw; z_hr <- tmp$z

  yHRcl_bayes <- colMeans(out_nHDP$highres_mean)
  yHRcl <- as.numeric(yHR)
  for(zk in unique(z_hr)){
    ind <- which(z_hr == zk)
    yHRcl[ind] <-mean(yHRcl[ind])
  }

  tmp <- get_partition(out_nHDP$lowres)
  pw_lr <- tmp$pw; z_lr <- tmp$z
  ## this is because yLR = NULL
  lowres_mean1 <- get_lowres_mean(out_nHDP$lowres, yLR)
  yLRcl_bayes <- colMeans(lowres_mean1)
  ##
  yLRcl <- as.numeric(yLR)
  for(zk in unique(z_lr)){
    ind <- which(z_lr == zk)
    yLRcl[ind] <-mean(yLRcl[ind])
  }

  resultHR["nHDP1_0",] <- write_resultsHR_var(yHR, yHRcl, yHRcl_bayes, yHRcl_true, mHR,
                                       z_hr, zHR, out_nHDP$highres, out_nHDP$highres_mean, yHR_outsample)
  resultLR["nHDP1_0",] <- write_resultsLR_var(yLR, yLRcl, yLRcl_bayes, yLRcl_true, mLR, Gjs_means,
                                       z_lr, zLR, out_nHDP$lowres, lowres_mean1, yLR_outsample)
  accept["nHDP1_0",] <- (out_nHDP1$accept[nchains,] + out_nHDP2$accept[nchains,])/2
  zHR_est["nHDP1_0",] <- z_hr
  zLR_est["nHDP1_0",] <- z_lr
}

# nHDP1_0_sampleeta
if(method == 5){
  ## sample eta
  ptm <- proc.time()
  out_nHDP1 <- paralleltemp(Y, YLR = NULL, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = FALSE, one_vs_manyHR_bool = FALSE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            grid_eta = grid, prior_eta = prior, 
                            kH = k0, alphaH = alpha0, betaH = beta0, 
                            kL = k0, alphaL = alpha0_LR, betaL = beta0_LR, 
                            seed = seed1, alpha1 = 1.0, alpha2 = 0.0)
  out_nHDP2 <- paralleltemp(Y, YLR = NULL, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = FALSE, one_vs_manyHR_bool = FALSE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            grid_eta = grid, prior_eta = prior, 
                            kH = k0, alphaH = alpha0, betaH = beta0, 
                            kL = k0, alphaL = alpha0_LR, betaL = beta0_LR, 
                            seed = seed2, alpha1 = 1.0, alpha2 = 0.0)
  ptm <- proc.time() - ptm
  timing["nHDP1_0_sampleeta"] <- ptm[3]
  index <- seq(1,save+1, by = skip)
  out_nHDP <- list()
  out_nHDP$highres <- rbind(out_nHDP1$highres[index,], out_nHDP2$highres[index,])
  out_nHDP$lowres <- rbind(out_nHDP1$lowres[index,], out_nHDP2$lowres[index,])
  out_nHDP$highres_mean <- rbind(out_nHDP1$highres_mean[index,], out_nHDP2$highres_mean[index,])
  out_nHDP$lowres_mean <- rbind(out_nHDP1$lowres_mean[index,], out_nHDP2$lowres_mean[index,])
  index <- 1:nrow(out_nHDP$lowres)

  tmp <- get_partition(out_nHDP$highres) # sara's method
  pw_hr <- tmp$pw; z_hr <- tmp$z

  yHRcl_bayes <- colMeans(out_nHDP$highres_mean)
  yHRcl <- as.numeric(yHR)
  for(zk in unique(z_hr)){
    ind <- which(z_hr == zk)
    yHRcl[ind] <-mean(yHRcl[ind])
  }

  tmp <- get_partition(out_nHDP$lowres)
  pw_lr <- tmp$pw; z_lr <- tmp$z
  ## this is because yLR = NULL
  lowres_mean1 <- get_lowres_mean(out_nHDP$lowres, yLR)
  yLRcl_bayes <- colMeans(lowres_mean1)
  ##
  yLRcl <- as.numeric(yLR)
  for(zk in unique(z_lr)){
    ind <- which(z_lr == zk)
    yLRcl[ind] <-mean(yLRcl[ind])
  }

  resultHR["nHDP1_0_sampleeta",] <- write_resultsHR_var(yHR, yHRcl, yHRcl_bayes, yHRcl_true, mHR,
                                       z_hr, zHR, out_nHDP$highres, out_nHDP$highres_mean, yHR_outsample)
  resultLR["nHDP1_0_sampleeta",] <- write_resultsLR_var(yLR, yLRcl, yLRcl_bayes, yLRcl_true, mLR, Gjs_means,
                                       z_lr, zLR, out_nHDP$lowres, lowres_mean1, yLR_outsample)
  accept["nHDP1_0_sampleeta",] <- (out_nHDP1$accept[nchains,] + out_nHDP2$accept[nchains,])/2
  zHR_est["nHDP1_0_sampleeta",] <- z_hr
  zLR_est["nHDP1_0_sampleeta",] <- z_lr
}

# HDPoracle
if(method == 6){
  ptm <- proc.time()
  out_nHDP1 <- paralleltemp_noLR(Y, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = FALSE, one_vs_manyHR_bool = FALSE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = k0, alphaH = alpha0, betaH = beta0, 
                            seed = seed1, gammaLR_init = zLR)
  out_nHDP2 <- paralleltemp_noLR(Y, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = FALSE, one_vs_manyHR_bool = FALSE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = k0, alphaH = alpha0, betaH = beta0, 
                            seed = seed2, gammaLR_init = zLR)
  ptm <- proc.time() - ptm
  timing["HDPoracle"] <- ptm[3]
  index <- seq(1,save+1, by = skip)
  out_nHDP <- list()
  out_nHDP$highres <- rbind(out_nHDP1$highres[index,], out_nHDP2$highres[index,])
  out_nHDP$lowres <- rbind(out_nHDP1$lowres[index,], out_nHDP2$lowres[index,])
  out_nHDP$highres_mean <- rbind(out_nHDP1$highres_mean[index,], out_nHDP2$highres_mean[index,])
  out_nHDP$lowres_mean <- rbind(out_nHDP1$lowres_mean[index,], out_nHDP2$lowres_mean[index,])
  index <- 1:nrow(out_nHDP$lowres)

  tmp <- get_partition(out_nHDP$highres) # sara's method
  pw_hr <- tmp$pw; z_hr <- tmp$z

  yHRcl_bayes <- colMeans(out_nHDP$highres_mean)
  yHRcl <- as.numeric(yHR)
  for(zk in unique(z_hr)){
    ind <- which(z_hr == zk)
    yHRcl[ind] <-mean(yHRcl[ind])
  }

  z_lr <- zLR

  ## this is because yLR = NULL
  lowres_mean1 <- get_lowres_mean(out_nHDP$lowres, yLR)
  yLRcl_bayes <- colMeans(lowres_mean1)
  ##
  yLRcl <- as.numeric(yLR)
  for(zk in unique(z_lr)){
    ind <- which(z_lr == zk)
    yLRcl[ind] <-mean(yLRcl[ind])
  }

  resultHR["HDPoracle",] <- write_resultsHR_var(yHR, yHRcl, yHRcl_bayes, yHRcl_true, mHR,
                                       z_hr, zHR, out_nHDP$highres, out_nHDP$highres_mean, yHR_outsample)
  resultLR["HDPoracle",] <- write_resultsLR_var(yLR, yLRcl, yLRcl_bayes, yLRcl_true, mLR, Gjs_means,
                                       z_lr, zLR, out_nHDP$lowres, lowres_mean1, yLR_outsample)
  accept["HDPoracle",] <- (out_nHDP1$accept[nchains,] + out_nHDP2$accept[nchains,])/2
  zHR_est["HDPoracle",] <- z_hr
  zLR_est["HDPoracle",] <- z_lr
}

# HDPone
if(method == 7){
  ptm <- proc.time()
  out_nHDP1 <- paralleltemp_noLR(Y, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = TRUE, one_vs_manyHR_bool = FALSE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = k0, alphaH = alpha0, betaH = beta0, 
                            seed = seed1)
  out_nHDP2 <- paralleltemp_noLR(Y, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = TRUE, one_vs_manyHR_bool = FALSE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = k0, alphaH = alpha0, betaH = beta0, 
                            seed = seed2)
  ptm <- proc.time() - ptm
  timing["HDPone"] <- ptm[3]
  index <- seq(1,save+1, by = skip)
  out_nHDP <- list()
  out_nHDP$highres <- rbind(out_nHDP1$highres[index,], out_nHDP2$highres[index,])
  out_nHDP$lowres <- rbind(out_nHDP1$lowres[index,], out_nHDP2$lowres[index,])
  out_nHDP$highres_mean <- rbind(out_nHDP1$highres_mean[index,], out_nHDP2$highres_mean[index,])
  out_nHDP$lowres_mean <- rbind(out_nHDP1$lowres_mean[index,], out_nHDP2$lowres_mean[index,])
  index <- 1:nrow(out_nHDP$lowres)

  tmp <- get_partition(out_nHDP$highres) # sara's method
  pw_hr <- tmp$pw; z_hr <- tmp$z

  yHRcl_bayes <- colMeans(out_nHDP$highres_mean)
  yHRcl <- as.numeric(yHR)
  for(zk in unique(z_hr)){
    ind <- which(z_hr == zk)
    yHRcl[ind] <-mean(yHRcl[ind])
  }

  z_lr <- rep(1, n_ct)

  ## this is because yLR = NULL
  lowres_mean1 <- get_lowres_mean(out_nHDP$lowres, yLR)
  yLRcl_bayes <- colMeans(lowres_mean1)
  ##
  yLRcl <- as.numeric(yLR)
  for(zk in unique(z_lr)){
    ind <- which(z_lr == zk)
    yLRcl[ind] <-mean(yLRcl[ind])
  }

  resultHR["HDPone",] <- write_resultsHR_var(yHR, yHRcl, yHRcl_bayes, yHRcl_true, mHR,
                                       z_hr, zHR, out_nHDP$highres, out_nHDP$highres_mean, yHR_outsample)
  resultLR["HDPone",] <- write_resultsLR_var(yLR, yLRcl, yLRcl_bayes, yLRcl_true, mLR, Gjs_means,
                                       z_lr, zLR, out_nHDP$lowres, lowres_mean1, yLR_outsample)
  accept["HDPone",] <- (out_nHDP1$accept[nchains,] + out_nHDP2$accept[nchains,])/2
  zHR_est["HDPone",] <- z_hr
  zLR_est["HDPone",] <- z_lr
}

# HDPmany
if(method == 8){
  ptm <- proc.time()
  out_nHDP1 <- paralleltemp_noLR(Y, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = FALSE, one_vs_manyHR_bool = FALSE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = k0, alphaH = alpha0, betaH = beta0, 
                            seed = seed1)
  out_nHDP2 <- paralleltemp_noLR(Y, Mapping, InvMapping, n, n2, burnin = burnin, adjust = adjust, 
                            save = save, gibbs = gibbs, 
                            nchains = nchains, temps = temps,
                            print_bool = FALSE, 
                            one_vs_manyLR_bool = FALSE, one_vs_manyHR_bool = FALSE,
                            eta_CT = eta_CT, eta_TD = eta_TD, eta_LR = eta_LR,
                            kH = k0, alphaH = alpha0, betaH = beta0, 
                            seed = seed2)
  ptm <- proc.time() - ptm
  timing["HDPmany"] <- ptm[3]
  index <- seq(1,save+1, by = skip)
  out_nHDP <- list()
  out_nHDP$highres <- rbind(out_nHDP1$highres[index,], out_nHDP2$highres[index,])
  out_nHDP$lowres <- rbind(out_nHDP1$lowres[index,], out_nHDP2$lowres[index,])
  out_nHDP$highres_mean <- rbind(out_nHDP1$highres_mean[index,], out_nHDP2$highres_mean[index,])
  out_nHDP$lowres_mean <- rbind(out_nHDP1$lowres_mean[index,], out_nHDP2$lowres_mean[index,])
  index <- 1:nrow(out_nHDP$lowres)

  tmp <- get_partition(out_nHDP$highres) # sara's method
  pw_hr <- tmp$pw; z_hr <- tmp$z

  yHRcl_bayes <- colMeans(out_nHDP$highres_mean)
  yHRcl <- as.numeric(yHR)
  for(zk in unique(z_hr)){
    ind <- which(z_hr == zk)
    yHRcl[ind] <-mean(yHRcl[ind])
  }

  z_lr <- 1:n_ct

  ## this is because yLR = NULL
  lowres_mean1 <- get_lowres_mean(out_nHDP$lowres, yLR)
  yLRcl_bayes <- colMeans(lowres_mean1)
  ##
  yLRcl <- as.numeric(yLR)
  for(zk in unique(z_lr)){
    ind <- which(z_lr == zk)
    yLRcl[ind] <-mean(yLRcl[ind])
  }

  resultHR["HDPmany",] <- write_resultsHR_var(yHR, yHRcl, yHRcl_bayes, yHRcl_true, mHR,
                                       z_hr, zHR, out_nHDP$highres, out_nHDP$highres_mean, yHR_outsample)
  resultLR["HDPmany",] <- write_resultsLR_var(yLR, yLRcl, yLRcl_bayes, yLRcl_true, mLR, Gjs_means,
                                       z_lr, zLR, out_nHDP$lowres, lowres_mean1, yLR_outsample)
  accept["HDPmany",] <- (out_nHDP1$accept[nchains,] + out_nHDP2$accept[nchains,])/2
  zHR_est["HDPmany",] <- z_hr
  zLR_est["HDPmany",] <- z_lr
}


filename <- paste0("results/sim/modeldata/nCT",n_ct,"nBGinCT",n_bg_in_ct,"_etas",alpha_LR,gamma_HR,alpha0_HR,"_save",save,"_dist",dists_str[opt_dist],"_met",method,"sim",sim,".rdata")
varlist <- c(varlist, "timing","resultLR","resultHR","startingpoint", "accept","zLR_est","zHR_est",
             "yHR", "mHR", "yLR", "mLR", "zHR","rests", "zLR", "zLR_HR", 
             "zLR_emp", "zLR_HR_emp", "propHR", "prop")
save(list = varlist, file = filename)
