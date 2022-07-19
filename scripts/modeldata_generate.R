require(mcclust.ext)
require(mcclust)
require(spdep)

source("R/functions_plot.R")
source("R/functions_analysis.R")

## sim_number should be ranged from 1 to 50.
## config 1 corresponds to etas = 3, 5, 3; config 2 corresponds to etas = 1, 1, 1
args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
sim <- as.character(sim_number)
config <- as.numeric(args[2])

cat(sim_number, "\n")
cat(sim, "\n")
cat(config, "\n")

n_ct = 10                  # number of low-res units
n_bg_in_ct = 100           # we use the largest number we would consider
n_bg <- n_ct * n_bg_in_ct  # total number of high-res units
bg_in_ct <- t(array(1:n_bg, dim = c(n_bg_in_ct, n_ct)))

tmp_array <- array(NA, dim = c(n_bg,2))
for(i in 1:n_ct){
  for(x in bg_in_ct[i,]){
    tmp_array[x, ] <- c(x, i)
  }
}

nHR <- n_bg
nLR <- n_ct

TASK <- sim_number
set.seed(609 + TASK)
seed1 <- sample(1000000, 1)
seed2 <- sample(1000000, 1)

########
if(config == 1){
  alpha_LR = 3
  gamma_HR = 5   # G_0 ~ DP(gamma, H)    # eta_TD
  alpha0_HR = 3  # G_i ~ DP(alpha0, G0)  # eta_CT
} else if(config == 2){
  alpha_LR = 1
  gamma_HR = 1
  alpha0_HR = 1
}

pLR <- stick_breaking(alpha_LR, nLR)
x <- rmultinom(nLR, size = 1, prob = pLR)
zLR <- apply(x, 2, function(x) which(x == 1))
zLR <- reorder_part(zLR)

kLR <- length(unique(zLR))
zLR_HR <- get_ids_HR_col(zLR, tmp_array[,2])

## New way to generate Gjs so that they are as diverse as possible
eps <- 0.8
Gjs <- stick_breaking_HDP(gamma_HR, alpha0_HR, nHR, 5*kLR)
if(kLR == 1){
  max_col <- max(which(Gjs != 0))
} else {
  max_col <- max(which(colSums(Gjs) != 0))
}
Gjs_restr <- Gjs[,1:max_col]
new_rest <- find_unique_dist(Gjs_restr, eps)
new_n <- length(unique(new_rest))
Gjs_restr_new <- array(0, dim = c(new_n, nHR))
for(j in 1:new_n){
  Gjs_restr_new[j,1:max_col] <- Gjs_restr[which(new_rest == j)[1],]
}
while(new_n < kLR){
  Gjs2 <- stick_breaking_HDP(gamma_HR, alpha0_HR, nHR, 5*kLR)
  Gjs <- rbind(Gjs_restr_new, Gjs2)
  if(kLR == 1){
    max_col <- max(which(Gjs != 0))
  } else {
    max_col <- max(which(colSums(Gjs) != 0))
  }
  Gjs_restr <- Gjs[,1:max_col]
  new_rest <- find_unique_dist(Gjs_restr, eps)
  new_n <- length(unique(new_rest))
  Gjs_restr_new <- array(0, dim = c(new_n, nHR))
  for(j in 1:new_n){
    Gjs_restr_new[j,1:max_col] <- Gjs_restr[which(new_rest == j)[1],]
  }
}
Gjs <- Gjs_restr_new[1:kLR,, drop = FALSE]

filename <- paste0("data/modeldata/nCT",n_ct,"nBGinCT",n_bg_in_ct,"_etas",alpha_LR,gamma_HR,alpha0_HR,"_sim",sim,".rdata")
varlist <- ls()
save(list = varlist, file = filename)
