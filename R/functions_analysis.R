## Compute prior for DP/PYP, and HDP

log_pyp_prior <- function(part, alpha, sigma){
  # for each cluster we need (1-sigma)_{n_i-1} = Gamma(n_i - sigma)/Gamma(1-sigma)
  n <- length(part)
  x <- as.numeric(table(part))
  l <- sum(lgamma(x-sigma)-lgamma(1-sigma))
  l <- l + sum(log(alpha + (0:(length(x)-1))*sigma))
  l <- l - sum(log(alpha + 0:(n-1)))
  return(l)
}

log_pyp_prior2 <- function(ncl, alpha, sigma){
  # for each cluster we need (1-sigma)_{n_i-1} = Gamma(n_i - sigma)/Gamma(1-sigma)
  n <- sum(ncl)
  x <- ncl
  l <- sum(lgamma(x-sigma)-lgamma(1-sigma))
  l <- l + sum(log(alpha + (0:(length(x)-1))*sigma))
  l <- l - sum(log(alpha + 0:(n-1)))
  return(l)
}


log_hdp_prior <- function(parts_mat, c, c0){
  d <- nrow(parts_mat)
  k <- ncol(parts_mat)
  Ns <- rowSums(parts_mat) # number of elements in each restaurant
  # abs(Stirling1(n,k)) # log(abs(Stirling1.all(200)))
  # Ls_mat has same dimension as parts_mat, but contains the Lij
  lprob <- k * log(c0) - sum(lgamma(c + Ns) - lgamma(c))
  # we should have the sum over Ls_mat
  starting_Ls_mat <- parts_mat*0 + 1
  if(any(parts_mat == 0)){
    tmp <-which(parts_mat==0, arr.ind = T)
    starting_Ls_mat[tmp] <- 0
  }
  lprob <- lprob + big_sum(starting_Ls_mat, parts_mat, c, c0)
  
  return(lprob)
}

big_sum <- function(Ls_mat, parts_mat, c, c0){
  if(any(Ls_mat > parts_mat)){
    return(0)
  } else {
    # cat(c(Ls_mat), "\n")
    lprob <- sum(Ls_mat) * log(c) - (lgamma(c0 + sum(Ls_mat)) - lgamma(c0))
    lprob <- lprob + sum(lgamma(colSums(Ls_mat)))
    lprob <- lprob + sum(apply(cbind(as.numeric(parts_mat), as.numeric(Ls_mat)), MARGIN = 1, function(x) log(abs(Stirling1(x[1],x[2])))))
    s <- sum(sapply(1:length(Ls_mat), function(i){add_this <- Ls_mat*0; add_this[i] <- 1; big_sum(Ls_mat + add_this, parts_mat, c, c0)}))
    return(lprob +  s)
  }
}

## Create a vector containing the results to be saved for the simulation analysis
write_resultsHR_var <- function(yHR, yHRcl, yHRcl_bayes, yHRcl_true, mHR,
                          z_hr, zHR, out_nHDP_highres, out_nHDP_highres_mean, yHR_outsample){
  names_str <- c("R2_est","R2_bayes","R2_orac",
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
  vec <- numeric(length(names_str))
  names(vec) <- names_str
  vec["R2_est"]     <- 1 - (mean((yHR - yHRcl)^2)/var(yHR))
  vec["R2_bayes"]   <- 1 - (mean((yHR - yHRcl_bayes)^2)/var(yHR))
  vec["R2_orac"]    <- 1 - (mean((yHR - yHRcl_true)^2)/var(yHR))
  vec["R2_bayes_mean"] <- mean(apply(out_nHDP_highres_mean, MARGIN = 1, FUN = function(x) {1-(mean((yHR - x)^2)/var(yHR))}))
  vec["R2_bayes_sd"]   <- sd(apply(out_nHDP_highres_mean, MARGIN = 1, FUN = function(x) {1-(mean((yHR - x)^2)/var(yHR))}))
  vec["R2_OS_est"]   <- 1 - (mean((yHR_outsample - yHRcl)^2)/var(yHR_outsample))
  vec["R2_OS_bayes"]   <- 1 - (mean((yHR_outsample - yHRcl_bayes)^2)/var(yHR_outsample))
  vec["R2_OS_orac"]   <- 1 - (mean((yHR_outsample - yHRcl_true)^2)/var(yHR_outsample))

  vec["MSE_est"]    <- mean((mHR - yHRcl)^2)
  vec["MSE_bayes"]  <- mean((mHR - yHRcl_bayes)^2)
  vec["MSE_orac"]   <- mean((mHR - yHRcl_true)^2)
  vec["MSE_bayes_mean"] <- mean(apply(out_nHDP_highres_mean, MARGIN = 1, FUN = function(x) {mean((mHR - x)^2)}))
  vec["MSE_bayes_sd"]   <- sd(apply(out_nHDP_highres_mean, MARGIN = 1, FUN = function(x) {mean((mHR - x)^2)}))
  
  vec["MSE_Gjs_est"]    <- NA
  vec["MSE_Gjs_bayes"]  <- NA
  vec["MSE_Gjs_orac"]   <- NA

  vec["VarY"]       <- var(yHR)
  vec["B"]          <- binder(z_hr,zHR)
  binders <- apply(out_nHDP_highres, MARGIN = 1, FUN = function(x) binder(x,zHR))
  vec["B_bayes"]    <- mean(binders)
  vec["B_bayes_sd"] <- sd(binders)
  binders <- apply(out_nHDP_highres, MARGIN = 1, FUN = function(x) binder(x,z_hr))
  vec["B_bayes2"]    <- mean(binders)
  vec["B_bayes_sd2"] <- sd(binders)

  vec["R"] <- arandi(z_hr, zHR, adjust = FALSE)
  vec["Radj"] <- arandi(z_hr, zHR, adjust = TRUE)
  Rs <- apply(out_nHDP_highres, MARGIN = 1, FUN = function(x) arandi(x,zHR, adjust = FALSE))
  Radjs <- apply(out_nHDP_highres, MARGIN = 1, FUN = function(x) arandi(x,zHR, adjust = TRUE))
  vec["R_bayes"]    <- mean(Rs)
  vec["R_bayes_sd"]    <- sd(Rs)
  vec["Radj_bayes"] <- mean(Radjs)
  vec["Radj_bayes_sd"]    <- sd(Radjs)

  vec["VI_est"] <- vi.dist(z_hr,zHR)
  vis <- apply(out_nHDP_highres, MARGIN = 1, FUN = function(x) vi.dist(x,zHR))
  vec["VI_bayes"] <- mean(vis)
  vec["VI_bayes_sd"] <- sd(vis)
  
  vec["K_est"] <- length(unique(z_hr))
  Ks <- apply(out_nHDP_highres, MARGIN = 1, FUN = function(x) length(unique(x)))
  vec["K_bayes"] <- mean(Ks)
  vec["K_bayes_sd"] <- sd(Ks)
  vec["K_true"] <- length(unique(zHR)) # before it was (typo) z_hr
  
  vec["lprior_est"] <- log_pyp_prior(z_hr, alpha = eta_TD, sigma = 0)
  vec["lprior_bayes"]   <- median(apply(out_nHDP_highres, MARGIN = 1, FUN = function(x) log_pyp_prior(x, alpha = eta_TD, sigma = 0)))
  vec["lprior_bayes_iqr"]   <- IQR(apply(out_nHDP_highres, MARGIN = 1, FUN = function(x) log_pyp_prior(x, alpha = eta_TD, sigma = 0)))

  return(vec)
}
write_resultsLR_var <- function(yLR, yLRcl, yLRcl_bayes, yLRcl_true, mLR, Gjs_means,
                          z_lr, zLR, 
                          out_nHDP_lowres, out_nHDP_lowres_mean, yLR_outsample){
  names_str <- c("R2_est","R2_bayes","R2_orac",
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
  vec <- numeric(length(names_str))
  names(vec) <- names_str
  vec["R2_est"]     <- 1 - (mean((yLR - yLRcl)^2)/var(yLR))
  vec["R2_bayes"]   <- 1 - (mean((yLR - yLRcl_bayes)^2)/var(yLR))
  vec["R2_orac"]    <- 1 - (mean((yLR - yLRcl_true)^2)/var(yLR))
  vec["R2_bayes_mean"] <- mean(apply(out_nHDP_lowres_mean, MARGIN = 1, FUN = function(x) {1-(mean((yLR - x)^2)/var(yLR))}))
  vec["R2_bayes_sd"]   <- sd(apply(out_nHDP_lowres_mean, MARGIN = 1, FUN = function(x) {1-(mean((yLR - x)^2)/var(yLR))}))
  vec["R2_OS_est"]   <- 1 - (mean((yLR_outsample - yLRcl)^2)/var(yLR_outsample))
  vec["R2_OS_bayes"]   <- 1 - (mean((yLR_outsample - yLRcl_bayes)^2)/var(yLR_outsample))
  vec["R2_OS_orac"]   <- 1 - (mean((yLR_outsample - yLRcl_true)^2)/var(yLR_outsample))

  vec["MSE_est"]    <- mean((mLR - yLRcl)^2)
  vec["MSE_bayes"]  <- mean((mLR - yLRcl_bayes)^2)
  vec["MSE_orac"]   <- mean((mLR - yLRcl_true)^2)
  vec["MSE_bayes_mean"] <- mean(apply(out_nHDP_lowres_mean, MARGIN = 1, FUN = function(x) {mean((mLR - x)^2)}))
  vec["MSE_bayes_sd"]   <- sd(apply(out_nHDP_lowres_mean, MARGIN = 1, FUN = function(x) {mean((mLR - x)^2)}))

  vec["MSE_Gjs_est"]    <- mean((Gjs_means - yLRcl)^2)
  vec["MSE_Gjs_bayes"]  <- mean((Gjs_means - yLRcl_bayes)^2)
  vec["MSE_Gjs_orac"]   <- mean((Gjs_means - yLRcl_true)^2)
  
  vec["VarY"]       <- var(yLR)
  vec["B"]       <- binder(z_lr,zLR)
  binders <- apply(out_nHDP_lowres, MARGIN = 1, FUN = function(x) binder(x,zLR)) ## distance between truth
  vec["B_bayes"]       <- mean(binders)
  vec["B_bayes_sd"]    <- sd(binders)
  binders <- apply(out_nHDP_lowres, MARGIN = 1, FUN = function(x) binder(x,z_lr)) ## distance between estimate
  vec["B_bayes2"]       <- mean(binders)
  vec["B_bayes_sd2"]    <- sd(binders)

  vec["R"] <- arandi(z_lr, zLR, adjust = FALSE)
  vec["Radj"] <- arandi(z_lr, zLR, adjust = TRUE)
  Rs <- apply(out_nHDP_lowres, MARGIN = 1, FUN = function(x) arandi(x,zLR, adjust = FALSE))
  Radjs <- apply(out_nHDP_lowres, MARGIN = 1, FUN = function(x) arandi(x,zLR, adjust = TRUE))
  vec["R_bayes"]    <- mean(Rs)
  vec["R_bayes_sd"]    <- sd(Rs)
  vec["Radj_bayes"] <- mean(Radjs)
  vec["Radj_bayes_sd"]    <- sd(Radjs)

  vec["VI_est"] <- vi.dist(z_lr,zLR)
  vis <- apply(out_nHDP_lowres, MARGIN = 1, FUN = function(x) vi.dist(x,zLR))
  vec["VI_bayes"] <- mean(vis)
  vec["VI_bayes_sd"] <- sd(vis)
  
  vec["K_est"] <- length(unique(z_lr))
  Ks <- apply(out_nHDP_lowres, MARGIN = 1, FUN = function(x) length(unique(x)))
  vec["K_bayes"] <- mean(Ks)
  vec["K_bayes_sd"] <- sd(Ks)
  vec["K_true"] <- length(unique(zLR)) # before it was (typo) z_lr

  vec["lprior_est"] <- log_pyp_prior(z_lr, alpha = eta_LR, sigma = 0)
  vec["lprior_bayes"]   <- median(apply(out_nHDP_lowres, MARGIN = 1, FUN = function(x) log_pyp_prior(x, alpha = eta_LR, sigma = 0)))
  vec["lprior_bayes_iqr"]   <- IQR(apply(out_nHDP_lowres, MARGIN = 1, FUN = function(x) log_pyp_prior(x, alpha = eta_LR, sigma = 0)))
  return(vec)
}

get_lowres_mean <- function(lowres, yLR){
  n <- ncol(lowres)
  Niter <- nrow(lowres)
  lowres_mean <- array(NA, dim = c(Niter, n))
  for(i in 1:Niter){
    x <- lowres[i,]
    x_unik <- unique(x)
    for(xi in x_unik){
      ind <- which(x == xi)
      lowres_mean[i, ind] <- mean(yLR[ind])
    }
  }
  return(lowres_mean)
}


## Various functions to  compute properties of partitions

nclusters <- function(zs){
  # for each partition returns the number of unique clusters
  # works both if there is only one (input is a vector) or if there are many (input is an array)
  if(is.null(dim(zs))) return(length(unique(zs)))
  else apply(zs, 1, function(x) length(unique(x)))
}
sizes_clusters <- function(zs){
  apply(zs, MARGIN = 1, FUN = function(z){
    as.numeric(table(z))
  })
}
size_clusters <- function(zs, x){
  apply(zs, MARGIN = 1, FUN = function(z){
    sum(z == z[x])
  })
}
pairwise <- function(z){
  # give one vector of cluster membership returns the pairwise matrix
  n <- length(z)
  mat <- matrix(data = 0, nrow = n, ncol = n)
  zk <- unique(z)
  for(x in zk){
    ind <- which(z == x)
    mat[ind, ind] <- 1
  }
  return(mat)
}
pairwise_vec <- function(z){
  mat <- pairwise(z)
  vec <- mat[lower.tri(mat)]
  return(vec)
}
pairwises <- function(zs, index = NULL){
  ## NOTE: mcclust::comp.psm is WAY FASTER
  cat("Use mcclust::comp.psm instead!")
  # given an array (or vector) of cluster memberships returns the average 
  n <- ncol(zs)
  if(!is.null(index)){
    zs <- zs[index,]
  }
  vec <- rowMeans(apply(zs, MARGIN = 1, FUN = pairwise))
  matrix(vec, nrow = n, ncol = n, byrow = T) # tanto e' simmetrica
}
pairwises_TS <- function(zs, index = NULL){
  n <- ncol(zs)
  if(!is.null(index)){
    zs <- zs[index,]
  }
  # each column is a vector
  TS <- apply(X = zs, MARGIN = 1, FUN = pairwise_vec)
  return(TS)
}

largest_cl <- function(z, i = 1){
  t <- table(z)
  # cl_names <- as.numeric(names(t))
  cl_freq <- as.numeric(t)
  tmp <- sort(cl_freq, decreasing = T, index.return = T)
  freq <- tmp$x[i]
  if(is.na(freq)) freq <- 0
  # lab <- cl_names[tmp$ix[i]]
  return(freq)
}

largest_cl_Y <- function(z, i = 1, Y = NULL){
  t <- table(z)
  cl_names <- as.numeric(names(t))
  cl_freq <- as.numeric(t)
  tmp <- sort(cl_freq, decreasing = T, index.return = T)
  lab <- cl_names[tmp$ix[i]]
  index <- which(z == lab)
  if(!is.null(Y)){
    m <- mean(Y[index])
    s <- sd(Y[index])
  }
  # return(list(mean = m, sd = s))
  return(m)
}

cl_Y <- function(z, cl, Y = NULL){
  index <- which(z == cl)
  if(!is.null(Y)){
    m <- mean(Y[index])
  }
  return(m)
}

Y_cl_mean_sorted <- function(z, Y, ncl = NULL){
  if(is.null(ncl)) ncl = length(z)
  z_un <- unique(z)
  K <- length(z_un)
  means <- rep(NA, ncl)
  for(k in 1:min(K, ncl)){
    zk <- z_un[k]
    ind <- which(z == zk)
    means[k] <- mean(Y[ind])
  }
  vec <- sort(means, decreasing = T, na.last = T)
  return(vec)
}

get_partition <- function(samples, hier = FALSE, K = NULL, Kmode = FALSE){
  # samples should be out$lowres[index,], or out$highres[index,]
  if(hier == FALSE){
    if(!is.null(K)){
      ind <- which(nclusters(samples) == K)
      pw_hr <- pairwises(samples[ind,])
      z_hr <- minVI(pw_hr, cls.draw = samples[ind,], method = "draws")$cl
    } else {
      if(Kmode){
        k <- as.numeric(names(which.max(table(nclusters(samples)))))
        ind <- which(nclusters(samples) == k)
        pw_hr <- pairwises(samples[ind,])
        z_hr <- minVI(pw_hr, cls.draw = samples[ind,], method = "draws")$cl
      } else {
        pw_hr <- mcclust::comp.psm(samples)
        # pw_hr <- pairwises(samples)
        z_hr <- minVI(pw_hr, cls.draw = samples, method = "draws")$cl
      }
    }
  } else {
    if(is.null(K)){
      K <- as.numeric(names(which.max(table(nclusters(samples)))))
    }
    pw_hr <- pairwises(samples)
    d <- dist(pw_hr)
    hc <- hclust(d)
    z_hr <- cutree(hc, k = K)
  }
  return(list(z = z_hr, pw = pw_hr))
}

binder <- function(c1,c2){
  n <- length(c1)
  M1 <- pairwise(c1)
  M2 <- pairwise(c2)
  sum(M1 * (1-M2) + (1-M1) * M2)/(n*(n-1))
}

## Functions to generate partitions

PYP <- function(n, alpha = 1, sigma = 0){
  # conditional sampling to get DP samples
  ids <- numeric(n)
  ids[1] <- 1
  k <- 1
  ids_un <- c(1)
  n_k <- c(1)
  if(n > 1){
    for(i in 2:n){
      # ids_un <- unique(ids[1:(i-1)])
      # k <- length(ids_un)
      pr <- c(n_k-sigma, alpha+sigma*k)
      pr <- pr/sum(pr)
      cl <- sample(x = k+1, size = 1, prob = pr)
      ids[i] <- cl
      if(cl <= k){
        n_k[cl] <- n_k[cl]+1
      } else {
        ids_un <- c(ids_un, cl)
        n_k <- c(n_k,1)
        k <- k + 1
      }
    }
  }
  return(ids)
}
stick_breaking <- function(alpha,k){
  v <- rbeta(n = k-1, shape1 = 1, shape2 = alpha)
  p <- c(v,1) * cumprod(1-c(0,v))
  return(p)
}
reorder_part <- function(z){
  z_new <- z
  z_un <- unique(z)
  count <- 1
  for(k in z_un){
    z_new[which(z == k)] <- count
    count <- count + 1
  }
  return(z_new)
}
## using the construction of Teh et al.
stick_breaking_HDP_j <- function(alpha0, beta, k){
  # k should be length(beta)
  temp <- (1 - cumsum(beta))
  temp[temp < 0] <- 0
  v_j <- rbeta(k, alpha0 * beta, alpha0 * temp)
  v_j <- v_j[-k]
  p_j <- c(v_j,1) * cumprod(1-c(0,v_j))
  return(p_j)
}
stick_breaking_HDP <- function(gamma, alpha0, k, nrest){
  beta <- stick_breaking(gamma,k)
  pis <- array(NA, dim = c(nrest, k))
  for(j in 1:nrest){
    pis[j,] <- stick_breaking_HDP_j(alpha0, beta, k)
  }
  return(pis)
}

find_unique <- function(prop){
  nrest <- nrow(prop)
  new_rest <- numeric(nrest)
  new_prop <- c()
  for(i in 1:nrest){
    p <- prop[i,]
    if(i == 1){
      new_rest[i] <- max(new_rest) + 1
      new_prop <- rbind(new_prop, p)
    } else {
      equal <- which(apply(new_prop, 1, function(x) all.equal(x, p)) == TRUE)
      if(length(equal) == 1){
        new_rest[i] <- as.numeric(equal)# new_rest[equal]
      } else {
        new_rest[i] <- max(new_rest) + 1
        new_prop <- rbind(new_prop, p)
      }
    }
  }
  return(new_rest)
}

find_unique_dist <- function(prop, eps){
  nrest <- nrow(prop)
  new_rest <- numeric(nrest)
  new_prop <- c()
  for(i in 1:nrest){
    p <- prop[i,]
    if(i == 1){
      new_rest[i] <- max(new_rest) + 1
      new_prop <- rbind(new_prop, p)
    } else {
      TVdist <- rowSums(abs(t(t(new_prop) - p)))/2
      similar <- which.min(TVdist)
      if(TVdist[similar] > eps){
        similar <- NULL
      }
      if(length(similar) == 1){
        new_rest[i] <- as.numeric(similar)# new_rest[equal]
      } else {
        new_rest[i] <- max(new_rest) + 1
        new_prop <- rbind(new_prop, p)
      }
    }
  }
  return(new_rest)
}

log_pyp_prior <- function(part, alpha, sigma){
  # for each cluster we need (1-sigma)_{n_i-1} = Gamma(n_i - sigma)/Gamma(1-sigma)
  n <- length(part)
  x <- as.numeric(table(part))
  l <- sum(lgamma(x-sigma)-lgamma(1-sigma))
  l <- l + sum(log(alpha + (0:(length(x)-1))*sigma))
  l <- l - sum(log(alpha + 0:(n-1)))
  return(l)
}

# parts_mat is a (d, k) matrix (d rows, k columns)
# with the number of elements in restaurant i, in cluster j

log_hdp_prior <- function(parts_mat, c, c0){
  d <- nrow(parts_mat)
  k <- ncol(parts_mat)
  Ns <- rowSums(parts_mat) # number of elements in each restaurant
  # abs(Stirling1(n,k)) # log(abs(Stirling1.all(200)))
  # Ls_mat has same dimension as parts_mat, but contains the Lij
  lprob <- k * log(c0) - sum(lgamma(c + Ns) - lgamma(c))
  # we should have the sum over Ls_mat
  starting_Ls_mat <- parts_mat*0 + 1
  if(any(parts_mat == 0)){
    tmp <-which(parts_mat==0, arr.ind = T)
    starting_Ls_mat[tmp] <- 0
  }
  lprob <- lprob + big_sum(starting_Ls_mat, parts_mat)
  
  return(lprob)
}

big_sum <- function(Ls_mat, parts_mat){
  if(any(Ls_mat > parts_mat)){
    return(0)
  } else {
    lprob <- sum(Ls_mat) * log(c) - (lgamma(c0 + sum(Ls_mat)) - lgamma(c0))
    lprob <- lprob + sum(lgamma(colSums(Ls_mat)))
    lprob <- lprob + sum(apply(cbind(as.numeric(parts_mat), as.numeric(Ls_mat)), MARGIN = 1, function(x) log(abs(Stirling1(x[1],x[2])))))
    s <- sum(sapply(1:length(Ls_mat), function(i){add_this <- Ls_mat*0; add_this[i] <- 1; big_sum(Ls_mat + add_this, parts_mat)}))
    return(lprob +  s)
  }
}