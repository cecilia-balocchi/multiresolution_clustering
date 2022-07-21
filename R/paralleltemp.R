dyn.load(paste("src/main", .Platform$dynlib.ext, sep = ""))
paralleltemp <- function(Y, YLR = NULL, Mapping, InvMapping, n, n2, burnin, adjust, save, gibbs, 
                     nchains, temps,
                     print_bool = FALSE, sampleLR_bool = TRUE, sampleHR_bool = TRUE,
                     one_vs_manyLR_bool = TRUE, one_vs_manyHR_bool = TRUE,
                     eta_LR = 1, eta_CT = 1, eta_TD = 1, 
                     kH = 0.1, alphaH = 1, betaH = 1, 
                     kL = 0.1, alphaL = 1, betaL = 1, seed = NULL, 
                     grid_eta = NULL, prior_eta = NULL,
                     etaLR_bool = TRUE, etaCT_bool = TRUE, etaTD_bool = TRUE,
                     gammaLR_init = NULL, gammaHR_init = NULL,
                     alpha1 = 1.0, alpha2 = NULL){
  #' Run the nHDP MCMC algorithm using parallel tempering
  #' 
  #' @param Y A vector of high resolution data to be partitioned. Each observation corresponds to a HR unit.
  #' @param YLR A vector of low resolution data to be partitioned. Not required for the clustering to work. Each observation corresponds to a LR unit. If provided, a second likelihood term will be considered for the LR data, and it will be weighted using alpha2.
  #' @param Mapping A two column array where each row corresponds to (i_HR, i_LR), where i_HR is the index of a HR unit and i_LR is the index of the LR unit in which the former is contained. Should be 0-indexed.  
  #' @param InvMapping A two column array where the j-th row corresponds to (r_HR, j_LR), where r_HR corresponds to the row of Mapping such that Mapping[r_HR,1] = j, and j_LR is the index of the LR unit containing i.
  #' @param n An integer, specifying the length of Y, i.e. the number of high resolution units
  #' @param n2 An integer, specifying the length of YLR, i.e. the number of low resolution units
  #' @param burnin An integer, specifying how many initial MCMC samples should be discarded as burn-in
  #' @param adjust Discarded, should be equal to zero.
  #' @param save An integer, specifying how many MCMC iterations should be run and saved after the burn-in stage.
  #' @param gibbs An integer, specifying how many gibbs iteration should be performed for the Split-Merge Metropolis Hasting proposal
  #' @param nchains An integer, specifying how many chains should be ran for the parallel tempering algorithm.
  #' @param temps A vector of nchains positive temperatures (greater 0 and less than or equal than 1) for the parallel tempering algorithm. Should be ordered to be increasing.
  #' @param print_bool A boolean, specifying if additional messages should be printed. Default is FALSE.
  #' @param sampleLR_bool A boolean, specifying if the LR partition should be sampled or fixed. Default is TRUE (should be sampled).
  #' @param sampleHR_bool A boolean, specifying if the HR partition should be sampled or fixed. Default is TRUE (should be sampled).
  #' @param one_vs_manyLR_bool A boolean, specifying if the LR partition should be initialized with one cluster (TRUE) or n2 clusters (FALSE). Default is TRUE.
  #' @param one_vs_manyHR_bool A boolean, specifying if the HR partition should be initialized with one cluster (TRUE) or n2 clusters (FALSE). Default is TRUE.
  #' @param eta_LR A positive number, specifying the concentration hyper-parameter for the DP prior distribution on the LR parition. Default is 1.
  #' @param eta_CT A positive number, specifying the concentration hyper-parameter for the HDP prior distribution on the HR parition. CT stands for 'costumers-tables', i.e. the first level of the HDP distribution. Default is 1.
  #' @param eta_TD A positive number, specifying the concentration hyper-parameter for the HDP prior distribution on the HR parition. TD stands for 'tables-dishes', i.e. the second level of the HDP distribution. Default is 1.
  #' @param kH A positive number, specifying a hyper-parameter for the base measure of the nHDP, i.e. k_0 in equation (1) of the paper. Default is 0.1.
  #' @param alphaH A positive number, specifying the first hyper-parameter for the Inv-Gamma prior on sigma^2, i.e. beta_0 in equation (1) of the paper. Default is 1.
  #' @param betaH A positive number, specifying the second hyper-parameter for the Inv-Gamma prior on sigma^2, i.e. beta_1 in equation (1) of the paper. Default is 1.
  #' @param kL Discarded. A positive number, specifying a hyper-parameter for the low resolution base measure of the nHDP, when YLR is provided. Default is 0.1.
  #' @param alphaL Discarded. A positive number, specifying the first hyper-parameter for the  Inv-Gamma prior on the low resolution sigma^2, when YLR is provided. Default is 1.
  #' @param betaL Discarded. A positive number, specifying the second hyper-parameter for the  Inv-Gamma prior on the low resolution sigma^2, when YLR is provided. Default is 1.
  #' @param seed A positive number, specifying the seed to run the algorith. Default is NULL.
  #' @param grid_eta A vector of positive numbers. The grid of values among which the various eta hyper-parameters should be sampled. Default is NULL. If not provided, the eta hyper-parameters are considered fixed.
  #' @param prior_eta A vector of prior probabilities assigned to each value of grid_eta (length of grid_eta should be equal to length of prior_eta).
  #' @param etaLR_bool A boolean, specifying if eta_LR should be sampled.
  #' @param etaCT_bool A boolean, specifying if eta_CT should be sampled.
  #' @param etaTD_bool A boolean, specifying if eta_TD should be sampled.
  #' @param gammaLR_init A vector of n2 positive integers, specifying the cluster assignemnts of a custom partition, with which to initialize the low resolution partition.
  #' @param gammaHR_init A vector of n positive integers, specifying the cluster assignemnts of a custom partition, with which to initialize the high resolution partition.
  #' @param alpha1 A positive number between 0 and 1, to weight (temper) the high resolution data likelihood. Default is 1. Can be speficied if low resolution data is provided.
  #' @param alpha2 Discarded. A positive number between 0 and 1, to weight (temper) the low resolution data likelihood. Default is NULL. Can be speficied if low resolution data is provided.

  if(is.null(grid_eta)){
    # if not specified, it means we do not want to sample eta, then we choose a grid with one value only
    grid_eta = 1
  }
  ngrid = length(grid_eta)
  if(is.null(prior_eta)){
    prior_eta = rep(1, ngrid)/ngrid
  }

  print <- ifelse(print_bool, 1, 0)
  sampleLR <- ifelse(sampleLR_bool, 1, 0)
  sampleHR <- ifelse(sampleHR_bool, 1, 0)
  one_vs_manyLR <- ifelse(one_vs_manyLR_bool, 1, 0)
  one_vs_manyHR <- ifelse(one_vs_manyHR_bool, 1, 0)
  etaLRb <- ifelse(etaLR_bool, 1, 0)
  etaCTb <- ifelse(etaCT_bool, 1, 0)
  etaTDb <- ifelse(etaTD_bool, 1, 0)
  if(is.null(seed)){
    seedp = integer(1)
  } else {
    seedp = as.integer(seed)
  }
  gammapLR = integer(1)
  gammapHR = integer(1)
  KpLR = 0; KpHR = 0
  if(!is.null(gammaLR_init)){
    gammapLR = as.integer(gammaLR_init)
    KpLR = as.integer(length(unique(gammaLR_init)))
  }
  if(!is.null(gammaHR_init)){
    gammapHR = as.integer(gammaHR_init)
    KpHR = as.integer(length(unique(gammaHR_init)))
  }
  if(is.null(YLR)){
    YLRp = integer(1)
  } else {
    YLRp = as.double(YLR)
  }
  if(is.null(alpha2)){
    alpha2p = integer(1)
  } else {
    alpha2p = as.double(alpha2)
  }
  out <- .C("mainfunc", Y_input = as.double(Y), YLR_input = as.double(YLR),
            Mapping_input = as.double(Mapping), InvMapping_input = as.double(InvMapping), n_ptr = as.integer(n), 
            iter_burnin_ptr = as.integer(burnin), iter_adjust_ptr = as.integer(adjust), 
            iter_save_ptr = as.integer(save), gibbs_iter_ptr = as.integer(gibbs),
            print_bool_ptr = as.integer(print), sampleLR_bool_ptr = as.integer(sampleLR), sampleHR_bool_ptr = as.integer(sampleHR),
            one_vs_many_LR_ptr = as.integer(one_vs_manyLR), one_vs_many_HR_ptr = as.integer(one_vs_manyHR),
            eta_LR_ptr = as.double(eta_LR), eta_CT_ptr = as.double(eta_CT), eta_TD_ptr = as.double(eta_TD), 
            ngrid_eta_ptr = as.integer(ngrid), grid_eta_input = as.double(grid_eta), prior_eta_input = as.double(prior_eta),
            sampleEtaLR_bool_ptr = as.integer(etaLRb), sampleEtaCT_bool_ptr = as.integer(etaCTb), sampleEtaTD_bool_ptr = as.integer(etaTDb),
            k_H_ptr = as.double(kH), alpha_H_ptr = as.double(alphaH), beta_H_ptr = as.double(betaH), 
            k_L_ptr = as.double(kL), alpha_L_ptr = as.double(alphaL), beta_L_ptr = as.double(betaL), 
            alpha1_ptr = as.double(alpha1), alpha2_ptr = as.double(alpha2p),
            seed_ptr = seedp, gammaLR_init_ptr = gammapLR, gammaHR_init_ptr = gammapHR, KLR_init_ptr = KpLR, KHR_init_ptr = KpHR,
            accept = double(3 * nchains + nchains - 1), MCMCchain_low = double(n2 * (1 + save)), MCMCchain_high = double(n * (1 + save)),
            MCMCchain_highRest = double(n * (1 + save)), MCMCchain_highTable = double(n * (1 + save)),
            MCMCchain_eta = double(3 * (1 + save)), switch_chain = double(1+save),
            postmean_low = double(n2 * (1 + save)), postmean_high = double(n * (1 + save)),
            nchains_ptr = as.integer(nchains), Ts = as.double(temps))
  if(nchains>1){
    out$accept0 <- out$accept[(3 * nchains + 1):(3 * nchains + nchains - 1)]
  } else {
    out$accept0 <- NULL
  }
  out$accept <- matrix(out$accept[1:(3 * nchains)], nrow = nchains, byrow = TRUE)
  out$lowres <- matrix(out$MCMCchain_low, nrow = save+1, byrow = TRUE)
  out$highres <- matrix(out$MCMCchain_high, nrow = save+1, byrow = TRUE)
  out$highRest <- matrix(out$MCMCchain_highRest, nrow = save+1, byrow = TRUE)
  out$highTable <- matrix(out$MCMCchain_highTable, nrow = save+1, byrow = TRUE)
  out$highResTab <- array(paste0(out$highRest,":",out$highTable), dim = c(save+1,n))
  out$highres_all <- array(paste0(out$highRest,":",out$highTable,":",out$highres), dim = c(save+1,n))
  out$eta <- matrix(out$MCMCchain_eta, nrow = save+1, byrow = TRUE)
  out$lowres_mean <- matrix(out$postmean_low, nrow = save+1, byrow = TRUE)
  out$highres_mean <- matrix(out$postmean_high, nrow = save+1, byrow = TRUE)
  return(out)
}
# Note: sapply(strsplit("0:0:1", ":"), FUN = "[", 2) first splits "0:0:1" into c("0" "0" "1"), then selects the second element



