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



