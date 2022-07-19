rm(list = ls())

Nsim <- 50
n_ct = 10            # number of low-res units

possible_dist_strs <- c("2", "25", "3", "4")
possible_ms <- c(5, 10, 25, 50)
possible_methods_strs <- c("km_nHDP", "km", "nDP", 
                           "nHDP1_0","nHDP1_0_sampleeta", 
                           "nHDPheur", "HDPoracle", "HDPone", "HDPmany")
possible_etas <- c("353","111")
used_methods <- c(1,2,3)
methods_strs <- possible_methods_strs[used_methods]
used_ms <- c(2,4)
ms <- possible_ms[used_ms]
dist_strs <- possible_dist_strs[c(2,4)]
used_etas <- c(1,2)
etas <- possible_etas[used_etas]

## Create the big array to store all results
configs <- c("km_nHDP", "km", "nDP", 
             "nHDP1_0","nHDP1_0_sampleeta", 
             "nHDPheur", "HDPoracle", "HDPone", "HDPmany")
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
result <- array(NA, dim = c(length(configs), length(colnames_results), Nsim))
colnames(result) <- colnames_results
rownames(result) <- configs

get_string <- function(sim, m, d_str, met, etas){
  c("results/sim/modeldata/nCT",n_ct,"nBGinCT",m,"_etas",etas,"_save",10000,"_dist", d_str,"_met",met,"sim",sim,".rdata")
}
save_list <- c("get_string")

for(i_d in 1:length(dist_strs)){
  for(i_eta in 1:length(etas)){
    for(i_m in 1:length(ms)){
      d_str <- dist_strs[i_d]
      m <- ms[i_m]
      eta_str <- etas[i_eta]
      
      resultsHR <- result
      resultsLR <- result
      for(met in 1:length(methods_strs)){
        actualNsim <- 0
        for(sim in 1:Nsim){
          filestring <- paste0(get_string(sim, m, d_str, used_methods_inarr[met], etas = eta_str), collapse = "")
          if(file.exists(filestring)){
            tmp <- load(filestring)
            
            resultsHR[methods_strs[met],,sim] <- resultHR[methods_strs[met],]
            resultsLR[methods_strs[met],,sim] <- resultLR[methods_strs[met],]
            
            
            actualNsim <- actualNsim + 1
          } else {
            cat(d_str, " ", m, " ", eta_str, " ", sim, " ", methods_strs[met], "\n")
          }
        }
      }
      assign(paste0("resultsHR_m",m,"_d",d_str,"sim",eta_str),resultsHR)
      assign(paste0("resultsLR_m",m,"_d",d_str,"sim",eta_str),resultsLR)
      save_list <- c(save_list, paste0("resultsHR_m",m,"_d",d_str,"sim",eta_str), 
                     paste0("resultsLR_m",m,"_d",d_str,"sim",eta_str))
    }
  }
}

save(list = save_list, file = "results/summaries/modeldata_all.rdata")

