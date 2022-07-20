rm(list= ls())
tmp <- load("results/summaries/modeldata_all.rdata") # All summaries are combined here

### To replicate the plots in the paper, we should run this for config = "353" and "111"
### n_ct is fixed = 10.
config <- "353"
# config <- "111"

d <- "25"
Nsim <- 50
possible_methods_strs <- c("km_nHDP", "km", "nDP", 
                           "nHDP1_0","nHDP1_0_sampleeta", 
                           "nHDPheur", "HDPoracle", "HDPone", "HDPmany")
used_methods <- c(1,2,3)
models <- possible_methods_strs[used_methods]
mcmc_models <- c(1,3)
non_mcmc_models <- 2

## First do m = 10 (= nBGinCT)
m <- 10
resultsHR <- get(paste0("resultsHR_m",m,"_d",d,"sim",config))
resultsLR <- get(paste0("resultsLR_m",m,"_d",d,"sim",config))

tmp_data <- data.frame(model = NULL, 
                       R2 = NULL, RMSE = NULL, VI_HR = NULL,
                       Radj_HR = NULL, Kratio_HR = NULL,
                       MSE_Gjs = NULL, VI_LR = NULL, 
                       Radj_LR = NULL, Kratio_LR = NULL)

for(i in mcmc_models){
  tmp1 <- data.frame(model = rep(models[i], Nsim), 
                     R2 = resultsHR[models[i], "R2_OS_bayes",],
                     RMSE = sqrt(resultsHR[models[i], "MSE_bayes",]),
                     VI_HR = resultsHR[models[i], "VI_est",],
                     Radj_HR = resultsHR[models[i], "Radj_bayes",],
                     Kratio_HR = resultsHR[models[i], "K_bayes",]/resultsHR[models[i], "K_true",],
                     R2_LR = resultsLR[models[i], "R2_OS_bayes",],
                     RMSE_Gjs = sqrt(resultsLR[models[i], "MSE_Gjs_bayes",]),
                     VI_LR = resultsLR[models[i], "VI_est",],
                     Radj_LR = resultsLR[models[i], "Radj_bayes",],
                     Kratio_LR = resultsLR[models[i], "K_bayes",]/resultsLR[models[i], "K_true",])
  tmp_data <- rbind(tmp_data, tmp1)
}

tmp_data$model[which(tmp_data$model == "km_nHDP")] <- "nHDP"

## for KM it's different
for(i in non_mcmc_models){ 
  tmp1 <- data.frame(model = rep(models[i], Nsim), 
                     R2 = resultsHR[models[i], "R2_OS_est",],
                     RMSE = sqrt(resultsHR[models[i], "MSE_est",]),
                     VI_HR = resultsHR[models[i], "VI_est",],
                     Radj_HR = resultsHR[models[i], "Radj",],
                     Kratio_HR = resultsHR[models[i], "K_est",]/resultsHR[models[i], "K_true",],
                     R2_LR = resultsLR[models[i], "R2_OS_est",],
                     RMSE_Gjs = sqrt(resultsLR[models[i], "MSE_Gjs_est",]),
                     VI_LR = resultsLR[models[i], "VI_est",],
                     Radj_LR = resultsLR[models[i], "Radj",],
                     Kratio_LR = resultsLR[models[i], "K_est",]/resultsLR[models[i], "K_true",])
  tmp_data <- rbind(tmp_data, tmp1)
}

tmp_data$model <- factor(tmp_data$model , levels=c("nHDP","km","nDP"))
tmp_data_n10 <- tmp_data

## Now we do m = 50 (= nBGinCT)
m <- 50
resultsHR <- get(paste0("resultsHR_m",m,"_d",d,"sim",config))
resultsLR <- get(paste0("resultsLR_m",m,"_d",d,"sim",config))

tmp_data <- data.frame(model = NULL, 
                       R2 = NULL, RMSE = NULL, VI_HR = NULL,
                       Radj_HR = NULL, Kratio_HR = NULL,
                       MSE_Gjs = NULL, VI_LR = NULL, 
                       Radj_LR = NULL, Kratio_LR = NULL)

for(i in mcmc_models){
  tmp1 <- data.frame(model = rep(models[i], Nsim), 
                     R2 = resultsHR[models[i], "R2_OS_bayes",],
                     RMSE = sqrt(resultsHR[models[i], "MSE_bayes",]),
                     VI_HR = resultsHR[models[i], "VI_est",],
                     Radj_HR = resultsHR[models[i], "Radj_bayes",],
                     Kratio_HR = resultsHR[models[i], "K_bayes",]/resultsHR[models[i], "K_true",],
                     R2_LR = resultsLR[models[i], "R2_OS_bayes",],
                     RMSE_Gjs = sqrt(resultsLR[models[i], "MSE_Gjs_bayes",]),
                     VI_LR = resultsLR[models[i], "VI_est",],
                     Radj_LR = resultsLR[models[i], "Radj_bayes",],
                     Kratio_LR = resultsLR[models[i], "K_bayes",]/resultsLR[models[i], "K_true",])
  tmp_data <- rbind(tmp_data, tmp1)
}

tmp_data$model[which(tmp_data$model == "km_nHDP")] <- "nHDP"

## for KM it's different
for(i in non_mcmc_models){ 
  tmp1 <- data.frame(model = rep(models[i], Nsim), 
                     R2 = resultsHR[models[i], "R2_OS_est",],
                     RMSE = sqrt(resultsHR[models[i], "MSE_est",]),
                     VI_HR = resultsHR[models[i], "VI_est",],
                     Radj_HR = resultsHR[models[i], "Radj",],
                     Kratio_HR = resultsHR[models[i], "K_est",]/resultsHR[models[i], "K_true",],
                     R2_LR = resultsLR[models[i], "R2_OS_est",],
                     RMSE_Gjs = sqrt(resultsLR[models[i], "MSE_Gjs_est",]),
                     VI_LR = resultsLR[models[i], "VI_est",],
                     Radj_LR = resultsLR[models[i], "Radj",],
                     Kratio_LR = resultsLR[models[i], "K_est",]/resultsLR[models[i], "K_true",])
  tmp_data <- rbind(tmp_data, tmp1)
}

tmp_data$model <- factor(tmp_data$model , levels=c("nHDP","km","nDP"))
tmp_data_n50 <- tmp_data



png(paste0("figures/boxplots_paper_modeldata",config,"_",d,"_text2.png"), width = 8, height = 4, units = "in", res = 200)
par(mar = c(2,3,2,1), mfrow = c(2,4), cex.axis = 1.2, cex.main = 1.3,oma=c(0,3,3,0))
boxplot(tmp_data_n10$RMSE ~ tmp_data_n10$model, main = expression(paste("RMSE ",theta[jk])))
boxplot(tmp_data_n10$VI_HR ~ tmp_data_n10$model, main = expression(paste("VI ",gamma^H)))
boxplot(tmp_data_n10$RMSE_Gjs ~ tmp_data_n10$model, main = expression("RMSE G"[j]))
boxplot(tmp_data_n10$VI_LR ~ tmp_data_n10$model, main = expression(paste("VI ",gamma^L)))

boxplot(tmp_data_n50$RMSE ~ tmp_data_n50$model)#, main = expression(paste("RMSE ",theta[jk])))
boxplot(tmp_data_n50$VI_HR ~ tmp_data_n50$model)#, main = expression(paste("VI ",gamma^H)))
boxplot(tmp_data_n50$RMSE_Gjs ~ tmp_data_n50$model)#, main = expression("RMSE G"[j]))
boxplot(tmp_data_n50$VI_LR ~ tmp_data_n50$model)#, main = expression(paste("VI ",gamma^L)))

mtext(text="High-resolution",side=3,line=0,outer=TRUE, adj = 0.25, font = 2)
mtext(text="Low-resolution",side=3,line=0,outer=TRUE, adj = 0.8, font = 2)
mtext(text=expression(bold(paste("n"["\u2113"]," = 10"))),side=2,line=0,outer=TRUE, adj = 0.8)
mtext(text=expression(bold(paste("n"["\u2113"]," = 50"))),side=2,line=0,outer=TRUE, adj = 0.25)
dev.off()
