require(compiler)

summarize_nDP <- function(dataset, burnin, skip, save,
              eta_LR = 1, eta_HR = 1,
              seed){
  set.seed(seed)
  J<-max(dataset[,1])
  Nt<-nrow(dataset)
  alpha<-eta_LR                                   # total mass parameter for LR 
  beta<-eta_HR                                    # total mass parameter for HR
  agam<-3                                         # alpha hyperparameter of Inverse gamma distribution (G0)
  bgam<-5                                         # beta hyperparameter of Inverse gamma distribution (G0)
  mu<-0                                           # mu hyperparameter of Inverse gamma distribution (G0)
  lambda<-0.01                                    # lambda hyperparameter of Inverse gamma distribution (G0)
  ni<-2*agam                                      # parameter of student-t distribution
  desvio<-sqrt((2*bgam*(1+lambda))/(ni*lambda))   # parameter of student-t distribution
  Ij<-numeric(J)
  for (i in 1:J) Ij[i]<-sum(dataset[,1]==i)
  #
  # Initialize Sj, Zij and theta
  #
  Sj<-rep(1,J)                                    # distributional cluster membership (LR cluster for LR units)
  Zij<-rep(1,Nt)                                  # observational cluster membership (HR cluster for HR units)
  K<-length(table(Sj))                            # current number of LR clusters
  Lk<-numeric(K)                                  # how many HR clusters are contained in each LR cluster
  tauj<-numeric()
  for (j in 1:length(Sj)) tauj<-c(tauj,rep(Sj[j],Ij[j]))       # LR cluster for HR units
  for (k in 1:K) Lk[k]<-length(table(Zij[which(tauj==k)]))     # number of HR clusters in each LR cluster
  
  theta.samp<-posteriori.theta(dataset, Zij, tauj, lambda, agam, bgam, mu)
  #
  indrejtotal<-0                                     # 0 if accept, 1 if reject
  probacetotal<-1                                    # vector of paccept
  #
  amostrasfin <- save                                  # MCMC sample size after removing burn in and jumps
  # burnin <- 2000                                        # burn in size
  saltos <- skip                                          # jumps size
  AmostrasTotal<-burnin+amostrasfin*saltos           # number of iterations to be run
  #
  enableJIT(3)
  #
  # sample a group from 1 to J
  # sample.tau.xi is the key
  # candidato[[6]] is the new LR cluster value for j (update Sj, tauj)
  # candidato[[7]] is the new HR cluster values for the HR units in j (updates Zij)
  # why sampling again Zij with posteriori.Zij2 ?
  #    I think it's because after the update in LR (which changes HR too), you want to 
  #    separately change HR. However, why doing that only when LR is accepted??

  ## NOTE: theta.samp starts with LR cluster 1, and lists the parameters 
  ## for the HR clusters in 1, ordered by number of HR cluster, and so on.
  ## So if Lk = 1,3,4,5 the first cluster is for LR cluster 1, then cluster 1,2,3 for LR cluster 2 etc.

  ## added by me:
  Sj_samples <- array(NA, dim = c(length(Sj), amostrasfin))
  Zij_samples <- array(NA, dim = c(length(Zij), amostrasfin))
  K_samples <- rep(NA, amostrasfin)
  theta.samp_list <- list()
  Lk_list <- list()

  for (int in 1:AmostrasTotal){
    #
    ######### update Sj and Zij by a MH step of one j
    #
    j<-sample(1:J,1)
    #
    # cat('\n', int, j, K)
    candidato<-sample.tau.xi(dataset, Ij, Sj, Zij, tauj, K, Lk, theta.samp, j, alpha, beta, ni, mu, desvio)
    paccept<-prob.accept(dataset, candidato[[8]], Ij, Sj, Zij, tauj, Lk, theta.samp, j, candidato[[3]], candidato[[2]], candidato[[6]], candidato[[7]], lambda, mu, agam, bgam, candidato[[4]], candidato[[5]], candidato[[10]], candidato[[11]]) [[1]]
    probacetotal<-c(probacetotal,paccept)
    aux2<-runif(1)
    if (aux2<paccept){
      indrejtotal<-c(indrejtotal,0)
      Sj[j]<-candidato[[6]]
      tauj[(cumsum(Ij)[j]-Ij[j]+1):cumsum(Ij)[j]]<-rep(candidato[[6]],Ij[j])
      Zij[(cumsum(Ij)[j]-Ij[j]+1):cumsum(Ij)[j]]<-candidato[[7]]
      while (length(table(Sj))<max(Sj)){ # exclude empty distributional clusters
        categr<-as.numeric(as.character(data.frame(table(Sj))[,1]))
        categd<-seq(1:length(table(Sj)))
        dif<-which(categr!=categd)
        for (i in 1:length(Sj)) if (Sj[i]>dif[1]) Sj[i]<-Sj[i]-1
        tauj<-numeric()
        for (n in 1:length(Sj)) tauj<-c(tauj,rep(Sj[n],Ij[n]))
      } 
      K<-length(table(Sj))
      Lk<-numeric(K)                    # how many HR clusters are contained in each LR cluster
      for (k in 1:K){
        observ<-which(tauj==k)              # which HR units are in LR cluster k?
        Lk[k]<-length(table(Zij[observ]))         # how many HR clusters are contained in each LR cluster
        while (Lk[k]<max(Zij[observ])){           # exclude empty observational clusters by changing the labels in Zij
          categr<-as.numeric(as.character(data.frame(table(Zij[observ]))[,1]))
          categd<-seq(1:Lk[k])
          dif<-which(categr!=categd)
          for (i in observ) if (Zij[i]>dif[1]) Zij[i]<-Zij[i]-1
        }
      }
      theta.samp<-posteriori.theta(dataset, Zij, tauj, lambda, agam, bgam, mu)
      # resamples the Zij for the observations in the cluster now assigned to j
      # in sample.tau.xi updates Zij only for observations in j, not for all the observation with same tauj
      Zij<-posteriori.Zij2(dataset, Zij, tauj, K, theta.samp, ni, mu, desvio, Lk, beta, j, Sj)
      Lk<-numeric(K)
      for (n in 1:K) Lk[n]<-length(table(Zij[which(tauj==n)]))
      theta.samp<-posteriori.theta(dataset, Zij, tauj, lambda, agam, bgam, mu)
    }
    if (aux2>=paccept) indrejtotal<-c(indrejtotal,1)
    #
    if (int>burnin & int%%saltos==0){
      Sj_samples[, (int-burnin)/saltos] <- Sj
      Zij_samples[, (int-burnin)/saltos] <- Zij
      theta.samp_list[[(int-burnin)/saltos]] <- theta.samp
      K_samples[(int-burnin)/saltos] <- K
      Lk_list[[(int-burnin)/saltos]] <- Lk
    }
  }

  # Sj_samples contains in each column the distributional cluster membership (LR cluster for LR units)
  # tauj becomes the LR cluster for HR units
  # Zij_samples contains in each column the observational cluster membership (HR cluster for HR units)
  # Lk contains in each element the number of HR clusters contained in each LR cluster
  # theta.samp_list is (tot_number_HR_cl x 2) (for mu and sigma2)

  lowres <- t(Sj_samples)
  highres <- array(NA, dim = c(amostrasfin, length(Zij)))
  highres_mean <- highres
  highres_sigma2 <- highres
  for(iter in 1:amostrasfin){
    Sj <- Sj_samples[,iter]                             # LR cluster for LR units
    tauj<-numeric(); for (j in 1:length(Sj)) tauj<-c(tauj,rep(Sj[j],Ij[j]))       # LR cluster for HR units
    Zij <- Zij_samples[,iter]                           # HR cluster for HR units
    Lk <- Lk_list[[iter]]                             # number of HR clusters contained in each LR cluster
    theta.samp <- theta.samp_list[[iter]]                     # values of parameters for each HR cluster

    # this gives you unique clusters for the HR units
    # from here you can get the matrix of pairwise clustering, then get the matrix of pairwise probabilities.
    highres_factor <- interaction(tauj, Zij)
    # to sort the interaction clusters correctly (like the theta.samp) you should transform them to characters
    unique_cl <- sort(unique(as.character(highres_factor)))
    highres_num <- match(as.character(highres_factor), unique_cl)

    highres[iter, ] <- highres_num
    #these are the thetas associated to each HR unit
    highres_mean[iter, ] <- theta.samp[highres_num,1]
    highres_sigma2[iter, ] <- theta.samp[highres_num,2]
  }

  return(list(lowres = lowres, highres = highres, highres_mean = highres_mean, highres_sigma2 = highres_sigma2))
}



############
# sample from conditional a posteriori distribution of theta for all k and l
############  
posteriori.theta<-function(dataset, Zij, tauj, lambda, agam, bgam,mu){
  mlk<-aggregate(rep(1,nrow(dataset)), by = list(Zij, tauj), FUN = "sum")[,3]
  Ymeanlk<-aggregate(dataset[,2], by = list(Zij, tauj), FUN = "mean")[,3]
  Yvarlk<-aggregate(dataset[,2], by = list(Zij, tauj), FUN = "var")[,3]*(mlk-1)
  mupost<-((mlk*Ymeanlk)+(lambda*mu))/(lambda+mlk)
  lambdapost<-lambda+mlk
  agampost<-agam+(mlk/2)
  Yvarlk[is.na(Yvarlk)] <- 0
  bgampost<-bgam+(Yvarlk+((mlk*lambda*(Ymeanlk-mu)**2)/(lambda+mlk)))/2
  theta.samp<-rinvgamma(mupost,lambdapost,agampost,bgampost)
  return(theta.samp)
}
#
############
# build a proposal of new cluster for each Sj and realocate its Zij
############
sample.tau.xi<-function(dataset, Ij, Sj, Zij, tauj, K, Lk, theta.samp, j, alpha, beta, ni, mu, desvio){
  nk<-numeric(K)
  for (k in 1:K) 
    nk[k]<-sum(Sj[-j]==k)
  observ<-which(dataset[,1]==j)
  Yclusterj<-dataset[observ,2]
  Zclusterj<-Zij[observ]
  nk[Sj[j]]<-0 # we force Sj candidate to be different of the current value
  nk<-c(nk,alpha)
  probSj<-nk/sum(nk)
  Sprop<-rDiscreta(probSj)
  Sold<-Sj[j]
  Zprop<-numeric(Ij[j])
  Zold<-numeric(Ij[j])
  lprior<-0
  lfunctrans<-0
  lpriorold<-0
  lfunctransold<-0
  taujcand<-tauj
  taujcand[(cumsum(Ij)[j]-Ij[j]+1):cumsum(Ij)[j]]<-rep(K+1,Ij[j])
  nj<-sum(Sj[-j]==Sold)
  if (Sprop <= K){
    for (i in 1:Ij[j]){
      l<-1:Lk[Sprop]
      Zvector<-c(Zij[which(tauj==Sprop)],Zprop[which(Zprop<=Lk[Sprop] & Zprop>0)])
      mlk<-table(Zvector)
      probZij<-mlk*dnorm(Yclusterj[i],theta.samp[(cumsum(Lk)[Sprop]-Lk[Sprop]+l),1],sqrt(theta.samp[(cumsum(Lk)[Sprop]-Lk[Sprop]+l),2]))
      if (length(which(Zprop>Lk[Sprop]))>0){
        mlk2<-table(Zprop[which(Zprop>Lk[Sprop])])
        probZij2<-mlk2*dstudentt(Yclusterj[i],ni,mu,desvio)
        mlk<-c(mlk,mlk2)
        probZij<-c(probZij,probZij2)
      }
      mlk<-c(mlk,beta)
      probZij<-c(probZij,mlk[length(mlk)]*dstudentt(Yclusterj[i],ni,mu,desvio))   
      probZij<-probZij/sum(probZij)
      mlk<-mlk/sum(mlk)
      Zprop[i]<-rDiscreta(probZij)
      lprior<-lprior+log(mlk[Zprop[i]])
      lfunctrans<-lfunctrans+log(probZij[Zprop[i]])
      if (nj > 0){
        l<-1:Lk[Sold]
        Zvector<-c(Zij[which(taujcand==Sold)],Zold[which(Zold<=Lk[Sold] & Zold>0)])
        mlk<-rep(0,Lk[Sold])
        need<-data.frame(table(Zvector))
        categr<-as.numeric(as.character(need[,1]))
        mlk[categr]<-as.numeric(need[,2])
        probZij<-mlk*dnorm(Yclusterj[i],theta.samp[(cumsum(Lk)[Sold]-Lk[Sold]+l),1],sqrt(theta.samp[(cumsum(Lk)[Sold]-Lk[Sold]+l),2]))
        if (length(which(mlk==0))==0){
          mlk<-c(mlk,beta)
          probZij<-c(probZij,mlk[length(mlk)]*dstudentt(Yclusterj[i],ni,mu,desvio))
        } else {
          empty<-which(mlk==0)
          if (mlk[Zclusterj[i]]==0){
            mlk[Zclusterj[i]]<-beta
            probZij[Zclusterj[i]]<-mlk[Zclusterj[i]]*dstudentt(Yclusterj[i],ni,mu,desvio)
          } else {
            mlk[empty[1]]<-beta
            probZij[empty[1]]<-mlk[empty[1]]*dstudentt(Yclusterj[i],ni,mu,desvio)
          }
        } 
        probZij<-probZij/sum(probZij)
        mlk<-mlk/sum(mlk)
        Zold[i]<-Zclusterj[i]
        lpriorold<-lpriorold+log(mlk[Zclusterj[i]])
        lfunctransold<-lfunctransold+log(probZij[Zclusterj[i]])
      } else {
        Zold[1]<-Zclusterj[1]
        Zvector<-c(Zij[which(taujcand==Sold)],Zold[which(Zold <= Lk[Sold] & Zold>0)])
        mlk<-rep(0,Lk[Sold])
        need<-data.frame(table(Zvector))
        categr<-as.numeric(as.character(need[,1]))
        if(length(categr)>0) 
          mlk[categr]<-as.numeric(need[,2])
        probZij<-mlk*dstudentt(Yclusterj[i],ni,mu,desvio)
        if (length(which(mlk==0))==0){
          mlk<-c(mlk,beta)
          probZij<-c(probZij,mlk[length(mlk)]*dstudentt(Yclusterj[i],ni,mu,desvio))
        } else { 
          empty<-which(mlk==0)
          if (mlk[Zclusterj[i]]==0){
            mlk[Zclusterj[i]]<-beta
            probZij[Zclusterj[i]]<-mlk[Zclusterj[i]]*dstudentt(Yclusterj[i],ni,mu,desvio)
          } else {
            mlk[empty[1]]<-beta
            probZij[empty[1]]<-mlk[empty[1]]*dstudentt(Yclusterj[i],ni,mu,desvio)
          }
        } 
        probZij<-probZij/sum(probZij)       
        mlk<-mlk/sum(mlk)
        Zold[i]<-Zclusterj[i]
        lpriorold<-lpriorold+log(mlk[Zclusterj[i]])
        lfunctransold<-lfunctransold+log(probZij[Zclusterj[i]])
      }
    }
  }       
  if (Sprop > K){
    Zprop[1]<-1
    if(Ij[j]>1){
      for (i in 2:Ij[j]){
        mlk<-table(Zprop[which(Zprop>0)])
        probZij<-mlk*dstudentt(Yclusterj[i],ni,mu,desvio)
        mlk<-c(mlk,beta)
        probZij<-c(probZij,mlk[length(mlk)]*dstudentt(Yclusterj[i],ni,mu,desvio))   
        probZij<-probZij/sum(probZij)
        mlk<-mlk/sum(mlk)
        Zprop[i]<-rDiscreta(probZij)
        lprior<-lprior+log(mlk[Zprop[i]])
        lfunctrans<-lfunctrans+log(probZij[Zprop[i]])
        if (nj > 0){
          l<-1:Lk[Sold]
          Zvector<-c(Zij[which(taujcand==Sold)],Zold[which(Zold<=Lk[Sold] & Zold>0)])
          mlk<-rep(0,Lk[Sold])
          need<-data.frame(table(Zvector))
          categr<-as.numeric(as.character(need[,1]))
          mlk[categr]<-as.numeric(need[,2])
          probZij<-mlk*dnorm(Yclusterj[i],theta.samp[(cumsum(Lk)[Sold]-Lk[Sold]+l),1],sqrt(theta.samp[(cumsum(Lk)[Sold]-Lk[Sold]+l),2]))
          if (length(which(mlk==0))==0){
            mlk<-c(mlk,beta)
            probZij<-c(probZij,mlk[length(mlk)]*dstudentt(Yclusterj[i],ni,mu,desvio))
          } else {
            empty<-which(mlk==0)
            if (mlk[Zclusterj[i]]==0){
              mlk[Zclusterj[i]]<-beta
              probZij[Zclusterj[i]]<-mlk[Zclusterj[i]]*dstudentt(Yclusterj[i],ni,mu,desvio)
            } else {
              mlk[empty[1]]<-beta
              probZij[empty[1]]<-mlk[empty[1]]*dstudentt(Yclusterj[i],ni,mu,desvio)
            }
          } 
          probZij<-probZij/sum(probZij)
          mlk<-mlk/sum(mlk)
          Zold[i]<-Zclusterj[i]
          lpriorold<-lpriorold+log(mlk[Zclusterj[i]])
          lfunctransold<-lfunctransold+log(probZij[Zclusterj[i]])
        } else {
          Zold[1]<-Zclusterj[1]
          Zvector<-c(Zij[which(taujcand==Sold)],Zold[which(Zold <= Lk[Sold] & Zold>0)])
          mlk<-rep(0,Lk[Sold])
          need<-data.frame(table(Zvector))
          categr<-as.numeric(as.character(need[,1]))
          if(length(categr)>0) 
            mlk[categr]<-as.numeric(need[,2])
          probZij<-mlk*dstudentt(Yclusterj[i],ni,mu,desvio)
          if (length(which(mlk==0))==0){
            mlk<-c(mlk,beta)
            probZij<-c(probZij,mlk[length(mlk)]*dstudentt(Yclusterj[i],ni,mu,desvio))
          } else { 
            empty<-which(mlk==0)
            if (mlk[Zclusterj[i]]==0){
              mlk[Zclusterj[i]]<-beta
              probZij[Zclusterj[i]]<-mlk[Zclusterj[i]]*dstudentt(Yclusterj[i],ni,mu,desvio)
            } else {
              mlk[empty[1]]<-beta
              probZij[empty[1]]<-mlk[empty[1]]*dstudentt(Yclusterj[i],ni,mu,desvio)
            }
          } 
          probZij<-probZij/sum(probZij)       
          mlk<-mlk/sum(mlk)
          Zold[i]<-Zclusterj[i]
          lpriorold<-lpriorold+log(mlk[Zclusterj[i]])
          lfunctransold<-lfunctransold+log(probZij[Zclusterj[i]])
        }
      }
    }
    
  }
  Sjcand<-Sj
  Sjcand[j]<-Sprop
  taujcand<-tauj
  taujcand[(cumsum(Ij)[j]-Ij[j]+1):cumsum(Ij)[j]]<-rep(Sprop,Ij[j])
  Zijcand<-Zij
  if(length(Zprop) != length((cumsum(Ij)[j]-Ij[j]+1):cumsum(Ij)[j])){
    warning("length(Zprop) different from Ij[j] in sample.tau.xi.")
  }
  Zijcand[(cumsum(Ij)[j]-Ij[j]+1):cumsum(Ij)[j]]<-Zprop
  return(list(Sjcand, taujcand, Zijcand, lprior, lfunctrans, Sprop, Zprop, Yclusterj, Zclusterj,lpriorold,lfunctransold))
}
#
############
# sample from conditional a posteriori distribution of Zij for just one j
############  
posteriori.Zij2<-function(dataset, Zij, tauj, K, theta.samp, ni, mu, desvio, Lk, beta, j, Sj){
  k<-Sj[j]
  observ<-which(tauj==k)
  Yclusterk<-dataset[observ,2]
  Zclusterk<-Zij[observ]
  cont<-1
  for (i in observ){
    probZij<-numeric(Lk[k])
    contmlk<-Zclusterk[-cont]
    l<-seq(1:Lk[k])
    mlk<-rep(0,Lk[k])
    if(sum(contmlk<=Lk[k]) > 0){
      need<-data.frame(table(contmlk[which(contmlk<=Lk[k])])) ## here we get error if no one is true in which.
      # they want to exclude the new clusters that have been created in this function (Zij[i]<-rDiscreta(probZij) if a new cluster is chosen)
      # this happens at the end of observ (if the previous ones all get assigned to new clusters)
      categr<-as.numeric(as.character(need[,1]))
      mlk[categr]<-as.numeric(need[,2])
    } else {# else mlk remains == 0
      warning("avoided disaster.")
    }
    probZij<-mlk*dnorm(Yclusterk[cont],theta.samp[(cumsum(Lk)[k]-Lk[k]+l),1],sqrt(theta.samp[(cumsum(Lk)[k]-Lk[k]+l),2]))
    probZij<-c(probZij,beta*dstudentt(Yclusterk[cont],ni,mu,desvio))    
    probZij<-probZij/sum(probZij)
    Zij[i]<-rDiscreta(probZij)
    Zclusterk<-Zij[observ]
    cont<-cont+1
  }
  while (length(table(Zclusterk))<max(Zclusterk)){ # exclude empty clusters
    categr<-as.numeric(as.character(data.frame(table(Zclusterk))[,1]))
    categd<-seq(1:length(table(Zclusterk)))
    dif<-which(categr!=categd)
    for (i in 1:length(Zclusterk)) 
      if (Zclusterk[i]>dif[1]) 
        Zclusterk[i]<-Zclusterk[i]-1
  }
  Zij[observ]<-Zclusterk
  return(Zij)
}
#
############
# sample from a Normal Inverse gamma distribituion (mu, lambda, alpha, beta)
############
rinvgamma<-function(media.mi,lambda,alpha,beta){
	n<-length(media.mi)
	sigma2<-1/(rgamma(n,alpha,beta))
	mu<-rnorm(n,media.mi,sqrt(sigma2/lambda))
	sample<-cbind(mu,sigma2)
	return(sample)}
#
############
# sample from a Multinomial distribituion (p)
############
rDiscreta<-function(p){
 u<-runif(1)
 P<-cumsum(p)
 val<-sum(P<u)+1
 return(val)}
#
############
# calculate the density function of a student t distribituion (yij, ni, mu, sigma)
############
dstudentt<-function(yij,ni,mu,desvio){
 dens<-(gamma((ni+1)/2)/gamma(ni/2))*((1+(1/ni)*((yij-mu)/desvio)**2)**(-(ni+1)/2))*(1/(sqrt(pi*ni)*desvio))
 return(dens)}
#
############
# calculate the marginalized likelihood of Yj given (Xij ,tauj) and specific l
############
dmarglikeli<-function(Yjl,lambda,agam,bgam,mu){
 mlkj<-length(Yjl)
 media<-mean(Yjl)
 varia<-var(Yjl)*(mlkj-1)
 if (is.na(varia)) varia<-0
 logdens<-(lgamma(mlkj/2+agam)-lgamma(agam))+((-mlkj/2)*log(2*pi))+(0.5*(log(lambda)-log(mlkj+lambda)))+(agam*log(bgam))+((-(mlkj/2+agam))*log(bgam+(varia/2)+((lambda*mlkj*((media-mu)**2))/(2*(mlkj+lambda)))))
 return(logdens)} 
#
############
# calculate the acceptance rate of the pair Sj and all Zij
############  
prob.accept<-function(dataset, Yclusterj, Ij, Sj, Zij, tauj, Lk, theta.samp, j, Zijcand, taujcand, Sprop, Zprop, lambda, mu, agam, bgam, lprior, lfunctrans, lpriorold, lfunctransold){
  # computes the likelihood of old (current) and proposal
  # also puts everything together and computes acceptance
  marglikel<-0
  marglikelold<-0
  #
  Yclusterk<-dataset[which(tauj==Sprop),2]
  Zclusterk<-Zij[which(tauj==Sprop)]
  Lsprop<-length(table(c(Zclusterk,Zprop)))
  if (length(Zclusterk)==0){
    mlk<-rep(0,Lsprop)
    Ymeanlk<-rep(0,Lsprop)
    Yvarlk<-rep(0,Lsprop)}
  if (length(Zclusterk)>0){
    mlk<-table(Zclusterk)
    Ymeanlk<-aggregate(Yclusterk, by = list(Zclusterk,rep(1,length(Zclusterk))), FUN = "mean")[,3]
    Ymeanlk[is.na(Ymeanlk)] <- 0
    Yvarlk<-aggregate(Yclusterk, by = list(Zclusterk,rep(1,length(Zclusterk))), FUN = "var")[,3]*(mlk-1)
    Yvarlk[is.na(Yvarlk)] <- 0
    dif<-Lsprop-length(mlk)
    if (dif > 0){
      mlk<-c(mlk,rep(0,dif))
      Ymeanlk<-c(Ymeanlk,rep(0,dif))
      Yvarlk<-c(Yvarlk,rep(0,dif))}}
  mupost<-((mlk*Ymeanlk)+(lambda*mu))/(lambda+mlk)
  lambdapost<-lambda+mlk
  agampost<-agam+(mlk/2)
  bgampost<-bgam+(Yvarlk+((mlk*lambda*(Ymeanlk-mu)**2)/(lambda+mlk)))/2
  for (l in 1:Lsprop){
    Yclusterjl<-Yclusterj[Zprop==l]
    if (length(Yclusterjl)>0) marglikel<-marglikel+dmarglikeli(Yclusterjl,lambdapost[l],agampost[l],bgampost[l],mupost[l])}
  #
  Yclusterk<-dataset[which(taujcand==Sj[j]),2]
  Zclusterk<-Zijcand[which(taujcand==Sj[j])]
  if (length(Zclusterk)==0){
    mlk<-rep(0,Lk[Sj[j]])
    Ymeanlk<-rep(0,Lk[Sj[j]])
    Yvarlk<-rep(0,Lk[Sj[j]])}
  if (length(Zclusterk)>0){
    mlk<-table(Zclusterk)
    Ymeanlk<-aggregate(Yclusterk, by = list(Zclusterk,rep(1,length(Zclusterk))), FUN = "mean")[,3]
    Ymeanlk[is.na(Ymeanlk)] <- 0
    Yvarlk<-aggregate(Yclusterk, by = list(Zclusterk,rep(1,length(Zclusterk))), FUN = "var")[,3]*(mlk-1)
    Yvarlk[is.na(Yvarlk)] <- 0
    dif<-Lk[Sj[j]]-length(mlk)
    if (dif > 0){
      mlk<-c(mlk,rep(0,dif))
      Ymeanlk<-c(Ymeanlk,rep(0,dif))
      Yvarlk<-c(Yvarlk,rep(0,dif))}}
  mupost<-((mlk*Ymeanlk)+(lambda*mu))/(lambda+mlk)
  lambdapost<-lambda+mlk
  agampost<-agam+(mlk/2)
  bgampost<-bgam+(Yvarlk+((mlk*lambda*(Ymeanlk-mu)**2)/(lambda+mlk)))/2
  Zold<-Zij[(cumsum(Ij)[j]-Ij[j]+1):cumsum(Ij)[j]]
  for (l in 1:Lk[Sj[j]]){
    Yclusterjl<-Yclusterj[Zold==l]
    if (length(Yclusterjl)>0) marglikelold<-marglikelold+dmarglikeli(Yclusterjl,lambdapost[l],agampost[l],bgampost[l],mupost[l])} 
  #
  paccept<-exp(marglikel+lprior+lfunctransold-(marglikelold+lpriorold+lfunctrans))
  return(list(paccept,lprior,lpriorold,lfunctrans,lfunctransold,marglikel,marglikelold))
}
############
# sample from conditional a posteriori distribution of Zij (for all distributional clusters) not used
############	
posteriori.Zij<-function(dataset, Zij, tauj, K, theta.samp, ni, mu, desvio, Lk, beta){
	for (k in 1:K){
		observ<-which(tauj==k)
		Yclusterk<-dataset[observ,2]
		Zclusterk<-Zij[observ]
		cont<-1
		for (i in observ){
			probZij<-numeric(Lk[k])
			contmlk<-Zclusterk[-cont]
			l<-seq(1:Lk[k])
			mlk<-rep(0,Lk[k])
			need<-data.frame(table(contmlk[which(contmlk<=Lk[k])]))
			categr<-as.numeric(as.character(need[,1]))
			mlk[categr]<-as.numeric(need[,2])
			probZij<-mlk*dnorm(Yclusterk[cont],theta.samp[(cumsum(Lk)[k]-Lk[k]+l),1],sqrt(theta.samp[(cumsum(Lk)[k]-Lk[k]+l),2]))
			probZij<-c(probZij,beta*dstudentt(Yclusterk[cont],ni,mu,desvio))		
			probZij<-probZij/sum(probZij)
			Zij[i]<-rDiscreta(probZij)
			Zclusterk<-Zij[observ]
			cont<-cont+1}
		while (length(table(Zclusterk))<max(Zclusterk)){ # exclude empty clusters
			categr<-as.numeric(as.character(data.frame(table(Zclusterk))[,1]))
			categd<-seq(1:length(table(Zclusterk)))
			dif<-which(categr!=categd)
			for (i in 1:length(Zclusterk)) if (Zclusterk[i]>dif[1]) Zclusterk[i]<-Zclusterk[i]-1}
		Zij[observ]<-Zclusterk}
	return(Zij)}
#





