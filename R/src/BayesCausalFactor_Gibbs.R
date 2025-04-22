library(truncnorm)
library(BNPmix)

library(statmod)
library(matlab)
library(R.matlab)
library(MASS)
library(mvtnorm)



causal_rf_control_stand_mix <- function(nrun = 30000, burn = 20000, thin = 1,
                                        nu = 3, nus = 3,
                                        a1 = 1.1, b1 = 1,
                                        a2 = 1.1, b2 = 1,
                                        a1s = 1.1, b1s = 1,
                                        a2s = 2.1, b2s = 1,
                                        apsi = 1, bpsi = 0.3,
                                        abeta = 0, bbeta = 5,
                                        sigma_mu = 1,  #HypPar var for \mu  l_it
                                        mu_alpha = 0, sigma_alpha = 5 #HypPar for \alpha regression (weights)
)
{
  return(list(nrun = nrun, burn = burn, thin = thin,
              nu = nu, nus = nus, a1 = a1, b1 = b1, a2 = a2, b2 = b2,
              a1s = a1s, b1s = b1s, a2s = a2s, b2s = b2s, apsi = apsi, bpsi = bpsi, 
              sigma_mu=sigma_mu, mu_alpha=mu_alpha, sigma_alpha=sigma_alpha,
              abeta = abeta, bbeta = bbeta))
}


causal_rfm_stand_mixfree_reg <- function(Y_t, j_t, X_reg = NULL, X_weight = NULL, R_cluster=10,
                                     trace = TRUE, nprint = 1000,
                                     outputlevel = 1, seed = 1, control = list(...), ...)
{
  set.seed(seed)
  
  if(is.null(X_reg)){
    for (s in 1:2){
      X_reg[[s]] <- ones(1, n_t[s])
    }
  }
  
  if(is.null(X_weight)){
    for (s in 1:2){
      X_weight[[s]] <- ones(1, n_t[s])
    }
  }
  
  #### preparation data ####
  p <- dim(Y_t[[1]])[2]             # dim variables
  d_X_r <- dim(X_reg[[1]])[1]       # dim covariates in the regression
  d_X_w <- dim(X_weight[[1]])[1]    # dim covariates in the mixture weights
  n_t <- c()                        # sample size
  
  #### setting up priors and initialize ####
  control <- do.call("causal_rf_control_stand_mix", control)
  nrun <- control$nrun
  thin <- control$thin
  burn <- control$burn
  sp  <- (nrun - burn) / thin
  a_psi_t <- b_psi_t <-  c()    # gamma hyperparameters for Psi_t
  nu_t <- c()                # gamma hyperpameters for omega_t
  a1_t <- b1_t <- c()      # gamma hyperparameters for delta_1
  a2_t <- b2_t <- c()      #gamma hyperparameters for delta_l^1
  sigma_mu <- c()               #HypPar var for \mu  l_it
  mu_alpha <- c()               #HypPar for \apha regretion - weights (mean)
  sigma_alpha <- c()            #HypPar for \apha regretion - weights (var)
  abeta_t <- bbeta_t <- c()     # hyperparameter for beta (regression parameters)
  
  #### hyper parameters priors ####
  apsi <- control$apsi
  bpsi<- control$bpsi
  nu <- control$nu
  nus <- control$nus
  a1 <- control$a1
  b1 <- control$b1
  a2 <- control$a2
  b2 <-  control$b2
  a1s <- control$a1s
  b1s <- control$b1s
  a2s <- control$a2s
  b2s <-  control$b2s
  sigma_mu <- control$sigma_mu        
  mu_alpha <- control$mu_alpha         
  sigma_alpha <- control$sigma_alpha
  abetat <- control$abeta
  bbetat <- control$bbeta
  
  #### initial values ####
  psi_t <- Psi_t <- Lambda_t <- l_t <- Meta_t <-  Veta_t <-  list()
  omegajh_t <- delta_t <- tauh_t <- Plam_t <- list()
  alpha_weights <- MN_allocation <- pi_prior <- list()      ## mixture: weights
  S_allocation <- n_clusters <- Z_latent <- list()   ## mixture: cluster allocation
  mu_litj <- mu_litj_Sall <- list()     ## mixture: atoms
  Y_centr <- sd_Y <- Y_stand <- Beta_t <- list()
  
  #### output ####
  Lambdaout <- psiout <- l_tout  <- SigmaLambda <- list()
  Y_obs_lt <- Y_mis_lt <- Y_obs_lt_temp <- Y_mis_lt_temp <- list()
  Y_obs_all <- Y_mis_all <- list()
  S_allocation_post <- mu_litj_post <- Beta_tout <- sd_Y_out <- list()
  S_allocation_mis <- l_mis <- l_tout_mis <- list()
  
  # check chains:
  lt_chains <- lambda_lt_chains <- lamda_lt_post <- list()
  
  #### useful functions ----
  # probit stick-breaking weights 
  pi_function=function(X,alpha){
    a_t=rbind(pnorm(alpha%*%X),1)
    compl_at=rbind(1,apply(1-a_t[-R_cluster,],2,cumprod))
    return(a_t*compl_at)
  }
  # gamma for Z variable
  a_weights_function=function(pi_j){
    compl_sum_pi_j=rbind(1,1-apply(pi_j[-R_cluster,],2,cumsum))
    phi_a=pi_j/compl_sum_pi_j
    phi_a[is.infinite(phi_a)]<-0.999
    phi_a[is.na(phi_a)]<-0.999
    phi_a[phi_a<=0] <- 0.001
    phi_a[phi_a>=1] <- 0.999
    return(qnorm(phi_a))
  }
  # point estimate partition
  estimation_partition<-function(S_part){
    # using Variation of Information as loss function
    return(partition.BNPdens(list(clust=t(S_part)),dist = "VI")$partitions[1,])
  }
  
  for (s in 1:2){
    n_t[s] <- nrow(Y_t[[s]])
    a_psi_t[s] <- apsi
    b_psi_t[s] <- bpsi
    nu_t[s] <- nus
    a1_t[s] <- a1s
    b1_t[s] <- b1s
    a2_t[s] <- a2s
    b2_t[s] <- b2s
    abeta_t[s] <-abetat
    bbeta_t[s] <-bbetat
    
    #### initial values ####
    ###for covariance error terms
    psi_t[[s]] <- rgamma(p, shape = a_psi_t[s], scale = 1 / b_psi_t[s])
    Psi_t[[s]] <- diag(1 / psi_t[[s]])
    ###values for study-specific f.l.
    Lambda_t[[s]] <- zeros(p, j_t[s])
    l_t[[s]] <- matrix(rnorm(n_t[s] * j_t[s]), n_t[s], j_t[s])
    Meta_t[[s]] <- zeros(n_t[s], j_t[s])
    Veta_t[[s]] <- eye(j_t[s])
    ###prior for study-specific f.l
    omegajh_t[[s]] <- matrix(rgamma(p * j_t[s], shape = nu_t[s] / 2,
                                    scale = 2 / nu_t[s]), p, j_t[s])
    delta_t[[s]] <- rgamma(j_t[s], shape=c(a1_t[s], rep(a2_t[s], j_t[s] - 1)),
                           scale = c(b1_t[s], rep(b2_t[s], j_t[s] - 1)))
    tauh_t[[s]] <- cumprod(delta_t[[s]])
    Plam_t[[s]] <- matvec(omegajh_t[[s]], tauh_t[[s]])
    ### prior for treatment-specific f.l.: mixture
    S_allocation[[s]] <- matrix(sample(1:R_cluster, n_t[[s]]*j_t[s], replace=TRUE), n_t[[s]], j_t[s])
    #S_allocation[[s]] <- data_allocation[[s]]    ###### remove after check
    n_clusters[[s]] <- matrix(0, ncol=j_t[[s]], nrow=R_cluster)
    for (j in 1:j_t[[s]]){
      temp_table = table(S_allocation[[s]][,j])
      n_clusters[[s]][1:length(temp_table),j] <- temp_table
    }
    Z_latent[[s]] <- array(0, dim=c(n_t[[s]], j_t[s], R_cluster))
    mu_litj[[s]] <- matrix(rnorm(R_cluster*j_t[[s]], 0, sigma_mu), R_cluster, j_t[[s]])
    mu_litj_Sall[[s]] <- sapply(1:j_t[[s]], function(j) mu_litj[[s]][S_allocation[[s]][,j],j])
    alpha_weights[[s]] <- rep(list(matrix(0,R_cluster-1, d_X_w)),j_t[[s]])  #array(0, dim=c(R_cluster-1, d_X_w, j_t[[s]]))
    MN_allocation[[s]] <- list()
    ###beta-covariates regression
    Beta_t[[s]] <- matrix(rnorm(d_X_r*p, abeta_t, bbeta_t), p, d_X_r)
    
    
    #### save output ####
    Y_obs_lt[[s]] <- zeros(n_t[s], p)
    Y_mis_lt[[s]] <- zeros(n_t[s], p)
    Y_obs_lt_temp[[s]] <- zeros(n_t[s], p)
    Y_mis_lt_temp[[s]] <- zeros(n_t[s], p)
    Y_obs_all[[s]] <- array(0, dim=c(n_t[s], p, sp))
    Y_mis_all[[s]] <- array(0, dim=c(n_t[s], p, sp))
      Beta_tout[[s]] <- zeros(p, d_X_r)
      SigmaLambda[[s]] <- zeros(p, p)
      Lambdaout[[s]] <- zeros(p, j_t[s])
      #psiout[[s]] <- zeros(p, 1)
      l_tout[[s]] <- zeros(n_t[s], j_t[s])
      l_tout_mis[[s]] <- zeros(n_t[s], j_t[3-s])
      #check clusters
      S_allocation_post[[s]]  <- array(0, dim=c(sum(n_t), j_t[s], sp))
      mu_litj_post[[s]] <- array(0, dim=c(R_cluster, j_t[s], nrun))
      #n_clusters_post[[s]] <- array(0, dim=c(R_cluster, j_t[s], nrun))
      #check some chains
      lt_chains[[s]] <- array(0, dim=c(5, j_t[s], nrun))
      lambda_lt_chains[[s]] <- array(0, dim=c(5, p, nrun))
      lamda_lt_post[[s]] <- zeros(n_t[s], p)
      Beta_tout[[s]] <- zeros(p, d_X_r)
      sd_Y_out[[s]] <- zeros(p, 1)
  }
  
  ##################################
  # --- start posterior sampling ---
  ##################################
  
  for(r_it in 1:nrun)
  { 
    #### Step 1: treatment-specific latent factors
    l_t2 <- list()
    for (s in 1:2){
      # Y - intercept (X_reg = only 1s)
      Y_centr[[s]] <- Y_t[[s]]-t(Beta_t[[s]]%*%X_reg[[s]])
      
      # standardize Y
      sd_Y[[s]] <- apply(Y_centr[[s]],2,sd)
      Y_stand[[s]] <- matvec(Y_centr[[s]], 1/sd_Y[[s]])
      
      Lmsg1 <- vecmat(psi_t[[s]], Lambda_t[[s]])
      Veta1 <- diag(j_t[s])*1 + t(Lmsg1) %*% Lambda_t[[s]]
      S1 <- solve(qr.R(qr(chol(Veta1))))
      Veta11 <- tcrossprod(S1)
      Mean1_part <- Y_stand[[s]] %*% Lmsg1
      
      if (j_t[[s]]>1){
        Meta1 <- t(sapply(1:n_t[[s]], function(i) 
          (Mean1_part[i,] + diag(mu_litj[[s]][S_allocation[[s]][i,],])/(1))
          %*% Veta11))
        x1 <- matrix(rnorm(n_t[s] * j_t[s]), nrow = n_t[s], ncol = j_t[s])
        l_t[[s]] <- t(sapply(1:n_t[[s]], function(i) 
          Meta1[i,] + x1[i,] %*% t(S1)))
        
      } else {
        Meta1 <- (Mean1_part[,1] + mu_litj[[s]][S_allocation[[s]],1]/1)*(Veta11[1,1])
        x1 <- matrix(rnorm(n_t[s] * j_t[s]), nrow = n_t[s], ncol = j_t[s])
        l_t[[s]] <- Meta1 + x1 * S1[1,1]
      }
      
      l_t2[[s]] <- crossprod(l_t[[s]])
    }
    
    ### Step 2: specific factor loadings, with constraints
    for(s in 1:2){
      for(i in 1:p){
        Qlam <- diag(Plam_t[[s]][i,], j_t[s]) + psi_t[[s]][i] * l_t2[[s]]
        blam <- psi_t[[s]][i] * (t(l_t[[s]]) %*% Y_stand[[s]][, i])
        Llam <- t(chol(Qlam))
        zlam <- rnorm(j_t[s])
        vlam <- forwardsolve(Llam, blam)
        mlam <- backsolve(t(Llam), vlam)
        ylam <- backsolve(t(Llam), zlam)
        zlam4 <- zlam
        Lambda_t[[s]][i,] <- t(ylam + mlam)
      }
    }
    
    ### Step 3: omegajh for Lambda_t
    for(s in 1:2){
      for(h in 1:j_t[s]){
        omegajh_t[[s]][, h] <- rgamma(p, shape= (nu_t[[s]] + 1) / 2,
                                      rate = (nu_t[[s]] + tauh_t[[s]][h] * Lambda_t[[s]][,h]^2) / 2)
      }
      mat_t <- omegajh_t[[s]] * Lambda_t[[s]]^2
      ad <- a1_t[s] + 0.5 * p*j_t[s]
      bd <- b1_t[s] + 0.5 * (1 / delta_t[[s]][1]) * sum(tauh_t[[s]] * colSums(mat_t))
      delta_t[[s]][1] <- rgamma(1, shape = ad, scale = 1 / bd)
      tauh_t[[s]] <- cumprod(delta_t[[s]])
      if (j_t[s]>1){
        for(h in 2:j_t[s]){
          ad <- a2_t[s] + 0.5 * p * (j_t[s] - h + 1)
          bd <- b2_t[s] + 0.5 * (1 / delta_t[[s]][h]) * sum(tauh_t[[s]][h:j_t[s]] * colSums(as.matrix(mat_t[, h:j_t[s]])))
          delta_t[[s]][h] <- rgamma(1, shape = ad, scale = 1 / bd)
          tauh_t[[s]] <- cumprod(delta_t[[s]])
        }
      }
      Plam_t[[s]] <- matvec(omegajh_t[[s]], tauh_t[[s]])
    }
    
    ### Step 4: Psi_t
    for (s in 1:2){
      Ytil_t <- Y_stand[[s]] - (l_t[[s]] %*% t(Lambda_t[[s]]))
      psi_t[[s]] <- rgamma(p, shape = a_psi_t[s] + 0.5 * n_t[s], rate = b_psi_t[s] + 0.5 * colSums(Ytil_t^2))
      Psi_t[[s]] <- diag(1 / psi_t[[s]])
    }
    
    ### ---- Step 5: cluster allocation and weights estimation ----
    #5a: cluster-specific parameters (mu) ----
    for (s in 1:2){
      M_mu <- sapply(1:j_t[[s]], function(j) sapply(1:R_cluster, function(cl)
        sum(l_t[[s]][S_allocation[[s]][,j]==cl,j])/(2/R_cluster) ))
      M_mu[is.na(M_mu)] <- 0
      V_mu_solve <- 1/((n_clusters[[s]])/(2/R_cluster) + 1/sigma_mu)
      mu_litj[[s]] <- matrix(rnorm(R_cluster*j_t[[s]], c(V_mu_solve*M_mu), sqrt(c(V_mu_solve))), 
                             R_cluster, j_t[[s]])
      mu_litj_Sall[[s]] <- sapply(1:j_t[[s]], function(j) mu_litj[[s]][S_allocation[[s]][,j],j])
    }
    
    # 5b: cluster allocation -----
    for (s in 1:2){
      pi_prior[[s]] <- lapply(1:j_t[[s]], function(j) pi_function(X_weight[[s]],alpha_weights[[s]][[j]]))
      
      prob_temp <- lapply(1:j_t[[s]], function(j) 
        sapply(l_t[[s]][,j], function(i) 
          dnorm(i, mu_litj[[s]][,j], 1)))
      for(j in 1:j_t[[s]]){
        prob_temp[[j]][(prob_temp[[j]])<=0] <- 0.000001
      }
      
      MN_allocation[[s]] <- lapply(1:j_t[[s]], function(j) 
        sapply(1:n_t[[s]], function(i) 
          t(rmultinom(1,1, prob_temp[[j]][,i]*(pi_prior[[s]][[j]][,i])))))
      
      S_allocation[[s]] <- sapply(MN_allocation[[s]], function(i) matrix(1:R_cluster,nrow=1)%*%(i))
      n_clusters[[s]] <- sapply(MN_allocation[[s]], function(i) apply(i, 1, sum))
      
    }
    
    #5c: augmentation scheme ----
    for (s in 1:2){
      a_mean_Z <- lapply(1:j_t[[s]], function(j) a_weights_function(pi_prior[[s]][[j]]))
      
      Z_latent[[s]] <- lapply(1:j_t[[s]], function(j) 
        sapply(1:n_t[[s]], function(i) c(rtruncnorm(S_allocation[[s]][i,j]-1,b=0,
                                                    mean=a_mean_Z[[j]][1:(S_allocation[[s]][i,j]-1),i])*(S_allocation[[s]][i,j]>1),
                                         rtruncnorm(1,a=0, mean=a_mean_Z[[j]][S_allocation[[s]][i,j],i]), 
                                         rep(NA,R_cluster-S_allocation[[s]][i,j])) ))
    }
    
    #5d: weights parameters ----
    for (s in 1:2){
      Var_alpha_solve <- lapply(1:j_t[[s]], function(j)
        lapply(1:(R_cluster-1), function(r)
          solve(tcrossprod(matrix(X_weight[[s]][,S_allocation[[s]][,j]>=r], nrow=d_X_w))+diag(d_X_w)/sigma_alpha)))
      
      alpha_weights[[s]] <- lapply(1:j_t[[s]], function(j)
        t(sapply(1:(R_cluster-1), function(r)
          rmvnorm(1, Var_alpha_solve[[j]][[r]]%*%(t(t(X_weight[[s]][,S_allocation[[s]][,j]>=r]))%*%(Z_latent[[s]][[j]][r,S_allocation[[s]][,j]>=r])
                                                  +matrix(mu_alpha,nrow=d_X_w)/sigma_alpha),
                  Var_alpha_solve[[j]][[r]]))))
    }
    
    ### Step 6: beta regression
    for(s in 1:2){
      Y_star <- Y_t[[s]] - (l_t[[s]])%*%t(vecmat(sd_Y[[s]], Lambda_t[[s]]))
      prod_X <- (X_reg[[s]])%*%t(X_reg[[s]])
      
      V_beta <- lapply(1:p, function(i) chol2inv(chol(prod_X*(psi_t[[s]][i])/(sd_Y[[s]][i]^2) + (1/bbeta_t[s])*diag(d_X_r))))
      M_beta <- abeta_t[s]/bbeta_t[s] + (X_reg[[s]])%*%Y_star*(matrix(psi_t[[s]]/(sd_Y[[s]]^2),ncol=p, nrow=d_X_r, byrow=TRUE))
      
      xP1 <- matrix(rnorm(d_X_r*p), nrow = d_X_r, ncol = p)
      Beta_t[[s]] <- matrix(sapply(1:p, function(i) M_beta[,i]%*%V_beta[[i]] + xP1[,i] %*% chol(V_beta[[i]])),
                            ncol=d_X_r, nrow=p, byrow=TRUE)
    }
    
    ### Step IMPUTATION missing data
    for (s in 1:2){
      pi_lmis <- lapply(1:j_t[[s]], function(j) pi_function(X_weight[[s]],alpha_weights[[3-s]][[j]]))
      MN_lmis <- lapply(1:j_t[[s]], function(j) sapply(1:n_t[[s]], function(i) rmultinom(1,1,pi_lmis[[j]][,i])))
      S_allocation_mis[[s]] <- sapply(MN_lmis, function(i) matrix(1:R_cluster,nrow=1)%*%(i))
      
      if (j_t[[3-s]]>1){
        Meta1 <- t(sapply(1:n_t[[s]], function(i) 
          (diag(mu_litj[[3-s]][S_allocation_mis[[s]][i,],])/1)))
        x1 <- matrix(rnorm(n_t[s] * j_t[3-s]), nrow = n_t[s], ncol = j_t[3-s])
        l_mis[[s]] <- t(sapply(1:n_t[[s]], function(i) 
          Meta1[i,] + x1[i,]))
        
      } else {
        Meta1 <- mu_litj[[3-s]][S_allocation_mis[[s]],]
        x1 <- matrix(rnorm(n_t[s] * j_t[3-s]), nrow = n_t[s], ncol = j_t[3-s])
        l_mis[[s]] <- Meta1 + x1
      }
      
    }
    
    #### ---- STORAGE ----
    if(outputlevel==3){
      for(s in 1:2){
        #n_clusters_post[[s]][ , , r_it] <- n_clusters[[s]]
        mu_litj_post[[s]][ , , r_it] <- mu_litj[[s]]
        lt_chains[[s]][ , , r_it] <- l_t[[s]][1:5,]
        lambda_lt_chains[[s]][ , , r_it] <- l_t[[s]][1:5,] %*% t(Lambda_t[[s]])
      }
    }
    
    if(r_it > burn){
      neff <- (r_it - burn) / thin
      teff <- (nrun - burn) / thin
      for (s in 1:2){
        X_beta_obs <- t(Beta_t[[s]]%*%X_reg[[s]])
        X_beta_mis <- t(Beta_t[[3-s]]%*%X_reg[[s]])
        
        Y_obs_lt_temp[[s]] <- X_beta_obs + l_t[[s]] %*% t(vecmat(sd_Y[[s]], Lambda_t[[s]])) + 
          matrix(rnorm(n_t[s]*p),n_t[s],p)*matrix(sd_Y[[s]]/sqrt(psi_t[[s]]),n_t[s],p, byrow = TRUE)
        Y_mis_lt_temp[[s]] <- X_beta_mis + l_mis[[s]] %*% t(vecmat(sd_Y[[3-s]], Lambda_t[[3-s]])) +
          matrix(rnorm(n_t[s]*p),n_t[s],p)*matrix(sd_Y[[3-s]]/sqrt(psi_t[[3-s]]),n_t[s],p, byrow = TRUE)
        
        Y_obs_lt[[s]] <- Y_obs_lt[[s]] + Y_obs_lt_temp[[s]]/teff
        Y_mis_lt[[s]] <- Y_mis_lt[[s]] + Y_mis_lt_temp[[s]]/teff
        
      }
      if(outputlevel==1)
      {
        for(s in 1:2){
          #Lambdaout[[s]][, , neff] <- Lambda_t[[s]]
          #siout[[s]][, , neff] <- 1 / psi_t[[s]]
          #l_tout[[s]] <- l_tout[[s]] + l_t[[s]] / teff
        }
      }
      if(outputlevel==2){
        for(s in 1:2){
          Beta_tout[[s]] <- Beta_tout[[s]] + Beta_t[[s]] / teff
          Y_obs_all[[s]][ , , neff] <- Y_obs_lt_temp[[s]]
          Y_mis_all[[s]][ , , neff] <- Y_mis_lt_temp[[s]]
          SigmaLambda[[s]] <- SigmaLambda[[s]] + tcrossprod(Lambda_t[[s]]) / teff
        }
      }
      if(outputlevel==3){
        for(s in 1:2){
          Lambdaout[[s]] <- Lambdaout[[s]] + Lambda_t[[s]]  / teff
          #psiout[[s]] <- psiout[[s]] + (1 / psi_t[[s]]) / teff
          l_tout[[s]] <- l_tout[[s]] + l_t[[s]] / teff
          l_tout_mis[[s]] <- l_tout_mis[[s]] + l_mis[[s]] / teff
          SigmaLambda[[s]] <- SigmaLambda[[s]] + tcrossprod(Lambda_t[[s]]) / teff
          Beta_tout[[s]] <- Beta_tout[[s]] + Beta_t[[s]] / teff
          lamda_lt_post[[s]] <- lamda_lt_post[[s]] + (l_t[[s]] %*% t(Lambda_t[[s]]))/teff
          sd_Y_out[[s]] <- sd_Y_out[[s]] + sd_Y[[s]]/teff
        }
      }
    }
    if (trace & r_it %% nprint == 0) cat("r=",r_it,"/",nrun,"\n")
  }
  ##### Save and exit
  if(outputlevel == 1)  {
    l_tout <- NULL
    SigmaLambda <- NULL
  }
  if(outputlevel == 2)  {
    Lambdaout <- lt_chains <- NULL
    psiout <- lambda_lt_chains <- sd_Y_out <- NULL
    mu_litj_post <- lamda_lt_post <- NULL
    Lambda <- l_tout <- l_tout_mis <- NULL
  }
  if(outputlevel == 3)  {
    Lambdaout <- NULL
    psiout <- NULL
    Y_obs_all <- Y_mis_all <- NULL
  }
  if(outputlevel == 4)  {
    Lambdaout <- NULL
    psiout <- NULL
    l_tout <- NULL
    SigmaLambda <- NULL
  }
  out <- list(Lambda = Lambdaout, psi = psiout, l_t = l_tout, l_mis = l_tout_mis,
              SigmaLambda = SigmaLambda,
              Y_obs_lt = Y_obs_lt, Y_mis_lt = Y_mis_lt,
              Y_obs_all = Y_obs_all, Y_mis_all = Y_mis_all,
              mu_litj_post = mu_litj_post,
              lt_chains = lt_chains, lambda_lt_chains = lambda_lt_chains,
              lamda_lt_post = lamda_lt_post,
              Beta_tout = Beta_tout, sd_Y_out=sd_Y_out)
  return(structure(out,  class="causal_fa_mix"))
}

