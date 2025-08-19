control_function <- function(nrun = 10000, burn = 5000, thin = 1,
                            nu = 3, nus = 3,
                            a1 = 1.1, b1 = 1,
                            a2 = 1.1, b2 = 1,
                            a1s = 1.1 , b1s = 1,
                            a2s = 2.1, b2s = 1,
                            apsi = 1, bpsi = 0.3,
                            abeta = 0, bbeta = 5 #HypPar for Covariate regression
                            )
{
  return(list(nrun = nrun, burn = burn, thin = thin,
              nu = nu, nus = nus, a1 = a1, b1 = b1, a2 = a2, b2 = b2,
              a1s = a1s, b1s = b1s, a2s = a2s, b2s = b2s, apsi = apsi, bpsi = bpsi,
              abeta = abeta, bbeta = bbeta))
}


StandardFactor <- function(Y_t, X_t, j_t, trace = TRUE, nprint = 1000,
                             outputlevel = 1, seed = 1, control = list(...), ...)
{
  set.seed(seed)
  
  #### preparation data ####
  p <- dim(Y_t[[1]])[2]             # dim variables
  n_t <- c()                        # sample size
  X_t[[1]] <- rbind(1,X_t[[1]])
  X_t[[2]] <- rbind(1,X_t[[2]])
  d_cov <- dim(X_t[[1]])[1]         # dimension of covariates 
  
  #### setting up priors and initialize ####
  control <- do.call("control_function", control)
  nrun <- control$nrun
  thin <- control$thin
  burn <- control$burn
  sp  <- (nrun - burn) / thin
  a_psi_t <- b_psi_t <-  c()    # gamma hyperparameters for Psi_t
  nu_t <- c()                   # gamma hyperpameters for omega_t
  a1_t <- b1_t <- c()           # gamma hyperparameters for delta_1
  a2_t <- b2_t <- c()           #gamma hyperparameters for delta_l^1
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
  abetat <- control$abeta
  bbetat <- control$bbeta
  
  #### initial values ####
  psi_t <- Psi_t <- Lambda_t <- l_t <- Meta_t <-  Veta_t <-  list()
  omegajh_t <- delta_t <- tauh_t <- Plam_t <- list()
  Y_centr <- Beta_t <- list()
  
  #### output ####
  Lambdaout <- psiout <- l_tout  <- SigmaLambda <- list()
  Y_obs_lt <- Y_mis_lt <- list()
    #Y_obs_var <- Y_mis_var <- list()
  Beta_tout <- Lambda_l_out <- list()
  l_mis <- l_tout_mis <- list()
  Y_obs_all <- Y_mis_all <- list()
  
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
    ###beta-covariates regression 
    Beta_t[[s]] <- matrix(rnorm(d_cov*p, abeta_t, bbeta_t), p, d_cov)
    
    #### save output ####
    Y_obs_lt[[s]] <- array(0, dim=c(n_t[s], p, sp))
    Y_mis_lt[[s]] <- array(0, dim=c(n_t[s], p, sp))
    #Y_obs_var[[s]] <- array(0, dim=c(n_t[s], p, sp))
    #Y_mis_var[[s]] <- array(0, dim=c(n_t[s], p, sp))
    if(outputlevel == 1) {
      Lambdaout[[s]] <- array(0, dim=c(p, j_t[s], sp))
      psiout[[s]] <- array(0, dim=c(p, 1, sp))
      #l_tout[[s]] <- zeros(n_t[s], j_t[s])
    }
    if(outputlevel == 3) {
      #Lambdaout[[s]] <- zeros(p, j_t[s])
      #psiout[[s]] <- zeros(p, 1)
      l_tout[[s]] <- array(0, dim=c(n_t[s], j_t[s], sp))
      #SigmaLambda[[s]] <- zeros(p, p)
      Beta_tout[[s]] <- array(0, dim=c(p, d_cov, sp))
      Y_obs_all[[s]] <- array(0, dim=c(n_t[s], p, sp))
      Y_mis_all[[s]] <- array(0, dim=c(n_t[s], p, sp))
    }
    if(outputlevel == 2) {
      SigmaLambda[[s]] <- zeros(p, p)
      Lambdaout[[s]] <- zeros(p, j_t[s])
      #psiout[[s]] <- zeros(p, 1)
      l_tout[[s]] <- zeros(n_t[s], j_t[s])
      l_tout_mis[[s]] <- zeros(n_t[s], j_t[3-s])
      Beta_tout[[s]] <- zeros(p, d_cov)
      Lambda_l_out[[s]] <- zeros(p, n_t[s])
    }
    if(outputlevel == 4) {
      #psiout[[s]] <- zeros(p, 1)
      #l_tout[[s]] <- zeros(n_t[s], j_t[s])
    }
  }
  CE_mean <- zeros(p, sp)
  CE_ite <- zeros(n_t[1]+n_t[2], p)
  
  ##################################
  # --- start posterior sampling ---
  ##################################
  
  for(r in 1:nrun)
  { 
    #### Step 1: treatment-specific latent factors
    l_t2 <- list()
    for (s in 1:2){
      # Y - beta X
      Y_centr[[s]] <- Y_t[[s]]-t(Beta_t[[s]]%*%X_t[[s]])
      
      Lmsg1 <- vecmat(psi_t[[s]], Lambda_t[[s]])
      Veta1 <- diag(j_t[s]) + t(Lmsg1) %*% Lambda_t[[s]]
      T1 <- chol(Veta1)
      qrT1 <- qr(T1)
      R1 <- qr.R(qrT1)
      S1 <- solve(R1)
      Veta11 <- tcrossprod(S1)
      Meta1 <- Y_centr[[s]] %*% Lmsg1 %*% Veta11
      x1 <- matrix(rnorm(n_t[s] * j_t[s]), nrow = n_t[s], ncol = j_t[s])
      l_t[[s]] <- Meta1 + x1 %*% t(S1)
      l_t2[[s]] <- crossprod(l_t[[s]])
    }
    
    ### Step 2: specific factor loadings, with constraints
    for(s in 1:2){
      for(i in 1:p){
        Qlam <- diag(Plam_t[[s]][i,], j_t[s]) + psi_t[[s]][i] * l_t2[[s]]
        blam <- psi_t[[s]][i] * (t(l_t[[s]]) %*% Y_centr[[s]][, i])
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
      Ytil_t <- Y_centr[[s]] - (l_t[[s]] %*% t(Lambda_t[[s]]))
      psi_t[[s]] <- rgamma(p, shape = a_psi_t[s] + 0.5 * n_t[s], rate = b_psi_t[s] + 0.5 * colSums(Ytil_t^2))
      Psi_t[[s]] <- diag(1 / psi_t[[s]])
    }
    
    ### Step 5: beta regression
    for(s in 1:2){
      Y_star <- Y_t[[s]] - (l_t[[s]])%*%t(Lambda_t[[s]])
      prod_X <- (X_t[[s]])%*%t(X_t[[s]])
      
      V_beta <- lapply(1:p, function(i) chol2inv(chol(prod_X*(psi_t[[s]][i]) + (1/bbeta_t[s])*diag(d_cov))))
      M_beta <- abeta_t[s]/bbeta_t[s] + (X_t[[s]])%*%Y_star*(matrix(psi_t[[s]],ncol=p, nrow=d_cov, byrow=TRUE))
      
      xP1 <- matrix(rnorm(d_cov*p), nrow = d_cov, ncol = p)
      Beta_t[[s]] <- t(sapply(1:p, function(i) M_beta[,i]%*%V_beta[[i]] + xP1[,i] %*% chol(V_beta[[i]])))
    }
    
    ### Step IMPUTATION missing data
    for (s in 1:2){
      # I do not observe Y, so the treat-spec factor score just come from the prior
      l_mis[[s]] <- matrix(rnorm(n_t[s] * j_t[3-s]), nrow = n_t[s], ncol = j_t[3-s])
    }
    
    if(r > burn){
      neff <- (r - burn) / thin
      teff <- (nrun - burn) / thin
      for (s in 1:2){
        X_beta_obs <- t(Beta_t[[s]]%*%X_t[[s]])
        X_beta_mis <- t(Beta_t[[3-s]]%*%X_t[[s]])
        
        Y_obs_lt[[s]] <- X_beta_obs + l_t[[s]] %*% t(Lambda_t[[s]]) + 
          matrix(rnorm(n_t[s]*p),n_t[s],p)*matrix(1/sqrt(psi_t[[s]]),n_t[s],p, byrow = TRUE)
        Y_mis_lt[[s]] <- X_beta_mis + l_mis[[s]] %*% t(Lambda_t[[3-s]]) +
          matrix(rnorm(n_t[s]*p),n_t[s],p)*matrix(1/sqrt(psi_t[[3-s]]),n_t[s],p, byrow = TRUE)
        
        #Y_obs_var[[s]] <- X_beta_obs + mvrnorm(n_t[s], rep(0, p), 
        #                          tcrossprod(Lambda_t[[s]])+Psi_t[[s]])
        #Y_mis_var[[s]] <- X_beta_mis + mvrnorm(n_t[s], rep(0, p), 
        #                          tcrossprod(Lambda_t[[3-s]])+Psi_t[[3-s]])
      }
      if(outputlevel==1)
      {
        for(s in 1:2){
          Lambdaout[[s]][, , neff] <- Lambda_t[[s]]
          psiout[[s]][, , neff] <- 1 / psi_t[[s]]
          #l_tout[[s]] <- l_tout[[s]] + l_t[[s]] / teff
        }
      }
      if(outputlevel==2){
        for(s in 1:2){
          #Lambdaout[[s]] <- Lambdaout[[s]] + Lambda_t[[s]]  / teff
          #psiout[[s]] <- psiout[[s]] + (1 / psi_t[[s]]) / teff
          #l_tout[[s]] <- l_tout[[s]] + l_t[[s]] / teff
          #l_tout_mis[[s]] <- l_tout_mis[[s]] + l_mis[[s]] / teff
          SigmaLambda[[s]] <- SigmaLambda[[s]] + tcrossprod(Lambda_t[[s]]) / teff
          Beta_tout[[s]] <- Beta_tout[[s]] + Beta_t[[s]] / teff
          Lambda_l_out[[s]] <- Lambda_l_out[[s]] + (Lambda_t[[s]] %*% t(l_t[[s]])) / teff
        }
        CE_all <- rbind(Y_mis_lt[[1]]-Y_obs_lt[[1]],Y_obs_lt[[2]]-Y_mis_lt[[2]])
        CE_mean[ , neff] <- apply(CE_all,2,mean) 
        CE_ite <- CE_ite + CE_all/teff
      }
      if(outputlevel==3){
        for(s in 1:2){
          #Lambdaout[[s]] <- Lambdaout[[s]] + Lambda_t[[s]]  / teff
          #psiout[[s]] <- psiout[[s]] + (1 / psi_t[[s]]) / teff
          #l_tout[[s]][, , neff] <- l_t[[s]]
          #SigmaLambda[[s]] <- SigmaLambda[[s]] + tcrossprod(Lambda_t[[s]]) / teff
          Beta_tout[[s]][ , , neff] <- Beta_t[[s]] 
          Y_obs_all[[s]][ , , neff] <- Y_obs_lt[[s]]
          Y_mis_all[[s]][ , , neff] <- Y_mis_lt[[s]]
        }
      }
    }
    if (trace & r %% nprint == 0) cat("r=",r,"/",nrun,"\n")
  }
  ##### Save and exit
  if(outputlevel == 1)  {
    l_tout <- NULL
    CE_mean <- CE_ite <- NULL
  }
  if(outputlevel == 3)  {
    Lambdaout <- l_tout <- l_mis <- NULL
    psiout <- NULL
    CE_mean <- CE_ite <- NULL
  }
  if(outputlevel == 2)  {
    Lambdaout <- l_tout <- l_tout_mis <- l_mis <- NULL
    psiout <- NULL
    Y_obs_all <- Y_mis_all <- NULL
  }
  if(outputlevel == 4)  {
    Lambdaout <- NULL
    psiout <- NULL
    l_tout <- NULL
    SigmaLambda <- NULL
    CE_mean <- CE_ite <- NULL
  }
  out <- list(Lambda = Lambdaout, psi = psiout, l_t = l_tout, l_mis = l_tout_mis,
              SigmaLambda = SigmaLambda,
              CE_mean = CE_mean, CE_ite = CE_ite,
              Y_obs_lt = Y_obs_lt, Y_mis_lt = Y_mis_lt,
              Y_obs_all = Y_obs_all, Y_mis_all = Y_mis_all,
              #Y_obs_var = Y_obs_var, Y_mis_var = Y_mis_var,
              Beta_tout = Beta_tout, Lambda_l_out = Lambda_l_out)
  return(structure(out,  class="causal_fa"))
}

