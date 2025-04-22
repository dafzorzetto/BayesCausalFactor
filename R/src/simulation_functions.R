# common factors simulation
funct_common_factors<-function(p,k, no_zeros=0.6, perc_nozeros){
  set.seed(1)
  
  PH <- as.vector(zeros(p, k))
  noZEROc <- (p*k)*perc_nozeros                               # number of no zero entries
  studyc <- runif(noZEROc, no_zeros,1)                        # no zero values (values between 0.6 and 1)
  sign <- sample(x = length(studyc), size = (length(studyc)/2))
  studyc[sign] <- studyc[sign]*(-1)
  positionc <- sample(x = k*p, size = length(studyc))   # random position for the no zero entries
  PH[positionc] <- studyc
  Phi <- matrix(PH,p,k)
  
  return(Phi)
}

# treatment-specific factors simulation
funct_specific_factors<-function(p, j_t, no_zeros=0.6, perc_nozeros){
  set.seed(1)
  
  L <- noZERO <- studyc <- position <- Lambda_s <- list()
  
  for(s in 1:2){
    L[[s]] <- as.vector(zeros(p, j_t[s]))
    noZERO[[s]] <- (p*j_t[s])*perc_nozeros
    studyc[[s]] <- runif(noZERO[[s]], no_zeros,1)
    sign <- sample(x = length(studyc), size = (length(studyc)/2))
    studyc[[s]][sign] <- studyc[[s]][sign]*(-1)
    position[[s]] <- sample(x = p*j_t[s], size=length(studyc[[s]]))
    L[[s]][position[[s]]] <- studyc[[s]]
    Lambda_s[[s]] <- matrix(L[[s]], p, j_t[s])
  }
  
  return(Lambda_s)
}

# regretion for potential multivariate outocome Y(0), Y(1)
funct_regres_Y<-function(p, k, j_t, 
                         lim_nozeros_Phi=0.6, lim_nozeros_Lambda=0.6, 
                         range_beta=c(-1,1),
                         dim_x_i,n_s, X_t,
                         param_x=c(0,1),
                         perc_nozeros_phi=0.5, perc_nozeros_lambda=0.5,
                         clusters_allocation, cl_miss,
                         mean_cl=list(0,0), var_cl=list(1,1)){
  
  Psi_s <- Psi_mis <- beta_t <- Sigma_s <- Y_t <- f_t <- l_t <- list()
  Y_mis <- l_mis <- list()
  
  Phi <- funct_common_factors(p, k, no_zeros=lim_nozeros_Phi, perc_nozeros=perc_nozeros_phi)
  Lambda_s <- funct_specific_factors(p, j_t, no_zeros=lim_nozeros_Lambda, perc_nozeros=perc_nozeros_lambda)
  
  for(s in 1:2){
    # errors
    Psi_s[[s]] = diag(runif(p,0,1), p)    
    Psi_mis[[s]] = diag(runif(p,0,1), p) 
    
    #covariates 
    betas <- seq(range_beta[1]+(s-1)*0.5,range_beta[2]+(s-1)*0.2,
                 length.out=dim_x_i*p)
    beta_t[[s]] = matrix(betas, nrow=p, ncol=dim_x_i, byrow=TRUE)
    
    # common factors
    f_t[[s]] <- matrix(rnorm(n_s[s]*k), n_s[s], k)
    
    # treatment-specific factors
    l_t[[s]] <- sapply(1:j_t[s], function(j) rnorm(n_s[s], 
                                                   mean_cl[[s]][[j]][clusters_allocation[[s]][,j]], 
                                                   var_cl[[s]][[j]][clusters_allocation[[s]][,j]]))
    
    # regression 
    Y_t[[s]] <- t(X_t[[s]])%*%t(beta_t[[s]]) + f_t[[s]]%*%t(Phi) + l_t[[s]]%*%t(Lambda_s[[s]]) +
      mvrnorm(n_s[s], rep(0, p), Psi_s[[s]])
  }
  
  ## "missing" data
  for (s in 1:2){
    l_mis[[s]] <- sapply(1:j_t[3-s], function(j) rnorm(n_s[s],
                                                       mean_cl[[3-s]][[j]][cl_miss[[s]][,j]],
                                                       var_cl[[3-s]][[j]][cl_miss[[s]][,j]]))
    Y_mis[[s]] <- t((beta_t[[3-s]])%*%(X_t[[s]])) + f_t[[s]]%*%t(Phi) +
      l_mis[[s]]%*%t(Lambda_s[[3-s]]) +
      mvrnorm(n_s[s], rep(0, p), Psi_mis[[3-s]])
  }
  
  return(list(Phi = Phi, Lambda_s = Lambda_s, Psi_s = Psi_s, 
              beta_t = beta_t, Y_t = Y_t, Y_mis = Y_mis, f_t = f_t,
              l_t = l_t, l_mis = l_mis))
}

# simulation of all the variables: cofounders, treatment, potential outcome
funct_simulation<-function(dim_x_i=2, n_tot=60, 
                           p=6, k=3, j_t=rep(1,2),
                           lim_nozeros_Phi=0.6, lim_nozeros_Lambda=0.6,
                           perc_nozeros_phi=0.5, perc_nozeros_lambda=0.5,
                           mean_cl=list(0,0), var_cl=list(1,1),
                           seed){
  
  set.seed(seed)
  range_beta=c(-0.8,0.8)
  
  # confonders
  # number confounders: dim_x_i
  X <- matrix(rnorm(dim_x_i*n_tot, 0,1), nrow=dim_x_i, ncol=n_tot)
  
  # treatment
  # expit of a covariates function
  random_coef <- runif(dim_x_i,0.05,0.2)  
  reg_T <- random_coef%*%X
  logit_T <- exp(reg_T)/(1+exp(reg_T))
  Treat <- rbinom(n_tot,1,logit_T)
  
  # treatment-specific information
  n_s <- c(table(Treat))
  X_t  <- list()
  for (s in 0:1){
    X_t[[s+1]] <- X[ ,Treat==s]
  }
  
  # cluster allocation
  clusters_allocation <- cl_miss <- list(matrix(NA,n_s[1],j_t[1]),
                                         matrix(NA,n_s[2],j_t[2]))
  
  for (s in 0:1){
    for (j in 1:j_t[s+1]){
      n_cl=length(mean_cl[[s+1]][[j]])
      
      if(n_cl==1){
        cl_tr <- rep(1,n_tot)
        clusters_allocation[[s+1]][,j] <- cl_tr[Treat==s]
        cl_miss[[2-s]][,j] <- cl_tr[Treat==(1-s)]
      }
      if(n_cl==2){
        cl_tr <- rep(1,n_tot)
        cl_tr[(X[1,]+X[2,])>0] <- 2
        clusters_allocation[[s+1]][,j] <- cl_tr[Treat==s]
        cl_miss[[2-s]][,j] <- cl_tr[Treat==(1-s)]
      }
      if(n_cl==3){
        cl_tr <- rep(1,n_tot)
        cl_tr[X[1,]*2+(X[2,]-0.5)*1.5>0] <-2
        cl_tr[(X[1,]+0.5)*(1.5)+(X[3,]+0.3)*2<0] <-3
        clusters_allocation[[s+1]][,j] <- cl_tr[Treat==s]
        cl_miss[[2-s]][,j] <- cl_tr[Treat==(1-s)]
      }
    }
  }
  
  # outcome and factors
  Y_factors <- funct_regres_Y(p=p, k=k, j_t=j_t, 
                              lim_nozeros_Phi=lim_nozeros_Phi, 
                              lim_nozeros_Lambda=lim_nozeros_Lambda,
                              range_beta=range_beta,
                              dim_x_i=dim_x_i,n_s=n_s,
                              X_t=X_t,
                              perc_nozeros_phi=perc_nozeros_phi,
                              perc_nozeros_lambda=perc_nozeros_lambda,
                              clusters_allocation=clusters_allocation, cl_miss=cl_miss,
                              mean_cl=mean_cl, var_cl=var_cl)
  
  return(list(data=list( Y_t = Y_factors$Y_t, Treat = Treat, n_s = n_s, X_t = X_t),
              factors=list(Phi = Y_factors$Phi, Lambda_s = Y_factors$Lambda_s, 
                           Psi_s = Y_factors$Psi_s, f_t = Y_factors$f_t,
                           l_t = Y_factors$l_t),
              parameters=list(k = k, j_t = j_t, beta_t = Y_factors$beta_t,
                              clusters_allocation = clusters_allocation, 
                              mean_cl = mean_cl, var_cl = var_cl),
              mis_data=list(Y_mis = Y_factors$Y_mis,  l_mis = Y_factors$l_mis, cl_miss = cl_miss)))
  
}