#' @title
#' Estimation of BART and BCF
#'
#' @description
#' 2 function for the estimation of BART and BCF, respectively
#'
#' @param Y_t : outsome (multivariate) divided by treatment levels (control and treatment)
#' @param X_t : covariates divided by treatment levels (control and treatment)
#' @param n_t : number of units in each treatment level (control and treatment)
#'
#' @return
#' For BART and BCF functions:
#' - *tau* : vector of individual treatment effect (ITE)
#' - *Y_imp* : list of the imputed potential outcome
#'
#' @import bartCause
#' @import bcf

#########################################################################
# libraries
library(bartCause)
library(bcf)

#########################################################################
#    ---     BART  ----
#########################################################################

# function to estimate the BART and to extrapolate the quantities of interest
bart_estimation<-function(Y_t, X_t, n_t){
  
  # set seed for riproducibility
  set.seed(1)
  
  # ------   prearing variables   ------
  
  # main variables
  T_level=c(rep(0, n_t[1]),rep(1, n_t[2]))    # treatment
  X=rbind(t(X_t[[1]]), t(X_t[[2]]) )          # covariates-confounders
  Y_obs=rbind(Y_t[[1]], Y_t[[2]])             # observed outcome
  dim_Y = dim(Y_obs)[2]
  
  # ------   BART estimation   ------
  extract_bart_results<-function(each_Y){
    bart_fit <- bartCause::bartc(as.matrix(each_Y),as.matrix(T_level), as.matrix(X),
                                 n.samples = 1000, n.burn = 1000)

     # espcted values for Y -> Y imputed
     #Y_imp=list(Y_0=apply(bartCause::extract(bart_fit, type = "y.0"),2,mean),
    #            Y_1=apply(bartCause::extract(bart_fit, type = "y.1"),2,mean))
     
    Y_obs_chain = bart_fit$mu.hat.obs[1,,]
    Y_cf_chain = bart_fit$mu.hat.cf[1,,]
    
    tau <- cbind(Y_cf_chain[,T_level==0] - Y_obs_chain[,T_level==0],
                 Y_obs_chain[,T_level==1] - Y_cf_chain[,T_level==1])
    
    tau_Y <- colMeans(tau)
    tau_CE <- rowMeans(tau)
    
    #tau <- apply(tau,1,mean)
    #quantiles_CE = quantile(tau_CE, prob=c(0.025,0.05,0.5,0.95,0.975))
     
     return(list(tau_Y=tau_Y,
                 mean_CE = tau_CE))

  }
  # estimation
  bart_fit <- lapply(1:dim_Y, function(y) extract_bart_results(Y_obs[,y]))
  
  return(bart_fit)
}

#########################################################################
#    ---    BCF    ----
#########################################################################

# function to estimate the BCF and to extrapolate the quantities of interest
BCF_estimation<-function(Y_t, X_t, Treat){
  
  # set seed for riproducibility
  set.seed(1)
  
  # ------   prearing variables   ------
  
  # main variables
  n_t = c(dim(X_t[[1]])[2],dim(X_t[[2]])[2])
  T_level=c(rep(0, n_t[1]),rep(1, n_t[2]))    # treatment
  X=rbind(t(X_t[[1]]), t(X_t[[2]]) )          # covariates-confounders
  Y_obs=rbind(Y_t[[1]], Y_t[[2]])             # observed outcome
  dim_Y = dim(Y_obs)[2]
  
  # ------   BCF estimation   ------
  
  extract_bcf_results<-function(each_Y){
    #propensity score
    p.score <- glm(T_level ~ X,
                   family = binomial,
                   data = as.data.frame(cbind(T_level, X)))
    pihat <- predict(p.score, as.data.frame(X))
    
    # estimation
    bcf_fit <- bcf(each_Y, T_level, X, X, pihat,
                   nburn = 1000, nsim = 250)
    
    #temp_tau=colMeans(bcf_fit$tau)
    #tau_def=rep(NA,sum(n_t))
    #tau_def[Treat==0]=temp_tau[1:n_t[1]]
    #tau_def[Treat==1]=temp_tau[(n_t[1]+1):sum(n_t)]
    
    return(tau=colMeans(bcf_fit$tau))
  }
  
  # estimation
  bcf_fit_all <- sapply(1:dim_Y, function(y) extract_bcf_results(each_Y = Y_obs[,y]))
  
  print("one more sample done")
  return(list(bcf_fit_all))
}
