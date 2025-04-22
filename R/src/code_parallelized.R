
##### BayesCausalFactor model
BayesCausalFactor_parall<-function(data, factors){
  return(BayesCausalFactor(Y_t = data$data$Y_t,
                                      j_t = factors,
                                      X_weight = list(rbind(1,data$data$X_t[[1]]),
                                                      rbind(1,data$data$X_t[[2]])),
                                      X_reg = list(rbind(1,data$data$X_t[[1]]),
                                                   rbind(1,data$data$X_t[[2]])),                                  
                                      R_cluster=6,
                                      outputlevel = 3,
                                      nrun = 6000, burn = 4000,
                                      nprint = 1000))
}

##### StandardFactor model
# ex: causal_fa_Lt_X
StandardFactor_parall<-function(data){
  return(StandardFactor(Y_t = data$data$Y_t,
                        X_t = data$data$X_t,
                        j_t = data$parameters$j_t,
                        outputlevel = 2,
                        nrun = 5000, burn = 3500,
                        nprint = 1000))
}


##### BART model
parallel_funct_bart<-function(data){
  return(bart_estimation(Y_t = data$data$Y_t, 
                         X_t = data$data$X_t, 
                         n_t = c(sum(data$data$Treat==0),sum(data$data$Treat==1))))
}
