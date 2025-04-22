# Phi*Phi^T
phi_post <- function(simulation,post_models,n_set,
                     breaks=c(-2,-0.8,-0.3,0.3,0.8,2)){
  
  SigmaPhi_all <- array(0, dim=c(rep(dim(simulation$Phi)[1],2), n_sample))
  for (i in 1:n_sample){
    SigmaPhi_all[,,i] <- (post_models[[i]]$SigmaPhi)
  }
  SigmaPhi_median <- apply(SigmaPhi_all,c(1,2),median)
  
  par(mfrow=c(1,2), mar=c(5.1, 4.1,4.1,4.1))
  plot(tcrossprod(simulation$Phi), 
       main=paste0("Scenario ",n_set), breaks=breaks,
       xlab=" ", ylab=expression(paste(Phi, Phi^T,"- simulated")), 
       border=NA, digits=2, 
       text.cell=list(cex=0.45),
       key=list(side=4,  font=1, cex.axis=0.75))
  plot(SigmaPhi_median, 
       main=paste0("Scenario ",n_set), breaks=breaks,
       xlab=" ", ylab=expression(paste(Phi, Phi^T," - estimated")), 
       border=NA, digits=2, 
       text.cell=list(cex=0.45),
       key=list(side=4,  font=1, cex.axis=0.75))
}

# Lambda_t*Lambda_t^T
lambda_post <- function(simulation,post_models, n_set,
                        breaks=c(-2,-0.8,-0.3,0.3,0.8,2)){
  
  for(s in 1:2){
    
    if (s==1) {
      a_sim <- expression(paste(Lambda[0], Lambda[0]^T,"- simulated"))
      a_est <- expression(paste(Lambda[0], Lambda[0]^T,"- estimated"))
    } else {
      a_sim <- expression(paste(Lambda[1], Lambda[1]^T,"- simulated"))
      a_est <- expression(paste(Lambda[1], Lambda[1]^T,"- estimated"))
    }
    
    SigmaLambda_all <- array(0, dim=c(rep(dim(simulation$factors$Lambda_s[[s]])[1],2), n_sample))
    for (i in 1:n_sample){
      SigmaLambda_all[,,i] <- (post_models[[i]]$SigmaLambda[[s]])
    }
    SigmaLambda_median <- apply(SigmaLambda_all,c(1,2),median)
    
    par(mfrow=c(1,2), mar=c(5.1, 4.1,4.1,4.1))
    plot(tcrossprod(simulation$factors$Lambda_s[[s]]), 
         main=paste0("Scenario ",n_set), breaks=breaks,
         xlab=" ", ylab=a_sim, border=NA, digits=2, 
         text.cell=list(cex=0.45), 
         key=list(side=4,  font=2, cex.axis=0.75))
    plot(SigmaLambda_median, 
         main=paste0("Scenario ",n_set), breaks=breaks,
         xlab=" ", ylab=a_est, border=NA, digits=2, 
         text.cell=list(cex=0.45),
         key=list(side=4,  font=2, cex.axis=0.75))
  }
}

# causal effects
CE_comparison<-function(simulation,post_models, n_set,
                        breaks=c(-8,-0.6,-0.3,0.3,0.6,8)){
  
  CE_sim <- sapply(1:n_sample, function(i)
    apply(rbind(simulation[[i]]$Y_mis[[1]]-simulation[[i]]$Y_t[[1]],
                -simulation[[i]]$Y_t[[2]]+simulation[[i]]$Y_mis[[2]]),
          2,mean))
  CE_est_y <- sapply(1:n_sample, function(i)
    apply(rbind(post_models[[i]]$Y_mis_distr[[1]]-post_models[[i]]$Y_obs_distr[[1]],
                -post_models[[i]]$Y_obs_distr[[2]]+post_models[[i]]$Y_mis_distr[[2]]),
          2,mean))
  CE_est_direct <- sapply(1:n_sample, function(i)
    apply(rbind(post_models[[i]]$CE_ite[[1]],-post_models[[i]]$CE_ite[[2]]),
          2,mean))
  
  par(mfrow=c(1,3), mar=c(5.1, 4.1,4.1,4.1))
  plot(matrix(apply(CE_sim,1,mean),ncol=1), 
       main=paste0("Scenario ",n_set), breaks=breaks,
       xlab=" ", ylab="causal effects - simulated", 
       border=NA, digits=2, 
       text.cell=list(cex=0.75), 
       key=list(side=4,  font=1, cex.axis=0.75))
  plot(matrix(apply(CE_est_y,1,mean),ncol=1),  
       main=paste0("Scenario ",n_set), breaks=breaks,
       xlab=" ", ylab="causal effects - estimated with Y", 
       border=NA, digits=2, 
       text.cell=list(cex=0.75),
       key=list(side=4,  font=1, cex.axis=0.75))
  plot(matrix(apply(CE_est_direct,1,mean),ncol=1),  
       main=paste0("Scenario ",n_set), breaks=breaks,
       xlab=" ", ylab="causal effects - estimated", 
       border=NA, digits=2, 
       text.cell=list(cex=0.75),
       key=list(side=4,  font=1, cex.axis=0.75))
}

results_boxplot <- function(settings,n_settings, name_title, name_y, ylim_par){
  boxplot_dataset <- cbind(ind=settings, sim=c(rep(1:n_settings, each = n_sample)))
  colour_palette=c("#FAD401","#E09600","#E06100","#E02700","#FA31E5","#B600FA","#5D00F5","#4100EB")
  
  boxplot(ind ~ sim, data=boxplot_dataset, 
          main=name_title, ylab=name_y, xlab="Scenarios",
          col=colour_palette[1:n_settings],
          ylim=ylim_par)
  abline(h=0, col="green")
}

each_CE_plot<-function(bias_setting, dim_Y, name_title,name_y,ylim_par){
  
  n_units=dim(bias_setting)[1]/dim_Y
  indicator_Y<- model.matrix(~ factor(rep(1:dim_Y, each=n_units)) - 1) 
  bias<-sapply(1:(dim(bias_setting)[2]), function(r) (t(bias_setting[,r])%*%indicator_Y)/n_units)
    
  boxplot_dataset <- cbind(bias=c(bias), 
                           sim=rep(1:dim_Y, n_sample))
  
  boxplot(bias ~ sim, data=boxplot_dataset, 
          main=name_title, ylab=name_y, xlab="outcomes",
          col=rainbow(17)[16:(17-dim_Y)],
          ylim=ylim_par)
  abline(h=0, col="#FAD401")
}

each_CE_plot_mis<-function(bias_setting, dim_Y, name_title,name_y,ylim_par){
  
  n_units <- sapply(bias_setting, function(i) length(i)/dim_Y)
  indicator_Y <- sapply(n_units, function(i) model.matrix(~ factor(rep(1:dim_Y, each=i)) - 1))
  bias <- sapply(1:length(bias_setting), function(i) (t(bias_setting[[i]])%*%indicator_Y[[i]])/n_units[i])
  
  boxplot_dataset <- cbind(bias=c(bias), 
                           sim=rep(1:dim_Y, n_sample))
  
  boxplot(bias ~ sim, data=boxplot_dataset, 
          main=name_title, ylab=name_y, xlab="outcomes",
          col=rainbow(17)[16:(17-dim_Y)],
          ylim=ylim_par)
  abline(h=0, col="#FAD401")
}
