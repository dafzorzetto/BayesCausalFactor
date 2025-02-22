###### ---- PLOTS  ---- 

#########################################################

ce_ourmodel=model_00_paper$CE_across
ce_standFA=stand_FM_00$CE_across
bart_results=bart_00
bcf_results=bcf_00[[1]]
lims = c(0,50)
month='Y(1)-Y(0)'
gray_line=0



mean_our=apply(ce_ourmodel[, 1001:1500],1,median)
CI_our=apply(ce_ourmodel[, 1001:1500],1, quantile, probs=c(0.05,0.95))
mean_sFA=apply(ce_standFA[, 1001:1500],1,median)
CI_sFA=apply(ce_standFA[, 1001:1500],1, quantile, probs=c(0.05,0.95))
CI_bart=apply(bart_results,2, quantile, probs=c(0.05,0.95))
CI_bcf=apply(bcf_results,2, quantile, probs=c(0.05,0.95))
colors <- colorRampPalette(brewer.pal(11, "Spectral"))(26)

colors_paper <- c(rep("#694573",3), rep("#8E414E",3),
                  rep("#F99D5B",8), rep("#FDCE57",2),
                  rep("#BDDEF2",2), rep("#86BF6B",5),
                  rep("#A86A8C",2), rep("#72A3A8",2))

plot(mean_our,27:1, pch=16, col=colors_paper,
     xlim=c(min(min(CI_our),0)-0.2, min(max(CI_our),120)+0.45),
     ylab='', xlab=' ', #main=month,
     yaxt = "n",bty = "n", cex=1)
segments(CI_our[1,], 27:1, CI_our[2,],27:1, col=colors_paper, lwd = 1)
points(apply(bart_results, 2, mean),27:1+0.2, 
       col=colors_paper, cex=0.7, pch=8)
points(apply(bcf_results, 2, mean),27:1+0.4, 
       col=colors_paper, cex=0.7, pch=17)
segments(CI_bart[1,], 27:1+0.2, CI_bart[2,],27:1+0.2, 
         col=colors_paper, lty=3)
segments(CI_bcf[1,], 27:1+0.4, CI_bcf[2,],27:1+0.4, 
         col=colors_paper, lty=2)
text(par("usr")[1]+0.1,(27:1)+0.2, cex=0.9, col=colors_paper,
     names_chemicals, srt = 0, xpd = TRUE, adj = 1)
abline(v=0, col="gray")

legend(0.6,18, c("Our model","BART","BCF"),
       pch=c(19,8,17), 
       col=1, 
       lty=c(1,3,2), cex=0.9)

plot(mean_our,27:1, pch=16, col=1,
     xlim=c(min(min(CI_our),0)-0.2, min(max(CI_our),120)+0.5),
     ylab='', xlab=' ', #main=month,
     yaxt = "n",bty = "n", cex=1)
segments(CI_our[1,], 27:1, CI_our[2,],27:1, col=1, lwd = 1)
points(apply(bart_results, 2, mean),27:1+0.2, 
       col=1, cex=0.7, pch=8)
points(apply(bcf_results, 2, mean),27:1+0.4, 
       col=1, cex=0.7, pch=17)
segments(CI_bart[1,], 27:1+0.2, CI_bart[2,],27:1+0.2, 
         col=1, lty=3)
segments(CI_bcf[1,], 27:1+0.4, CI_bcf[2,],27:1+0.4, 
         col=1, lty=2)
text(par("usr")[1]+0.1,(27:1)+0.2, cex=0.9, col=colors_paper,
     names_chemicals, srt = 0, xpd = TRUE, adj = 1)
abline(v=0, col="gray")

legend(0.6,21.5, c("Our model","BART","BCF"),
       pch=c(19,8,17), 
       col=c(1,1,1),
       lty=c(1,3,2), cex=0.85)
legend(0.56,18, c("Alkaline-earth metals", "Alkalini metals",
                 "Transition metals", "metalloids","Other metals",
                 "Nonmetals","Halogens","Organics"),
       #col=c(1,1,1,rep(NA,9)), 
       text.col=c(unique(colors_paper)), 
       cex=0.83)


plot(mean_our,27:1, pch=16, col=2,
     xlim=c(min(min(CI_our),0)-0.2, min(max(CI_our),120)+0.45),
     ylab='', xlab=' ', #main=month,
     yaxt = "n",bty = "n", cex=1)
segments(CI_our[1,], 27:1, CI_our[2,],27:1, col=2, lwd = 1)
points(apply(bart_results, 2, mean),27:1+0.2, 
       col=3, cex=0.7, pch=8)
points(apply(bcf_results, 2, mean),27:1+0.4, 
       col=4, cex=0.7, pch=17)
segments(CI_bart[1,], 27:1+0.2, CI_bart[2,],27:1+0.2, 
         col=3, lty=6)
segments(CI_bcf[1,], 27:1+0.4, CI_bcf[2,],27:1+0.4, 
         col=4, lty=2)
text(par("usr")[1]+0.1,(27:1)+0.2, cex=0.9, col='black',
     names_chemicals, srt = 0, xpd = TRUE, adj = 1)
abline(v=0, col="gray")

legend(0.6,18, c("Our model","BART","BCF"),
       pch=c(19,8,17), 
       col=c(2,3,4), 
       lty=c(1,6,2), cex=0.85, bty = "n")


palette_simulations <- c('#673F73',"#30688E","#3CB97D","#FDD225")
plot(mean_our,27:1, pch=16, col=palette_simulations[1],
     xlim=c(min(min(CI_our),0)-0.2, min(max(CI_our),120)+0.45),
     ylab='', xlab=' ', #main=month,
     yaxt = "n",bty = "n", cex=1)
segments(CI_our[1,], 27:1, CI_our[2,],27:1, col=palette_simulations[1], lwd = 1)
points(mean_sFA,27:1+0.2, 
       col=palette_simulations[2], cex=0.7, pch=8)
segments(CI_sFA[1,], 27:1+0.2, CI_sFA[2,],27:1+0.2, 
         col=palette_simulations[2], lty=1)
points(apply(bart_results, 2, mean),27:1+0.4, 
       col=palette_simulations[3], cex=0.7, pch=17)
segments(CI_bart[1,], 27:1+0.4, CI_bart[2,],27:1+0.4, 
         col=palette_simulations[3], lty=1)
points(apply(bcf_results, 2, mean),27:1+0.6, 
       col=palette_simulations[4], cex=0.7, pch=4)
segments(CI_bcf[1,], 27:1+0.6, CI_bcf[2,],27:1+0.6, 
         col=palette_simulations[4], lty=1)
text(par("usr")[1]+0.1,(27:1)+0.2, cex=0.9, col='black',
     names_chemicals, srt = 0, xpd = TRUE, adj = 1)
abline(v=0, col="gray")

legend(0.6,18, c("Our model","standardFA","BART","BCF"),
       pch=c(19,8,17,4), 
       col=palette_simulations, 
       lty=1, cex=0.85, bty = "n")

#########################################################

#########################################################

selection_factors_varimax <- function(sigma_matrix,treat_level, th_var){
  
  eigenPhi <- eigen(sigma_matrix)
  d_value <- eigenPhi$values 
  choiceK <- d_value / sum(d_value)
  
  print(round(choiceK,3))
  print(round(cumsum(choiceK),3))
  
  k <- length(which(choiceK>th_var))
  loadK <- varimax(eigenPhi$vectors[,1:k])$loadings
  prova<- loadK[1:27,]
  rownames(prova) <- rep(" ",27)
  def_matrix <- sapply(1:k, function(i) sign(prova[which.max(abs(prova[,i])),i])*prova[,i])
  
  max_loadings <- sapply(1:k, function(i) which.max(def_matrix[,i])) 
  def_matrix <- def_matrix[,order(max_loadings)]
  colnames(def_matrix) <- rep(" ",k)
  
  par(mfrow=c(1,1),mar=c(3,6.4,4,4.1)) #mar=c(3, 6,4,4.5))
  breaks=c(-0.9,-0.6,-0.45,-0.25,-0.15,0.15,0.25,0.45,0.6,0.9)
  #breaks=c(-0.9,-0.6,-0.45,-0.2,0.2,0.45,0.6,0.9)
  cols <- brewer.pal(length(breaks)-1, "RdBu")
  colors_paper <- c(rep("#694573",3), rep("#8E414E",3),
                    rep("#F99D5B",8), rep("#FDCE57",2),
                    rep("#BDDEF2",2), rep("#86BF6B",5),
                    rep("#A86A8C",2), rep("#72A3A8",2))
  plot(def_matrix, 
       main=treat_level, breaks=breaks,
       xlab=" ", ylab=" ", yaxt = "n",xaxt = "n", bty = "n",
       border=NA, digits=2, 
       text.cell=list(cex=0.45),
       col=cols,
       key=list(side=4,  font=1, cex.axis=0.75))
  for (i in 1:27){
    axis(2, at = i, labels = names_chemicals[28-i], 
         col.axis=colors_paper[28-i], las = 2, cex.axis = 0.75)
  }
  axis(1, at = 1:k, labels = c(1:k), las = 1, cex.axis = 0.75)
}


selection_factors_varimax(sigma_matrix = model_00_paper$sigma[[1]], 
                          treat_level="w/o wildfire",
                          th_var=0.08)
selection_factors_varimax(sigma_matrix = model_00_paper$sigma[[2]], 
                          treat_level="w/ wildfire",
                          th_var=0.08)
