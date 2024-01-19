setwd('~/Documents/IMB/Johanna/')
library(readxl)

# dataset
feats <- read.csv(file = 'data/3did_Interactome_ProtCID_AF_IUPred_metrics_mc_set_for_ML_without_synthetic_constructs.csv')

whole <- read.csv("data/3did_Interactome_ProtCID_metrics_inter_hetero_subset_without_synthetic_constructs_representativity.csv")

{
pdf(file = 'model/figs/representativity.pdf')
par(mfrow=c(4,3),cex=.7)

#interchain fraction
for (x in names(whole)[3:ncol(whole)]) {
  pop <- whole[,x]
  sam <- feats[,x]
  
  tmp <- sapply(1:1e5,function(x)mean(sample(pop,nrow(feats)),na.rm = T)) 
  hist(tmp,breaks=20, main=paste('p value (two-tailed):', 
                                 round(mean(abs(mean(pop,na.rm = T) - mean(sam,na.rm = T)) 
                                            < abs(mean(pop,na.rm = T)-tmp)),3)), 
       xlab = x)
  abline(v=mean(sam,na.rm=T),col='red')
  
  
}


dev.off()
}
