# Set the working directory to the specified path
setwd('~/Documents/IMB/Johanna/')

# Source a custom R script named 'funs.R' from the 'model/scripts/' directory
source('model/scripts//funs.R')

{ 
  
  # Read the 'REAL_NETWORK.txt' file in TSV format into the 'huri' data frame
  huri <- read.delim('model/HURI/REAL_NETWORK.txt', sep = '\t')
  
  # Rename the columns in the 'huri' data frame to 'gene_a' and 'gene_b'
  colnames(huri) <- c('gene_a', 'gene_b')
  
  # Map ENSG edges to Pfam DDI types using the 'edges_ensg2pfam' function
  huri <- edges_ensg2pfam(huri, ensg_pfam)
  
  # Calculate the number of unique PPIs per each RDDI type within each Zgroup
  huri_thrs <- sapply(1:4, function(gr)
    length(unique(huri$PPI[huri$RDDI_type %in% whole$RDDI_type[whole$Zgroup == gr]])))
  
}


# whole table
huri_zscore <- merge(huri, whole[,c('RDDI_type',"X3did.z.score","Zgroup")], by = "RDDI_type", all.x = T)
write.csv(huri_zscore, file = "model/results/huri_3did.csv", row.names = F)

# unique ddis
huri_uddi <- unique(huri_zscore[,c('RDDI_type', 'PPI', 'Zgroup', 'X3did.z.score')])

# PPIs making the threshold
sum(with(huri_uddi[!is.na(huri_uddi$Zgroup) & huri_uddi$Zgroup > 2,], table(table(PPI))))

# PPIs with more than one DDI
sum(with(huri_uddi[!is.na(huri_uddi$Zgroup) & huri_uddi$Zgroup > 2,], table(table(PPI)))[2:5])

# DDIs per PPI

{
pdf('model/figs/PPIvsDDI.pdf', width = 2, height = 2)
# Set plotting parameters for the graph
#par(cex = 1.5)
par(mar = c(3, 3, 1, 0), mgp = c(2, .5, 0), cex = .9)
barplot(with(huri_uddi[!is.na(huri_uddi$Zgroup) & huri_uddi$Zgroup > 2,], table(table(PPI))), 
        ylim=c(0,1500), ann=T, axes=F, xlab = 'n DDIs', ylab='PPIs')
axis(2,at=seq(0,1500,250))
dev.off()
}

{# PPIs per DDI
pdf('model/figs/DDIvsPPI.pdf', width = 2, height = 2)
# Set plotting parameters for the graph
#par(cex = 1.5)
par(mar = c(3, 3, 1, 0), mgp = c(2, .5, 0), cex = .9)

hist(with(huri_uddi[!is.na(huri_uddi$Zgroup) & huri_uddi$Zgroup > 2,], table(RDDI_type)), 
     ann=T, axes=F, ylab = 'DDI types', xlab='n PPIs', main='', breaks = 50, ylim = c(0,350))
axis(2,at=seq(0,350,50))
axis(1,at=seq(0,400,100))
dev.off()
}

median(with(huri_uddi[!is.na(huri_uddi$Zgroup) & huri_uddi$Zgroup > 2,], table(RDDI_type)))

x<- head(sort(with(huri_uddi[!is.na(huri_uddi$Zgroup) & huri_uddi$Zgroup > 2,], table(RDDI_type)), decreasing = T))
cat(paste(names(x), x, sep=": ", collapse='\n\n'))

# Read a CSV file into the 'rdata' data frame
rdata <- read.csv('model/results/HURI_random_DDI_types_thr.csv')

# Calculate the cumulative sum of the reversed table of 'whole$Zgroup'
(thrtab <- cumsum(rev(table(whole$Zgroup))))
  
# Calculate the cumulative proportion 
thrtab / nrow(whole)
  
  
#scale by Z score
{
  # Create a PDF file for the plot
  pdf('model/figs/huri_z.pdf')
  
  # Set the layout of plots to 2x2
  par(mfrow = c(2, 2), cex = 1.1)
  
  # Plot a histogram of the standardized PPI count for Z score in the first Z score group
  hist(scale(rdata$PPI_gr1), xlim = c(-10, 110), probability = TRUE,
       main = paste('Z score <= ', thrs[1]), xlab = 'Standardized PPI count')
  
  # Add a red vertical line with the value in the sample
  abline(v = (huri_thrs[1] - mean(rdata$PPI_gr1)) / sd(rdata$PPI_gr1), col = 'red')
  
  # Calculate the empirical p-value
  mean(rdata$PPI_gr1 > huri_thrs[1])
  
  # Repeat the above process for the other groups
  hist(scale(rdata$PPI_gr2), xlim = c(-10,110), probability = T,
     main=paste(thrs[1],'< Z score <= ', thrs[2]),xlab = 'Standardized PPI count')
  abline(v=(huri_thrs[2]-mean(rdata$PPI_gr2))/sd(rdata$PPI_gr2),col='red')
  mean(rdata$PPI_gr2 > huri_thrs[2])

  hist(scale(rdata$PPI_gr3), xlim = c(-10,110), probability = T,
     main=paste(thrs[2],'< Z score <= ', thrs[3]),xlab = 'Standardized PPI count')
  abline(v=(huri_thrs[3]-mean(rdata$PPI_gr3))/sd(rdata$PPI_gr3),col='red')
  mean(rdata$PPI_gr3 > huri_thrs[3])

  hist(scale(rdata$PPI_gr4), xlim = c(-10,110), probability = T,
     main=paste(thrs[3],'< Z score'),xlab = 'Standardized PPI count')
  abline(v=(huri_thrs[4]-mean(rdata$PPI_gr4))/sd(rdata$PPI_gr4),col='red')
  mean(rdata$PPI_gr4 > huri_thrs[4])
dev.off()
}

#scale by DDI
{
  
  # Create a PDF file for the plot
  pdf('model/figs/huri_ppi_ddi.pdf')
  
  # Set the layout of plots to 2x2
  par(mfrow = c(2, 2), cex = 1.1)
  
  # Calculate the value of PPIs per DDI in the real HURI for the different Z score groups
  stat <- huri_thrs / table(whole$Zgroup)
  
  # Calculate the proportion of PPIs per DDI for the first Z group in the random networks
  samp <- rdata$PPI_gr1 / table(whole$Zgroup)[1]
  
  # Plot a histogram of the proportion of PPIs per DDI for Z score <= thrs[1]
  hist(samp, xlim = c(0, 0.8), probability = TRUE,
       main = paste('3did Z <= ', thrs[1]), xlab = 'PPIs per DDI')
  
  # Add a red vertical line at 'stat[1]' and display the calculated Z score
  abline(v = stat[1], col = 'red')
  text(0.6, 4, paste('Z=', round((stat[1] - mean(samp)) / sd(samp), 1)))
  
  # Repeat the above process for the remaining Z score groups
  
  samp <- rdata$PPI_gr2/table(whole$Zgroup)[2]
  hist(samp, xlim = c(0,.8), probability = T,
     main=paste(thrs[1],'< 3did Z <= ', thrs[2]),xlab = 'PPIs per DDI')
  abline(v=stat[2],col='red')
  text(0.6,10,paste('Z=',round((stat[2] - mean(samp)) / sd(samp) ,1)))

  samp <- rdata$PPI_gr3/table(whole$Zgroup)[3]
  hist(samp, xlim = c(0,.8), probability = T,
     main=paste(thrs[2],'< 3did Z <= ', thrs[3]),xlab = 'PPIs per DDI')
  abline(v=stat[3],col='red')
  text(0.6,5,paste('Z=',round((stat[3] - mean(samp)) / sd(samp) ,1)))

  samp <- rdata$PPI_gr4/table(whole$Zgroup)[4]
  hist(samp, xlim = c(0,.8), probability = T,
     main=paste(thrs[3],'< 3did Z'),xlab = 'PPIs per DDI')
  abline(v=stat[4],col='red')
  text(0.6,10,paste('Z=',round((stat[4] - mean(samp)) / sd(samp) ,1)))

dev.off()
}
