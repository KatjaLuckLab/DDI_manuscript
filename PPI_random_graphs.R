# Set the working directory to the specified path
setwd('~/Documents/IMB/Johanna/')

# Source a custom R script named 'funs.R' from the 'model/scripts/' directory
source('model/scripts/funs.R')

# Load the 'doParallel' library for parallel computing
library(doParallel)

# Threshold analysis

# Setup backend to use multiple processors
totalCores = detectCores()

{
  # Print the current timestamp
  timestamp()
  
  # Leave one core to avoid overloading your computer
  cluster <- makeCluster(totalCores[1] - 1)
  registerDoParallel(cluster)
  
  # Run a for loop in parallel to process files in the 'model/HURI/huri' directory
  rdata <- foreach(i = list.files('model/HURI/huri', full.names = TRUE), combine = rbind) %dopar% {
    # Read a random network file and map ENSG edges to Pfam DDI types using the 'edges_ensg2pfam' function
    random <- edges_ensg2pfam(read.delim(i, sep = '\t'), ensg_pfam)
    
    # Calculate the average number of unique PPIs per RDDI type within each Zgroup
    c(
      strsplit(i, '/')[[1]][4],
      sapply(1:4, function(gr)
        length(unique(random$PPI[random$RDDI_type %in% whole$RDDI_type[whole$Zgroup == gr]]))
      )
    )
  }
  
  # Stop the parallel cluster
  stopCluster(cluster)
  
  # Print the current timestamp (time taken for parallel processing)
  timestamp()  # took 16 minutes
}

# Combine the results from the parallel processing into a single data frame
rdata <- as.data.frame(do.call(rbind, rdata))

# Set column names for the resulting data frame
colnames(rdata) <- c('file', 'PPI_gr1', 'PPI_gr2', 'PPI_gr3', 'PPI_gr4')

# Write the resulting data frame to a CSV file
write.csv(rdata, file = 'model/results/HURI_random_DDI_types_thr.csv', row.names = FALSE)