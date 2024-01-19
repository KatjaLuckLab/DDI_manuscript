# Define a function to map ENSG edges to Pfam DDI types
edges_ensg2pfam <- function(x, ensg_pfam) {
  # Merge the input dataset 'x' with the 'ensg_pfam' dataset based on gene_a and gene_b columns
  x_mapped <- merge(
    merge(x, ensg_pfam[, c('ensembl_id', 'Pfam_ID')], by.x = 'gene_a', by.y = 'ensembl_id'),
    ensg_pfam[, c('ensembl_id', 'Pfam_ID')],
    by.x = 'gene_b',
    by.y = 'ensembl_id'
  )
  
  # Create a new column 'RDDI_type' by pasting the sorted 'Pfam_ID' values with underscores
  x_mapped$RDDI_type <- apply(x_mapped[, c('Pfam_ID.x', 'Pfam_ID.y')],
                              1,
                              function(y) paste(sort(y), collapse = '_')
  )
  
  # Create a new column 'PPI' by pasting the sorted 'gene_a' and 'gene_b' values with underscores
  x_mapped$PPI <- apply(x_mapped[, c('gene_a', 'gene_b')],
                        1,
                        function(y) paste(sort(y), collapse = '_')
  )
  
  return(x_mapped)
}

# Define a vector of threshold values
thrs <- c(2.3, 4.47, 5.3)

# Read the 'ENSG_Pfam' dataset from a CSV file
ensg_pfam <- read.csv('data/ENSG_Pfam.csv')

# Read the '3did_Interactome_ProtCID_metrics_inter_hetero_subset_without_synthetic_constructs_representativity' dataset from a CSV file
whole <- read.csv("data/3did_Interactome_ProtCID_metrics_inter_hetero_subset_without_synthetic_constructs_representativity.csv")

# Check that Pfam sorting is consistent in the 'whole' dataset
# Create a new column 'RDDI_type' by splitting the 'DDI.type' column and pasting sorted values with underscores
whole$RDDI_type <- sapply(strsplit(whole$DDI.type, '_'), function(y) paste(sort(y), collapse = '_'))

# Create a new column 'Zgroup' in the 'whole' dataset to split it in groups by Z score.
# The 'Zgroup' is determined based on whether the 'X3did.z.score' is greater than threshold values in 'thrs'
whole$Zgroup <- apply(sapply(thrs, function(x) whole$X3did.z.score > x), 1, sum) + 1
# Note: The 'Zgroup' is incremented by 1 to ensure that the group labels start from 1.