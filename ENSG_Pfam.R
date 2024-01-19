# map ENSG to Pfam
ensg_uniprot <- read.csv('data/ENSG_UniProt.csv')
uniprot_pfam <- read.csv('data/UniProt_Pfam.csv')
ensg_pfam <- merge(ensg_uniprot, uniprot_pfam,
                   by.x = 'uniprot_id',
                   by.y = 'UniProt_ID')[,c('ensembl_id', 'uniprot_id', 'hgnc_symbol',
                                           'hgnc_id', 'Pfam_ID', 'Start', 'End', 'DomainName')]
write.csv(unique.data.frame(ensg_pfam[,c('ensembl_id','Pfam_ID')]), file = 'data/ENSG_Pfam.csv', row.names = F)
