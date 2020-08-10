# Code created by Matthew Elliott
# Contact: melliot1@ucsc.edu (valid until 2023)
# Phone:  231-392-1263   (just in case)


#############
### Set Up Notebook
#############

library(plyr)
library(biomaRt)

#############
### Load Data
#############

  # We reharmnonize the wholeblood dataset using harmony ran method
  
  # We load in the datasets from CSV files:
  raw_guinea= read.table("~/test_data/guinea_rna.csv", header=TRUE , sep =",",  stringsAsFactors=FALSE) #quote = input$quote,
  raw_tanzi = read.table("~/test_data/tanzania_rna.csv", header=TRUE , sep =",",  stringsAsFactors=FALSE) 
  raw_mali  = read.table("~/test_data/mali_rna.csv", header=TRUE , sep =",",  stringsAsFactors=FALSE)
  
  # we create the south africa dataset using resutls from the unharmonized_counts file
  # The south africa data is counts, so the unharmonized counts should have th correctly formatted data
  raw_all  = read.table("~/test_data/unharmonized_counts.csv", header=TRUE , sep =",",  stringsAsFactors=FALSE)
  raw_sa = raw_all[ , c(2,111:355) ] 
  raw_guinea[1:10,1:10]
  
  # check that mali is actually counts per million
  raw_mali[1:10,1:10]
  dim(raw_mali)
  colSums(raw_mali[,2:81] )
  # IT's currently not summing to 1 million
  dim(raw_sa)
  raw_sa[1:10,1:10]

#############
### compare to unharmonized_data.csv 
#############
  #Compare data agianst the unharmonized counts file
  length( unique( raw_guinea[,1] ))
  length(  raw_guinea[,1]  )
  
  
  raw_all  = read.table("~/test_data/unharmonized_counts.csv", header=TRUE , sep =",",  stringsAsFactors=FALSE)
  
  # check raw_all against raw_guinea
  dim(raw_guinea)
  dim(raw_tanzi)
  names(raw_all)
  raw_all[1:10,1:4]
  raw_guinea[1:10,1:3]
  raw_guinea[ which( raw_guinea$hgnc=="AARS" ), 1:3]
  # The datasets are the same just in different order
  
  # The mali dataset is changed
  names(raw_mali)
  raw_mali[1:10,1:3]
  names(raw_all)
  raw_all[1:10,c(2,31,32)]
  
  
#############
###    Remove Duplicates
#############
  
  # Guinea - Duplicates
  #length(unique(raw_guinea$hgnc))
  #length(raw_guinea$hgnc)
  raw_guinea2 = ddply( raw_guinea,"hgnc", numcolwise(sum) )

  # check that it works
  #length(unique(raw_guinea$hgnc))
  #length(raw_guinea$hgnc)
  #raw_guinea$hgnc[ which( duplicated(raw_guinea$hgnc)  ) ]
  #raw_guinea[ which(raw_guinea$hgnc=="RPS27") ,]
  #raw_guinea2[ which(raw_guinea2$hgnc=="RPS27") ,]
  
  
  # Tanzania - No duplicates
  length(unique(raw_tanzi$GeneSymbol))
  length(raw_tanzi$GeneSymbo)
  raw_tanzi[1:10,1:10]
  
  # South Africa - No duplicates
  #length( raw_sa$hgnc )
  #length( unique(raw_sa$hgnc ))
  
  # Mali - No duplicates
  #length( raw_mali$hgnc )
  #length( unique(raw_mali$hgnc ))
  
  
  

#############
### TPM Function
#############

  #counts = raw_guinea2 # for testing
  #counts = raw_tanzi # for testing
  
  
  
  # TPM Function
  TPM_Converter <- function( counts ) {
    # setup: Get gene lengths
    human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    
    # Get gene lengths for guinea
    gene_coords=getBM(attributes=c("hgnc_symbol", "start_position","end_position"), filters="hgnc_symbol", values=counts[,1], mart=human)
    gene_coords$size=gene_coords$end_position - gene_coords$start_position
    gene_coords$end_position = NULL
    gene_coords$start_position = NULL
    #dim(gene_coords2)
    #dim(counts)
    dim(counts_small)
    
    # Only keep genes with corresponding gene lengths
    gene_coords2 = gene_coords[ !duplicated(gene_coords$hgnc_symbol)  ,] # remove duplicates from gene sizes
    counts_small =  counts[ counts[,1] %in%  gene_coords2$hgnc_symbol ,] # only keep genes with gene sizes
    counts_small = counts_small[order(counts_small[,1]),] # order alphabetically
    gene_coords3 =  gene_coords2[ gene_coords2$hgnc_symbol %in%  counts_small[,1]  ,] # only keep genes with gene sizes
    gene_coords3 = gene_coords3[order(gene_coords3$hgnc_symbol),] # order alphabetically
    print( all( gene_coords3$hgnc_symbol == counts_small[,1] ) )
    
    # Create counts varibale
    counts_gene =  sweep( counts_small[,-1], MARGIN = 1, STATS = gene_coords3$size, FUN = "/")
    scaling_factor=colSums(counts_gene)/1000000
    TPM = sweep( counts_gene, MARGIN = 2, STATS = scaling_factor, FUN = "/")
    TPM_final = cbind( counts_small[,1], TPM )
    names(TPM_final)[1] = "hgnc"
    return(TPM_final)
    #counts[1:6,1:10] #head(gene_coords2$size) #141403.840 / 1491100 #TPM2[1:6,1:10]
    #counts2 = round( counts ) #counts2[1:10,1:10]
  }
  
#############
### Convert Data to TPM
#############  

  # Convert Guinea, Tanzania, and South Africa
  tpm_guinea = TPM_Converter(raw_guinea2)
  tpm_tanzi = TPM_Converter( raw_tanzi )
  tpm_sa = TPM_Converter( raw_sa )
  
  # Make datasets look slightly different
  # Guinea: shuffle rows and rename first column and add a duplicate row
  rows <- sample(nrow(tpm_guinea))
  tpm_guinea2 = tpm_guinea[rows,]
  names(tpm_guinea2)[1]="gene_symbol"
  tpm_guinea2 = rbind( tpm_guinea2, tpm_guinea2[1,] ) # duplicate a row
  
  
  
  # Tanzania: shuffle rows 
  rows <- sample(nrow(tpm_tanzi))
  tpm_tanzi2 = tpm_tanzi[rows,]

  #  Rescale Mali
  scaling_factor=colSums(raw_mali[,-1])/1000000
  tpms = sweep( raw_mali[,-1], MARGIN = 2, STATS = scaling_factor, FUN = "/")
  tpm_mali = cbind( raw_mali[,1], tpms )
  names(tpm_mali)[1] = "hgnc"
  
  #tpm_mali[1:20,1:10]
  #dim(raw_mali)
  #olSums( tpm_mali[,2:81] )
  #head(raw_mali[,-1])
  
  #colSums(tpm_guinea[,2:5])
  #colSums(tpm_tanzi[,3:6])
  #colSums(tpm_sa[,33:35])
  
  # Check it works
  #tpm_guinea[1:10,1:10]
  #dim(TPM_final)
  #TPM_final[1:10,1:10]
  #colSums(TPM_final[,-1])
  #safety=TPM_final
  

#############
### Save Datasets
#############  
  
  write.csv( tpm_tanzi2, "tpm_data/tanzania_tpm.csv", row.names=FALSE)
  write.csv( tpm_guinea2, "tpm_data/guinea_tpm.csv", row.names=FALSE)
  write.csv( tpm_mali, "tpm_data/mali_tpm.csv", row.names=FALSE)
  write.csv( tpm_sa, "tpm_data/south_africa_tpm.csv", row.names=FALSE)
  
  


  



