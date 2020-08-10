# Code created by Matthew Elliott
# Contact: melliot1@ucsc.edu (valid until 2023)
# Phone:  231-392-1263   (just in case)



#############
### Load Data
#############

  # We reharmnonize the wholeblood dataset using only the unahrmonized_counts file and mali_rna.csv
  # This is the simplest way to do it
  
  # load data
  raw_mali  = read.table("~/test_data/mali_rna.csv", header=TRUE , sep =",",  stringsAsFactors=FALSE)
  raw_all  = read.table("~/test_data/unharmonized_counts.csv", header=TRUE , sep =",",  stringsAsFactors=FALSE)
  
  
#############
### TPM Function
#############

  library(biomaRt)

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
    #dim(counts_s)
    
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
### Convert raw_all
#############
  
  raw_all2 = raw_all[,-1]
  tpm_all =  TPM_Converter( raw_all2 )

  
#############
### Combine with Mali data
#############
  
  # format mali data
  to_keep= raw_mali$hgnc %in% tpm_all$hgnc # only keep same rows
  mali_small = raw_mali[ to_keep, ]
  mali_small = mali_small[order(mali_small$hgnc),] # put data in same order
  #all( mali_small$hgnc == tpm_all$hgnc)
  
  # make mali columns sum to a million
  scaling_factor=colSums(mali_small[,-1])/1000000
  tpms = sweep( mali_small[,-1], MARGIN = 2, STATS = scaling_factor, FUN = "/")

  # Add to tpm_all
  tpm_all2 = tpm_all
  tpm_all2[, 30:109] = tpms
  
  # bring back in ensembl ids
  matches = raw_all[ raw_all$hgnc %in% tpm_all2$hgnc , 1:2 ]
  #all( matches$hgnc == tpm_all2$hgnc )
  tpm_final = cbind( matches$ensembl, tpm_all2 )
  names(tpm_final)[1] = "ensembl"
  
  #checkit = colSums(tpm_final[,c(-1,-2)])
  #sum( checkit -  1000000)
  
  #tpm_final[1:2,]
  #dim(raw_all)
  #dim(tpm_all2)
  #tpm_all2[1:5,1:5]
  #raw_all[1:5,1:5]
  #names(tpm_all2)
  #length(names(tpm_all2)[30:109])
  #dim(tpms)
  
  
  
#############
### Save Data
#############
  
  # we test the transformation we use for reformatting data for combat seq in harmony
  checkit = write.csv( tpm_final, "tpm_data/unharmonized_tpm.csv", row.names=FALSE)
  
  

  
#
  
  

  
#############
### Test Stuff
#############
  
  # we want to see what happens when we try to reformat the unharmonized counts files.
  # we test the transformation we use for reformatting data for combat seq in harmony rna
  checkit =  read.table("~/tpm_data/unharmonized_tpm.csv", header=TRUE , sep =",",  stringsAsFactors=FALSE)
  #checkit[1:10,1:10]
  check = checkit[,c(-1,-2)]
  check[1:2,]
  #check2 = round(check)
  #check2[1:10,1:10]
  
  #median(check)
  
  
  scaling_factor=colSums(check)/1000000 /1000
  tpms = sweep( check, MARGIN = 2, STATS = scaling_factor, FUN = "/")
  tpms = as.matrix( round(tpms) )

  mode(tpms) <- "integer" 
  max(tpms)
  colSums(tpms)
  
  
  scaling_factor=colSums(check)/1000000 
  back = sweep( check, MARGIN = 2, STATS = scaling_factor, FUN = "/")
  colSums(back)
  
  
  #colSums(tpms)
  #tpms[1:20,1:3]
  
  #
  
  
  
  
  
  
  
  
  

  