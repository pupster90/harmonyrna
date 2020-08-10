# Code created by Matthew Elliott
# Contact: melliot1@ucsc.edu (valid until 2023)
# Phone:  231-392-1263   (just in case)

# We check how well the harmony RNA method works
# We take Zichengs dataset as a test and harmonize it various ways. 
# We then measure the difference in results
# Notes: I checked what happend if I ran combat-seq twice on counts. It produced identical results


########################
### Set Up notebook  ###
########################

# for formatting plots
attach(mtcars)

# NOTE: Remember to run the original combat_seq before doing this
addResourcePath("www", paste(getwd() , "/www", sep="") ) # Have shiny recognize folder with files
source("www/helper_seq.R") # Helper functions for Combat Seq
source("www/combat_seq.R")


# Load dataset
unharmonized_TPM = read.table( "tpm_data/Zicheng_TPM.csv", stringsAsFactors=FALSE, header=TRUE, sep=',')

# Batches derived from code in other file
batches = c( rep(1, 30), rep(2, 27) )



###########################
### Make Counts and TPM ###
###########################


#  Here we convert a TPM to a count file. Then we run both in combat-seq to see if we can get the same output
# Useful Tutorial: https://www.youtube.com/watch?v=TTUrtCY2k-w

# Get gene lengths
library(biomaRt)
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords=getBM(attributes=c("hgnc_symbol", "start_position","end_position"), filters="hgnc_symbol", values=unharmonized_TPM[,1], mart=human)
gene_coords$size=gene_coords$end_position - gene_coords$start_position
gene_coords$end_position = NULL
gene_coords$start_position = NULL
#length(unharmonized_TPM[,1])
#dim(gene_coords)

# Only keep genes with corresponding gene lengths
gene_coords2 = gene_coords[ !duplicated(gene_coords$hgnc_symbol)  ,] # remove duplicates from gene sizes
TPM2 =  unharmonized_TPM[ unharmonized_TPM[,1] %in%  gene_coords2$hgnc_symbol ,] # only keep genes with gene sizes
TPM2 = TPM2[order(TPM2$gene_ID),] # order alphabetically
gene_coords2 = gene_coords2[order(gene_coords2$hgnc_symbol),] # order alphabetically
#all( gene_coords2$hgnc_symbol == TPM2$gene_ID )

# Reformat TPM
#TPM2[1:10,1:10]
TPM3 = as.matrix(TPM2[,-1])
TPM3[ TPM3<0 ] = 0 

# Create counts varibale
counts =  sweep( TPM2[,-1], MARGIN = 1, STATS = gene_coords2$size, FUN = "*")
counts2 = round( counts ) #counts2[1:10,1:10]
counts3 =  as.matrix( counts2 )
mode(counts3) <- "integer" 
counts3[ counts3<0 ] = 0 # Remove Negative Values

#colSums(TPM2[-1])


######################
### Run Combat-Seq ###
######################

# Counts
start_time <- Sys.time()
harmonized_counts= ComBat_seq( counts3, batch=batches, group=NULL )
counts_time = Sys.time() - start_time

# TPM * 10
start_time <- Sys.time()
TPM_10 = round( TPM3 * 10 )
mode(TPM_10) <- "integer" 
harmonized_TPM_10 = ComBat_seq( TPM_10, batch=batches, group=NULL )
TPM_10_time = Sys.time() - start_time

# TPM * 100
start_time <- Sys.time()
TPM_100 = round( TPM3 * 100 )
mode(TPM_100) <- "integer" 
harmonized_TPM_100 = ComBat_seq( TPM_100, batch=batches, group=NULL )
TPM_100_time = Sys.time() - start_time

# TPM * 1000
start_time <- Sys.time()
TPM_1000 = round( TPM3 * 1000 )
mode(TPM_1000) <- "integer" 
harmonized_TPM_1000 = ComBat_seq( TPM_1000, batch=batches, group=NULL )
TPM_1000_time = Sys.time() - start_time


# TPM * 100000
start_time <- Sys.time()
TPM_100000 = round( TPM3 * 100000 )
mode(TPM_1000) <- "integer" 
harmonized_TPM_100000 = ComBat_seq( TPM_100000, batch=batches, group=NULL )
TPM_100000_time = Sys.time() - start_time


# TPM * 1000000
start_time <- Sys.time()
TPM_1000000 = round( TPM3 * 1000000 )
mode(TPM_1000) <- "integer" 
#harmonized_TPM_1000000 = ComBat_seq( TPM_1000000, batch=batches, group=NULL )
# Causes computer to crash




########################################
###  Convert Harmonization to TPM
########################################


# Counts
harmonized_counts_rpk = sweep( harmonized_counts, MARGIN = 1, STATS = gene_coords2$size, FUN = "/")
scaling_factor=colSums(harmonized_counts_rpk)/1000000
harmonized_TPM_counts = sweep( harmonized_counts_rpk, MARGIN = 2, STATS = scaling_factor, FUN = "/")
#harmonized_counts_TPM[1800:1810,1:6] #colSums(harmonized_counts_TPM)

# TPM * 10
scaling_factor=colSums(harmonized_TPM_10)/1000000
harmonized_TPM_10 = sweep( harmonized_TPM_10, MARGIN = 2, STATS = scaling_factor, FUN = "/")

# TPM * 100
scaling_factor=colSums(harmonized_TPM_100)/1000000
harmonized_TPM_100 = sweep( harmonized_TPM_100, MARGIN = 2, STATS = scaling_factor, FUN = "/")


# TPM * 1000
scaling_factor=colSums(harmonized_TPM_1000)/1000000
harmonized_TPM_1000 = sweep( harmonized_TPM_1000, MARGIN = 2, STATS = scaling_factor, FUN = "/")

# TPM * 100000
scaling_factor=colSums(harmonized_TPM_100000 )/1000000
harmonized_TPM_100000  = sweep( harmonized_TPM_100000 , MARGIN = 2, STATS = scaling_factor, FUN = "/")


# TPM * 1000000
#scaling_factor=colSums(harmonized_TPM_1000000)/1000000
#harmonized_TPM_1000000 = sweep( harmonized_TPM_1000000, MARGIN = 2, STATS = scaling_factor, FUN = "/")


#######################
### Compare Results: easy ###
#######################


max
which(gene_names == "R")
#row_num=45
row_num=6
gene_names[6]

# gene_coords[ which( min(gene_coords2$size) == gene_coords2$size) ,]
#which( gene_names == "LENG8" )
#row_num = 578

gene_names[row_num]

plot( harmonized_TPM_counts[row_num,], main="", ylab="Transcripts per Million", xlab="Sample Index", col="blue" , pch=19 )
plot( harmonized_TPM_10[row_num,], main="AGBL2 Expression TPM", ylab="Counts per Million", xlab="Sample Number", col =c("red","blue","orange") , pch=19  )
plot( harmonized_TPM_1000[row_num,], main="Leng8 Expression TPM", ylab="Transcripts per Million", xlab="Sample Number", col =c("red","blue","orange") , pch=19  )
plot( harmonized_TPM_1000000[row_num,], main="Expression", ylab="Counts per Million", xlab="Sample Number", col =c("red","blue","orange") , pch=19  )



#################################
### Compare Results: rigorous ###
#################################


#dev.new()
#op <- par(no.readonly = TRUE)
#dev.off()
#op
#par(mfrow(1,1) )

# Note: We show that the na's are caused when result is guessed 100% accurately
gene_names=  TPM2[,1]
diffs = abs( harmonized_TPM_10 - harmonized_TPM_counts )
sds =  apply( harmonized_TPM_counts , 1, sd )
normalized = sweep( diffs, MARGIN = 1, STATS = sds , FUN = "/")

all(diffs[ which(sds==0) ,] == 0 ) # 100% correct where sds leads to NA
length(diffs[ which(sds==0) ,])    # length of the instances
sum(is.na(normalized))             # is the same length as the NA's



#   TPM * 10
diffs = abs( harmonized_TPM_10 - harmonized_TPM_counts )
sds =  apply(  harmonized_TPM_counts , 1, sd )
normalized = sweep( diffs, MARGIN = 1, STATS = sds , FUN = "/")
normalized[which(is.na(normalized))] = 0
median( normalized )
mean( normalized )
sd(normalized)
length( which(normalized<.01) )/ length(normalized)
#length( which(normalized<.05) )/ length(normalized)
#length( which(normalized<1) )/ length(normalized)
#hist(normalized, breaks=200, main="Histogram of Normalized \n Variation  between Methods", xlab="Deviation in SD units")
#hist(normalized, xlim=c(0,3),ylim=c(0,2000), breaks=4000, main="Histogram of Normalized \n Variation  between Methods", xlab="Deviation in SD units")
#hist(normalized, xlim=c(4,8),ylim=c(0,1000), breaks=100, main="Histogram of Normalized \n Variation  between Methods", xlab="Deviation in SD units")


#   TPM * 100
diffs = abs( harmonized_TPM_100 - harmonized_TPM_counts )
sds =  apply( harmonized_TPM_counts, 1, sd )
normalized = sweep( diffs, MARGIN = 1, STATS = sds , FUN = "/")
normalized[which(is.na(normalized))] = 0
median( normalized )
mean( normalized )
sd(normalized)
length( which(normalized<.01) )/ length(normalized)
#length( which(normalized<.005) )/ length(normalized)


#   TPM * 1000
diffs = abs( harmonized_TPM_1000 - harmonized_TPM_counts )
sds =  apply( harmonized_TPM_counts, 1, sd )
normalized = sweep( diffs, MARGIN = 1, STATS = sds , FUN = "/")
normalized[which(is.na(normalized))] = 0
median( normalized )
mean( normalized )
sd(normalized)
length( which(normalized<.01) )/ length(normalized)
#length( which(normalized<.005) )/ length(normalized)



#   TPM * 100000
diffs = abs( harmonized_TPM_100000 - harmonized_TPM_counts )
sds =  apply( harmonized_TPM_counts, 1, sd )
normalized = sweep( diffs, MARGIN = 1, STATS = sds , FUN = "/")
normalized[which(is.na(normalized))] = 0
median( normalized )
mean( normalized )
sd(normalized)
length( which(normalized<.01) )/ length(normalized)




##########################
###   Simulate Plots   ###
##########################

# Note: We run everything above only for the TPM*10 example to create these plots

# Use nonparametrics to to create a PDF for simulating error
#plot( density(normalized), ylim=c(0,1),xlim=c(0,3) )
sim = hist(normalized, xlim=c(0,3),ylim=c(0,2000), breaks=4000, main="Histogram of Normalized \n Variation  between Methods", xlab="Deviation in SD units")

Plot_Simulation <- function( num = 6, color="red" ) {
  # Plot OG data
  gene= gene_names[num]
  data = harmonized_TPM_counts[num,]
  ylim = range(data)
  real_plot = plot( data,xaxt='n',  main=paste(gene,"Counts"), ylim=ylim,  ylab="Counts per Million", xlab="Sample Number", col = color, pch=19 )
  
  # Plot TPM data
  data = harmonized_TPM_1000[num,]
  real_plot = plot( data, yaxt='n',xaxt='n',  main=paste(gene,"TPM"), ylim=ylim,  ylab="Counts per Million", xlab="Sample Number", col = color, pch=19 )
  
  
  # plot simulated result
  sim_error = sample( sim$mids, length(data), replace=TRUE, prob=sim$density/sum(sim$density)  )
  sim_error2 = sim_error * sd(data) * sample(c(-1,1), length(data), replace=TRUE )
  sim_data = data + sim_error2
  plot( sim_data, yaxt='n',xaxt='n', main=paste(gene,"Simulated"), ylim=ylim, ylab="", xlab="", col = color, pch=19 )

}

Plot_Simulations <- function( nums  ) {
  par(mfrow=c(4,3),
      oma = c(5,4,0,0) + 0.3,
      mar = c(1,0,1,1) + 0.3
  )
  Plot_Simulation(nums[1], color="blue")
  Plot_Simulation(nums[2], color="purple")
  Plot_Simulation(nums[3], color="orange")
  Plot_Simulation(nums[4], color="red")
  par(mfrow=c(1,1))
}



Plot_Simulations( num=c(10,50,500,1170) )



#######################
### Create Mutant   ###
#######################

# Test mutant dataset:Make 1 file tpm, the other is counts.  aka: "Does TPM matter?"
# Yes, mutant does make a differrence
TPM4 = round( TPM3 * 1000 )
mode(TPM4) <- "integer" 
mutant = cbind( TPM4[,which(batches==1)] , counts3[,which(batches==2)] )
#head(mutant)

# Combat on Mutant 
harmonized_mutant= ComBat_seq( mutant, batch=batches, group=NULL )

#  Rescale Mutant 
#which( max(gene_coords2$size) == gene_coords2$size)
# gene_coords[6,]
harmonized_mutant_rpk = harmonized_mutant
#harmonized_mutant_rpk[,which(batches==2)] = sweep( harmonized_mutant[,which(batches==2)], MARGIN = 1, STATS = gene_coords2$size, FUN = "/")
scaling_factor=colSums(harmonized_mutant_rpk)/1000000
harmonized_mutant_TPM = sweep( harmonized_mutant_rpk, MARGIN = 2, STATS = scaling_factor, FUN = "/")
harmonized_mutant_TPM[1800:1810,1:6]

# Plot mutant
plot( harmonized_mutant_TPM[row_num,], main="Leng8 Expression Counts/TPM", ylab="", xlab="Sample Number", col =c("red","blue","orange") , pch=19  )



#########################
###   Scratch Paper   ###
#########################





##########################




