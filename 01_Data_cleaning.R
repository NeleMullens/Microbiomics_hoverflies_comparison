#For amplicon sequencing of V3-V4 16S rRNA genes (Bacteria), we will use DADA2 
#to infer exact sequencing variants from the data. 


#Wnhen needed, install pacakges first: 
#phyloseq and dada2
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
BiocManager::install("phyloseq")
BiocManager::install("dada2", version = "3.12")

#Load packages and dependencies
library(dada2)
library(phyloseq)
library(ggplot2)

# Begin DADA2 pipeline -------------------------------------------
## Start quality check --------------------------------------------

#First we read in the names of the fastq files, and perform some string
#manipulation to get lists of the forward and reverse fastq files in matched
#order

#Sort the reads (Forward and Reverse) so they are in the same order
fnFs <- sort(list.files("data", pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files("data", pattern="_2.fastq.gz", full.names = TRUE))

#check quality of reads with `rplotQualityProfile`
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])


## Filtering and trimming --------------------------------- 
#Extract sample names (take subset ('[') from basename, when the basename string
#is split at "_")
sample.names <- sapply(strsplit(basename(fnFs), "DNA_"), `[`, 1)
#Alternatively, use purr: 
##unlist(purrr::map(strsplit(basename(fnFs),"_"),purrr::pluck,1))

#Prepare filenames/path for filtered data
filtFs <- file.path("data/Filtered_data", paste0(sample.names, "F_filtered.fq.gz"))
filtRs <- file.path("data/Filtered_data", paste0(sample.names, "R_filtered.fq.gz"))
#necessary for the dereplication and merging step 
names(filtFs) <- sample.names 
names(filtRs) <- sample.names

#For the filtering, we’ll use standard filtering parameters: maxN=0 (DADA2
#requires no Ns), truncQ=2 and rm.phix=TRUE. maxEE will be set at 3, whereas
#standard this is set on 2. Here we will use a value of 3 to better compare with
#the Uparse pipeline (where a value of 3 is chosen as a standard). The maxEE
#parameter sets the maximum number of “expected errors” allowed in a read, which
#is a better filter than simply averaging quality scores.

#trimLeft is used to remove primers. Don't use external trimming programmes, as
#they remove primer mismatches, but due to the high number of unknown bacteria,
#you will have these in your dataset and you would thus remove actual biological
#variation

#r filtering_trimming #runtime of 2-3 hours?
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(17, 21), truncLen=c(260,210),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE) # On non-Windows set multithread=TRUEhead(out)
head(out)
dim(out)
plotQualityProfile(filtFs[1:6])
plotQualityProfile(filtRs[1:6])

## Estimate the error rate --------------------------------- 

#The DADA2 algorithm depends on a parametric error model (err) and every
#amplicon dataset has a different set of error rates. The learnErrors method
#learns the error model from the data, by alternating estimation of the error
#rates and inference of sample composition until they converge on a jointly
#consistent solution. As in many optimization problems, the algorithm must begin
#with an initial guess, for which the maximum possible error rates in this data
#are used (the error rates if only the most abundant sequence is correct and all
#the rest are errors).

#runtime around 45 minutes for each one
#>**Note** Parameter learning is computationally intensive, so by default the
#*learnErrors* function uses only a subset of the data (the first 1M reads). If
#the plotted error model does not look like a good fit, try increasing the
#nreads parameter to see if the fit improves.
errF <- learnErrors(filtFs, nread=1e6, multithread=TRUE)
errR <- learnErrors(filtRs, nread=1e6, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

## Dereplication  ---------------------------------
#Dereplication optinal?
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


## Sample Inference  --------------------------------- 
#With dereplication step done
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#If dereplication is skipped #without dereplication compute time still okay
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Inspect the returned dada-class object
dadaFs[[1]]


## Merge paired reads ---------------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])



## Construct sequence table --------------------------------------------------
#make a sequence table of the merged reads with the counts per sample
seqtab <- makeSequenceTable(mergers)

#Check table
dim(seqtab)
sum(seqtab)

#Inspect distribution of sequence length
table(nchar(getSequences(seqtab)))

#save the sequence table as an RDS file, change the name according to the run
saveRDS(seqtab, "outputs/seqtab_Microbiome_hoverflies.rds")


## Remove chimeras -----------------------------------------------------------

#this step takes around 1 hour to run
#seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim2 <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE, minFoldParentOverAbundance=2)
dim(seqtab.nochim2)
sum(seqtab.nochim2)

#table(nchar(getSequences(seqtab.nochim)))
table(nchar(getSequences(seqtab.nochim2)))

sum(seqtab.nochim2)/sum(seqtab)

saveRDS(seqtab.nochim2, "outputs/seqtab.nochim_All.rds")

## Taxonomy assignment ---------------------------------------------------------

taxa <- assignTaxonomy(seqtab.nochim2, "data/SILVA_138.1_SSURef_NR99_tax_silva.0.fa", multithread=FALSE)

input_taxa <- read.csv("outputs/Input_taxonomic_assignment.csv", header=TRUE, row.names = 1)
input_taxa<-as.data.frame(input_taxa)

input_taxa_matrix <- data.matrix(t(input_taxa), rownames.force = NA)
head(colnames(input_taxa_matrix))

taxa_assigned_extra <- assignTaxonomy(input_taxa_matrix, "data/SILVA_data_DNA.fa", multithread=TRUE)

#Show taxa in dataset
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

table<-as.data.frame(t(seqtab.nochim2))
write.table(table,file="outputs/Table_with_reads_and_taxonomic_data_updated.txt",sep="\t",col.names=TRUE,row.names=TRUE)

table<-as.data.frame(t(taxa))
write.table(table,file="outputs/Table_with_bacteria.txt",sep="\t",col.names=TRUE,row.names=TRUE)

# End of DADA2 pipeline --------------------------------------------------

# Combine taxonomy and data in one file -------------------------------
dim(seqtab.nochim2)
dim(taxa)

#transpose table
table<-as.data.frame(t(seqtab.nochim2))
table<-as.data.frame(t(input_taxa_matrix))


#remove rownames
rownames(table)<-NULL

#add extra columns containing taxa, for this, first transform the taxa object to a data frame
taxadf<-as.data.frame(taxa_assigned)
table$Kingdom<-taxadf$Kingdom
table$Phylum<-taxadf$Phylum
table$Class<-taxadf$Class
table$Order<-taxadf$Order
table$Family<-taxadf$Family
table$Genus<-taxadf$Genus

#add extra column containing sequences
table$seq<-as.vector(colnames(seqtab.nochim2))
table$seq<-as.vector(colnames(input_taxa_matrix))


#write table to output file
write.table(table,"outputs/Full_table_with_taxonmy_other_ref.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
#write.xlsx(table,"outputs/Full_table_with_taxonmy.xlsx",col.names=TRUE,row.names=FALSE)

# Remove contaminations from wet lab procedure --------------------------------

install.packages("devtools")
library(microDecon)

# transform Full_table_with_taxonmy.txt: add a first column with OTUs, and the negative control(s) should be next. Remove taxonomic data (or reduce to one column, at the end)
clean_table <- read.csv("outputs/Table_to_clean_negs.csv")
clean_tab <- as.data.frame(clean_table)
View(clean_tab)
# if error, check if there are empty columns and remove these through clean_tab$X <- NULL
contamdone<-decon(data=clean_tab, numb.blanks=5, numb.ind=c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4), taxa=F) #all 96 samples in groups of 4; based on location, sex, and treatment
contamdone2<-contamdone$decon.table
contamdone2  #with this step we can see which ASV have been removed due to too low counts in samples compared to the negative control and how many counts have been subtracted for the ASVs that have not been removed.

cleaned<-as.data.frame(contamdone2)
write.table(cleaned,file="outputs/Filtered_OTUs.txt",sep="\t",col.names=TRUE,row.names=TRUE)


#Rows and columns are swithed, so switch them around:
library(tibble)
switched <- as_tibble(t(cleaned[,-1]))
write.table(switched,file="outputs/Filtered_OTUs.txt",sep="\t",col.names=TRUE,row.names=TRUE)


#### manually filter out Mitochondria and Chloroplast from file!
#Output file: outputs/Complete_filtered_data.csv
#Output file without negative control: Complete_filtered_data_neg_removed.csv


# check for which samples to discard -----------------------------

disard_data <- read.csv("outputs/Samples_discard_check.csv", row.names=1)
hist(disard_data$Reads)
quantile(disard_data$Reads)
quantile(disard_data$Reads, probs=0.05)
#Remove all samples with reads below 852
# samples removed:
# KTC_BB_03 (218)
# SGS_HB_06 (265)
# KTS_BC_03 (451)
# KDB_BB_01 (458)
# CMB_BA_04 (517)
# HUC_HB_04 (527)
# RUB_BC_02 (558)
# KTW_BA_04 (612)
# HUB_HA_05 (630)
# MKS_BE_05 (637)
# results in outputs/Complete_filtered_data.csv

# Create a frequency table --------------------------------
prop_table <- read.csv("outputs/Complete_filtered_data.csv", row.names=1)
View(prop_table)
#The last seven columns are for taxonomic identification, we want to select all but that
samplenumber<-ncol(prop_table)-7
#The samples are all columbs starting from 1 to the number calculated above
samples<-colnames(prop_table)[1:samplenumber]


#Place the samples in a dataframe
freq_table <- data.frame(matrix(nrow = nrow(prop_table), ncol=samplenumber, 
                                dimnames = list(c(rownames(prop_table)), 
                                                c(colnames(prop_table[1:samplenumber])))))

#calculate frequencies
for (i in 1:samplenumber){
  freq_table[,i] <- prop_table[,i]/colSums(prop_table[1:samplenumber])[i]
}
# [,i] => select all elements from column i

write.table(freq_table, "outputs/freq_table.txt", sep="\t")
write.xlsx(freq_table, file = "outputs/freq_table.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)


#save(freq_table, "outputs/freq_table.RData" )

#add species information to the frequency table again
freq_table$Kingdom<-taxadf$Kingdom
freq_table$Phylum<-taxadf$Phylum
freq_table$Class<-taxadf$Class
freq_table$Order<-taxadf$Order
freq_table$Family<-taxadf$Family
freq_table$Genus<-taxadf$Genus

#write to new output table
write.table(freq_table,"outputs/frequency_table_complete.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.xlsx(freq_table, file = "outputs/frequency_table_taxonomy.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

--------------------------------------------------------------------------------
  ##                                                                            ##
  ##                  Data now ready for further analysis                       ##
  ##                                                                            ##
  --------------------------------------------------------------------------------
  
  
  
  
  
# calculate phylogenetic tree ------------------------------------
install.packages("GUniFrac")
install.packages("phangorn")
install.packages("magrittr")
install.packages("ade4")
install.packages('xlsx')
install.packages("ALDEx2")

library(GUniFrac)
library(phangorn)
library(magrittr)
library(ade4)
library(xlsx)
library("ALDEx2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("DECIPHER")
BiocManager::install("rPlant")


seqtab.nochim2 <- readRDS("outputs/Tussenfiles/seqtab.nochim_All.rds")

# Write all sequences to a fasta file with the ID's corresponding to the rownumber
uniquesToFasta(seqtab.nochim2,"Sequences_Nov2022.fasta",ids = (1:ncol(seqtab.nochim2)))

# Read the fasta file and put it in a variable
microbiome_phylo.fasta<-Biostrings::readAAStringSet("outputs/Input_tree_reconstruction.fasta")

#An alignment is performed and processors are set to five to spead up the proces 
# but not overload the server. If you have enough processors available you can 
# elevate the number
Aligned_seq<-DECIPHER::AlignSeqs(microbiome_phylo.fasta,processors = 5)

# The aligned sequences are written to a file to be used by FASTREE
# The file is written to the current working directory
Biostrings::writeXStringSet(Aligned_seq,"outputs/Microbiome_Syrph_Apis_aligned.fasta")




# names of file are now wrong:
install.packages("phylotools")
library("phylotools")

old_name <- c("seq1", "seq2", "seq3", "seq5", "Seq9", "seq12")
new_name <- c("Magnolia", "Ranunculus", "Carex", "Morus", "Ulmus", "Salix")

ref2 <- data.frame(old_name, new_name)

# old_name <- get.fasta.name("outputs/Microbiome_Syrph_Apis_aligned.fasta")
rename <- read.csv("outputs/Rename_fasta_file.csv", header=T)
ref_name <- data.frame(rename$ï..old_name, rename$new_name)
rename.fasta(infile = "outputs/Microbiome_Syrph_Apis_aligned.fasta", ref_table = ref_name, outfile = "outputs/Microbiome_Syrph_Apis.fasta")

## now we can continue with correct names
# use fasttree on the HPC, very fast
# Alternatively, use this tree with CIPRES for RAxML analysis, Geneious for consensus


# Normalize the data -----------------------------



### MetagenomeSeq --------------------------------
#address the effects of both normalization and under-sampling of microbial communities
#packages:
library(vegan)
library(dplyr)
library(ggplot2)
library(ggpubr)
BiocManager::install("metagenomeSeq")
library(metagenomeSeq)
library(readxl)

raw_data_mbiom <- read.csv("outputs/Complete_filtered_data.csv", row.names=1)
# samplenumber<-ncol(raw_data_mbiom)-7
# #The samples are all columbs starting from 1 to the number calculated above
# samples<-colnames(raw_data_mbiom)[1:samplenumber]

#Place the samples in a dataframe

drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "seq")
abund_dat <- raw_data_mbiom[ , !(names(raw_data_mbiom) %in% drops)]


#raw_mbiome <- as.matrix(abund_dat)
meta_mbiom <- read.csv("outputs/Metadata_hoverflies.csv", row.names=1)

#create a MRexperiment object, input annotated dataframe and OTU table
phenotypeData = AnnotatedDataFrame(meta_mbiom)
OTUdata = AnnotatedDataFrame(abund_dat)
obj = newMRexperiment(abund_dat,phenoData=phenotypeData,featureData=OTUdata)

#calculate the proper percentile by which to normalize counts
p = cumNormStatFast(obj)
# result is 0.5

# calculate the scaling factors
calc_scaling = cumNorm(obj, p = p)

# to save normalised metrix
mat = MRcounts(calc_scaling, norm = TRUE, log = TRUE)
# exportMat(mat, file = "outputs/MetagenomeSeq/Metadata_hoverflies_normalised.csv")

mat_df <- as.data.frame(mat)
#write.table(mat_df,file="outputs/MetagenomeSeq/Microbiome_hoverflies_normalised.txt")
write.csv(mat_df,file="outputs/MetagenomeSeq/MetagenomeSeq_normalised_all_data.csv")



### crl central log transformed --------------------------------
## alternative normalisation: crl central log transformed 
library("Hotelling")
library(ALDEx2)
raw_microbiome <- read.csv("outputs/Complete_filtered_data.csv", row.names=1)

drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "seq")
raw_microbiome <- raw_microbiome[ , !(names(raw_microbiome) %in% drops)]

raw_microbiome_clr <- aldex.clr(raw_microbiome, mc.samples = 500, denom=rclr, verbose=TRUE)

ald_h <-lapply(X=raw_microbiome_clr@analysisData, FUN=rowMeans)
ald_h <-as.data.frame(ald_h)
#write.table(ald_h,file="outputs/clr_transformed_data.txt",sep="\t",col.names=TRUE,row.names=TRUE)
write.csv(ald_h,file="outputs/CentralLogTransformed/clr_normalised_data.csv")

# add taxonomic info to clr normalisation
taxa_data <- 
  raw_microbiome %>%
  rownames_to_column(var = "otus") %>%
  select(otus,all_of(drops))

ald_h_taxa <-
  ald_h %>% 
  rownames_to_column(var = "otus") %>% 
  left_join(taxa_data)

write.csv(ald_h_taxa,file="outputs/clr_normalised_taxa_data.csv")

# visualise to see if join worked well
ald_h_taxa %>% select(drops) %>% visdat::vis_dat()
ald_h_taxa %>% select(drops) %>% visdat::vis_miss()



### normalisation with vegan ----------------

raw_data_mbiom <- read.csv("outputs/Complete_filtered_data.csv", row.names=1)
meta_mbiom <- read.csv("outputs/Metadata_hoverflies.csv", row.names=1)

drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "seq")
raw_data_mbiom <- raw_data_mbiom[ , !(names(raw_data_mbiom) %in% drops)]
dat <- as.matrix(raw_data_mbiom)

vegan_freq <- decostand(dat, "frequency")

write.table(vegan_freq, file="outputs/Vegan_frequency/vegan_normfreq.txt", sep="\t", quote=F, col.names=NA)


## --------------Redo everything with aggregated genera. First, generate genera file ------------

GENUS_abundances<-read.csv("outputs/Complete_filtered_data_genera.csv", row.names=1)
dim(GENUS_abundances)
GENUS_abundances$Genus<-as.factor(GENUS_abundances$Genus)
length(levels(GENUS_abundances$Genus))
GENUS_abundances$Family<-as.factor(GENUS_abundances$Family)
length(levels(GENUS_abundances$Family))
GENUS_abundances$Phylum<-as.factor(GENUS_abundances$Phylum)
length(levels(GENUS_abundances$Phylum))
# 22002 OTUs representing 843 genera from 248 families and 30 phyla (including NAs pooled for each family, clade or phylum)

# rowsum (pooling of genera: remove text columns from the matrix!)
ncol(GENUS_abundances)
colnames(GENUS_abundances)
temp_GENUS_abundances<-rowsum(GENUS_abundances[,1:182], GENUS_abundances$Genus, na.rm=TRUE)

#(transpose and) save
# temp_GENUS_abundances<-t(temp_GENUS_abundances)
dim(temp_GENUS_abundances)
# saveRDS(temp_GENUS_abundances, file="GENERA_aggregated.rds")
write.csv(temp_GENUS_abundances, file="outputs/GENERA_aggregated.csv")

# rowsum (pooling of genera: remove text columns from the matrix!)
ncol(GENUS_abundances)
colnames(GENUS_abundances)
temp_FAM_abundances<-rowsum(GENUS_abundances[,1:182], GENUS_abundances$Family, na.rm=TRUE)

#(transpose and) save
# temp_GENUS_abundances<-t(temp_GENUS_abundances)
dim(temp_GENUS_abundances)
# saveRDS(temp_GENUS_abundances, file="GENERA_aggregated.rds")
write.csv(temp_FAM_abundances, file="outputs/FAMILIES_aggregated.csv")

# rowsum (pooling of genera: remove text columns from the matrix!)
ncol(GENUS_abundances)
colnames(GENUS_abundances)
temp_PHY_abundances<-rowsum(GENUS_abundances[,1:182], GENUS_abundances$Phylum, na.rm=TRUE)

#(transpose and) save
# temp_GENUS_abundances<-t(temp_GENUS_abundances)
dim(temp_GENUS_abundances)
# saveRDS(temp_GENUS_abundances, file="GENERA_aggregated.rds")
write.csv(temp_PHY_abundances, file="outputs/PHYLA_aggregated.csv")

## proportional normalisation ---------------------------

# Create a frequency table --------------------------------
prop_table <- read.csv("outputs/GENERA_aggregated.csv", row.names=1)
View(prop_table)
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "seq")
prop_table <- prop_table[ , !(names(raw_microbiome) %in% drops)]
#The last seven columns are for taxonomic identification, we want to select all but that
samplenumber<-ncol(prop_table)
#The samples are all columbs starting from 1 to the number calculated above
samples<-colnames(prop_table)[1:samplenumber]


#Place the samples in a dataframe
freq_table <- data.frame(matrix(nrow = nrow(prop_table), ncol=samplenumber, 
                                dimnames = list(c(rownames(prop_table)), 
                                                c(colnames(prop_table[1:samplenumber])))))

#calculate frequencies
for (i in 1:samplenumber){
  freq_table[,i] <- prop_table[,i]/colSums(prop_table[1:samplenumber])[i]
}
# [,i] => select all elements from column i

library(xlsx)
write.xlsx(freq_table, file = "outputs/freq_table_genera.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)



## central log clr normalisation ---------------------------

raw_microbiome <- read.csv("outputs/GENERA_aggregated.csv", row.names=1)
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "seq")
raw_microbiome <- raw_microbiome[ , !(names(raw_microbiome) %in% drops)]

raw_microbiome_clr <- aldex.clr(raw_microbiome, mc.samples = 500, denom=rclr, verbose=TRUE)

ald_h <-lapply(X=raw_microbiome_clr@analysisData, FUN=rowMeans)
ald_h <-as.data.frame(ald_h)

write.xlsx(ald_h, file = "outputs/CentralLogTransformed/clr_normalised_genera.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)



## MetagenomeSeq normalisation ---------------------------
abund_dat <- read.csv("outputs/GENERA_aggregated.csv", row.names=1)
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "seq")
abund_dat <- abund_dat[ , !(names(abund_dat) %in% drops)]

#raw_mbiome <- as.matrix(abund_dat)
meta_mbiom <- read.csv("outputs/Metadata_hoverflies.csv", row.names=1)

#create a MRexperiment object, input annotated dataframe and OTU table
phenotypeData = AnnotatedDataFrame(meta_mbiom)
OTUdata = AnnotatedDataFrame(abund_dat)
obj = newMRexperiment(abund_dat,phenoData=phenotypeData,featureData=OTUdata)

#calculate the proper percentile by which to normalize counts
p = cumNormStatFast(obj)

# calculate the scaling factors
calc_scaling = cumNorm(obj, p = p)

# to save normalised metrix
mat = MRcounts(calc_scaling, norm = TRUE, log = TRUE)
# exportMat(mat, file = "outputs/MetagenomeSeq/Metadata_hoverflies_normalised.csv")

mat_df <- as.data.frame(mat)
write.xlsx(mat_df, file = "outputs/MetagenomeSeq_normalised_genera.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)


### normalisation with vegan
raw_data_mbiom <- read.csv("outputs/GENERA_aggregated.csv", row.names=1)
meta_mbiom <- read.csv("outputs/Metadata_hoverflies.csv", row.names=1)

drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "seq")
raw_data_mbiom <- raw_data_mbiom[ , !(names(raw_data_mbiom) %in% drops)]
dat <- as.matrix(raw_data_mbiom)

vegan_freq <- decostand(dat, "frequency")

write.xlsx(vegan_freq, file = "outputs/Vegan_frequency/vegan_normfreq_genera.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)