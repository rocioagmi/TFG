# To install and use Dada2, follow the online tutorial here: https://benjjneb.github.io/dada2/tutorial.html
# Load dada2
#install.packages("RCurl")
library(dada2); 
packageVersion("dada2")

# Set working directory 
setwd("ENTER PATH HERE")
path <- "ENTER PATH HERE"
list.files(path)

#Read files
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles 
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Filter and trim data using standard filtering parameters 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# On Windows set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,260),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, trimLeft = 20, multithread=FALSE) # On Windows set multithread=FALSE

# Learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Dereplication step 
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Apply dada2 algorithm (sample inference algorithm) - data must be dereplicated
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge paired reads
mergers <- NA
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE, justConcatenate=TRUE)

# Construct sequence table (amplicon sequence variant table) - higher version of OTU table prepared by traditional methods
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# Remove chimera sequences 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Assign taxonomy using reference databases (downloaded online)
taxa <- assignTaxonomy(seqtab.nochim, "ENTER PATH TO silva_nr_v128_train_set.fa", multithread=TRUE)
#taxa <- assignTaxonomy(seqtab.nochim, "ENTER PATH TO gg_13_8_train_set_97.fa.gz", multithread=TRUE)
#taxa <- assignTaxonomy(seqtab.nochim, "ENTER PATH TO rdp_train_set_14.fa.gz", multithread=TRUE)

#install.packages("tidyr")
library(tidyr)
library(data.table)
temp <- as.data.table(taxa)
temp_taxa <- paste(unlist(temp[,1]), unlist(temp[,2]), unlist(temp[,3]), unlist(temp[,4]), unlist(temp[,5]), unlist(temp[,6]), sep="_") 

### USE this code to output different ASV tables
# First generate taxonomy headers
kingdom <- unlist(temp[,1])
phylum <- paste(unlist(temp[,1]), unlist(temp[,2]), sep=";")
class <- paste(unlist(temp[,1]), unlist(temp[,2]), unlist(temp[,3]), sep=";")
order <- paste(unlist(temp[,1]), unlist(temp[,2]), unlist(temp[,3]), unlist(temp[,4]), sep=";")
family <- paste(unlist(temp[,1]), unlist(temp[,2]), unlist(temp[,3]), unlist(temp[,4]),  unlist(temp[,5]), sep=";")
genus <- paste(unlist(temp[,1]), unlist(temp[,2]), unlist(temp[,3]), unlist(temp[,4]), unlist(temp[,5]), unlist(temp[,6]), sep=";") 

# Next, generate ASV tables with labels according to each level of taxonomy
kingdom_table <- seqtab.nochim
kingdom_header <- c(kingdom)
colnames(kingdom_table) <- kingdom_header
write.csv(t(kingdom_table), "ENTER PATH TO SAVE CSV FILE")

phylum_table <- seqtab.nochim
phylum_header <- c(phylum)
colnames(phylum_table) <- phylum_header
write.csv(t(phylum_table), "ENTER PATH TO SAVE CSV FILE")

class_table <- seqtab.nochim
class_header <- c(class)
colnames(class_table) <- class_header
write.csv(t(class_table), "ENTER PATH TO SAVE CSV FILE")

order_table <- seqtab.nochim
order_header <- c(order)
colnames(order_table) <- order_header
write.csv(t(order_table), "ENTER PATH TO SAVE CSV FILE")

family_table <- seqtab.nochim
family_header <- c(family)
colnames(family_table) <- family_header
write.csv(t(family_table), "ENTER PATH TO SAVE CSV FILE")

genus_table <- seqtab.nochim
genus_header <- c(genus)
colnames(genus_table) <- genus_header
write.csv(t(genus_table), "ENTER PATH TO SAVE CSV FILE")

#PHYLOSEQ Analysis after Dada2 Analysis
#TO INSTAll
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')

install.packages("shiny") 
install.packages("shinythemes")
install.packages("networkD3")
install.packages("rmarkdown")
install.packages("markdown")
install.packages("png")
install.packages("ggplot2")
install.packages("httpuv")
install.packages("stringi")
install.packages("XML")
install.packages("IRanges")
install.packages("genefilter")

#Importing METADATA FILE IN CSV
#Change path to metadata file!
library(phyloseq)
metadata <- data.frame(read.csv(file = 'ENTER PATH TO METADATA CSV FILE', sep = ',', header = TRUE))
samples.out <- rownames(seqtab.nochim)
rownames(metadata) <- samples.out

#CREATING THE PHYLOSEQ OBJECTS
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxa))
save(ps, file = "SAVE .R File to be uploaded to Phyloseq GUI")

##Go to shiny-phyloseq
shiny::runGitHub("shiny-phyloseq","joey711")
