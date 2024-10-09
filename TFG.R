# BIOCONDUCTOR
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

# Biostrings
# BiocManager::install("Biostrings", force = TRUE)

# DADA2
# BiocManager::install("dada2", force = TRUE)

# phyloseq
# BiocManager::install("phyloseq")

# DESeq2
# BiocManager::install("DESeq2")

# microbiome
# BiocManager::install("microbiome")

# qiime2R
# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")

# ANCOMBC
# BiocManager::install("ANCOMBC") 

# microbiotaProcess
# BiocManager::install("MicrobiotaProcess")

# mongolite
install.packages("mongolite")


# Descarga automatizada de los datos de la BBDD

library(readr)
library(dplyr)
library(mongolite)

tabla <- read_delim("C:/ANTIGUA_D/TodoTFG/filereport_read_run_PRJEB34168_tsv.tsv", col_names = TRUE, delim ="\t")
head(tabla)

tabla_ordenada <- tabla[order(tabla$submitted_ftp), ]
head(tabla_ordenada)

enlaces <- strsplit(tabla_ordenada$submitted_ftp, ";")
# print(enlaces[[65]][1])

# Revisarlo da error porque no encuentra el servidor
for (i in 1:length(enlaces)){
  flt <- mongo("flights")
  muestras <- flt$import(gzcon(curl::curl(enlaces[i])))
}
