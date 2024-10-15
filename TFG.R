# BIOCONDUCTOR
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

# Biostrings
BiocManager::install("Biostrings", force = TRUE)

# DADA2
BiocManager::install("dada2", force = TRUE)

# phyloseq
BiocManager::install("phyloseq")

# DESeq2
BiocManager::install("DESeq2")

# microbiome
BiocManager::install("microbiome")

# qiime2R
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
  devtools::install_github("jbisanz/qiime2R")

# ANCOMBC
BiocManager::install("ANCOMBC") 

# microbiotaProcess
BiocManager::install("MicrobiotaProcess")

# mongolite
install.packages("mongolite")

# fastqcr
install.packages("fastqcr")

# svDialogs
install.packages("svDialogs")

library(readr)
library(dplyr)
library(Biostrings)
library(mongolite)
library(svDialogs)

# AUTOMATIZACIÓN DE LA DESCARGA DE LAS MUESTRAS DEL ENA
#   Esto sólo haría falta ejecutarlo una vez

directorio <- dlgInput(message = "Introduzca el directorio al tsv: ")$res

tabla <- read_delim(directorio, col_names = TRUE, delim ="\t")

tabla_ordenada <- tabla[order(tabla$submitted_ftp), ]

enlaces <- strsplit(tabla_ordenada$submitted_ftp, ";")

for (enl in enlaces){
  for (i in enl){
    destino <- strsplit(i, "/")[[1]]
    download.file(i,destino[length(destino)])
  }
}

# ABRIR LAS SECUENCIAS PARA EL PROCESAMIENTO

