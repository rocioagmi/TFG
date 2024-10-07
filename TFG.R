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
# install.packages("mongolite")


# Descarga automatizada de los datos de la BBDD

texto <- readLines("C:/Users/rocio/Documents/data/filereport_read_run_PRJEB34168_tsv.txt")
completo <- paste(texto, collapse = " ")

frases <- strsplit(completo, "\t")[[1]]
indices <- seq(17, length(frases), by = 9)
enlaces <- frases[indices]

print(enlaces)

flt <- mongo("flights")
flt$import(gzcon(curl::curl("enlace_descarga")))