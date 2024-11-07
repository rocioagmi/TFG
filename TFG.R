# BIOCONDUCTOR
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

# Biostrings
BiocManager::install("Biostrings", force = TRUE)

# DADA2
BiocManager::install("dada2")

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

# TAF
install.packages("TAF")

# ShortRead
BiocManager::install("ShortRead", force = TRUE)

library(readr)
library(dplyr)
library(Biostrings)
library(ShortRead)
library(mongolite)
library(svDialogs)
library(TAF)
library(dada2)




# DESCARGA DE LOS DATOS
source("FUNC/Descargas_ENA.R")
nAcceso <- dlgInput(message = "Introduzca el número de acceso al proyecto ENA: ")$res
descargas_ENA(nAcceso)

# ALMACENAMIENTO DE LOS DATOS EN UNA BASE DE DATOS MONGODB
datos_ENA <- mongo(collection = "datos_ENA", url = "mongodb://localhost:27017")

muestras <- list.files("C:/ANTIGUA_D/TFG/INPUT/DATA", pattern = "\\.fastq\\.gz$", full.names = TRUE)

for (i in muestras) {
  contenido <- list(nombre = basename(i), contenido = i)
  datos_ENA$insert(contenido) 
}

datos_ENA$find(sort = '{"nombre": 1}')
datos_ENA$disconnect()

# PREPROCESAMIENTO DE LOS DATOS
#   Separamos las muestras de personas MS de Healthy
muestrasMS <- sort(list.files("INPUT/DATA", pattern = "MS", full.names = TRUE))
muestrasHealty <- sort(list.files("INPUT/DATA", pattern = "Healthy", full.names = TRUE))

#   Separamos las hebras R1 y R2
MS_R1 <- muestrasMS[grepl("R1", muestrasMS)]
MS_R2 <- muestrasMS[grepl("R2", muestrasMS)]

Healty_R1 <- muestrasHealty[grepl("R1", muestrasHealty)]
Healty_R2 <- muestrasHealty[grepl("R2", muestrasHealty)]

#   Visualizamos la calidad de las muestras MS y Healthy
plotQualityProfile(MS_R1[1:2])
plotQualityProfile(MS_R2[1:2]) # la calidad baja mucho

plotQualityProfile(Healty_R1[1:2])
plotQualityProfile(Healty_R2[1:2]) # la calidad baja mucho

# Filtramos las secuencias para mejorar la calidad
directorio <- mkdir("INPUT/DATA/FILTRADAS")
nMS_R1 <- paste("INPUT/DATA", strsplit(basename(MS_R1), "/M"),sep = "/FILTRADAS")
nMS_R2 <- paste("INPUT/DATA", strsplit(basename(MS_R2), "/M"),sep = "/FILTRADAS")

nH_R1 <- paste("INPUT/DATA", strsplit(basename(Healty_R1), "/H"),sep = "/FILTRADAS")
nH_R2 <- paste("INPUT/DATA", strsplit(basename(Healty_R2), "/H"),sep = "/FILTRADAS")

#   Vamos a obtener por dónde cortar las secuencias según la posición en la que empieza a disminuir la caliadad por debajo de un humbral
fqMS_R1 <- lapply(MS_R1, readFastq)
# qMS_R1 <- lapply(fqMS_R1, quality)
trimmed_MSR1 <- list()
for (i in 1:length(fqMS_R1)){
  trimmed_MSR1[[i]] <- trimEnds(fqMS_R1[[i]], "I") 
}



filterAndTrim(MS_R1, nMS_R1, MS_R2, nMS_R2, truncLen = c(240,240),
              maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
              compress = TRUE, trimLeft = 20, multithread = FALSE)

plotQualityProfile(nMS_R1[1:2])
plotQualityProfile(nMS_R2[1:2])