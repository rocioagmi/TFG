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

# Rqc
BiocManager::install("Rqc")


library(readr)
library(dplyr)
library(Biostrings)
library(ShortRead)
library(mongolite)
library(svDialogs)
library(TAF)
library(dada2)
library(Rqc)



# DESCARGA DE LOS DATOS
source("FUNC/Descargas_ENA.R")
nAcceso <- dlgInput(message = "Introduzca el número de acceso al proyecto ENA: ")$res
descargas_ENA(nAcceso)


# PREPROCESAMIENTO DE LOS DATOS
#   Separamos las muestras de personas MS de Healthy
muestrasMS <- sort(list.files("INPUT/DATA", pattern = "MS", full.names = TRUE))
muestrasHealty <- sort(list.files("INPUT/DATA", pattern = "Healthy", full.names = TRUE))

#   Separamos las hebras R1 y R2
MS_R1 <- muestrasMS[grepl("R1", muestrasMS)]
MS_R2 <- muestrasMS[grepl("R2", muestrasMS)]

Healty_R1 <- muestrasHealty[grepl("R1", muestrasHealty)]
Healty_R2 <- muestrasHealty[grepl("R2", muestrasHealty)]

# Obtenemos un objeto con reusltados relacionados con la calidad de las secuencias
rqc("INPUT/DATA", ".fastq.gz", pair = c(1,1), openBrowser = TRUE)


#   Visualizamos la calidad de las muestras MS y Healthy
plotQualityProfile(MS_R1[1:2])
plotQualityProfile(MS_R2[1:2]) # la calidad baja mucho

plotQualityProfile(Healty_R1[1:2])
plotQualityProfile(Healty_R2[1:2]) # la calidad baja mucho


# Filtramos las secuencias para mejorar la calidad
directorio <- mkdir("INPUT/DATA/FILTRADAS")

source("FUNC/FiltrarMuestras.R")
filtradasM <- filtrarMuestras(MS_R1, MS_R2)
filtradasH <- filtrarMuestras(Healty_R1, Healty_R2)

filtradasMS <- sort(list.files("INPUT/DATA/FILTRADAS", pattern = "MS", full.names = TRUE))
filtradasHealty <- sort(list.files("INPUT/DATA/FILTRADAS", pattern = "Healthy", full.names = TRUE))

filtradasMS_R1 <- filtradasMS[grepl("R1", filtradasMS)]
filtradasMS_R2 <- filtradasMS[grepl("R2", filtradasMS)]

filtradasH_R1 <- filtradasHealty[grepl("R1", filtradasHealty)]
filtradasH_R2 <- filtradasHealty[grepl("R2", filtradasHealty)]

plotQualityProfile(filtradasMS_R1[12:14])
plotQualityProfile(filtradasMS_R2[12:14])

plotQualityProfile(filtradasH_R1[12:14])
plotQualityProfile(filtradasH_R2[12:14])

# Dereplicación para eliminar posibles replicados de ADN
derepM_R1 <- derepFastq(filtradasMS_R1, verbose = TRUE)
derepM_R2 <- derepFastq(filtradasMS_R2, verbose = TRUE)

derepH_R1 <- derepFastq(filtradasH_R1, verbose = TRUE) 
derepH_R2 <- derepFastq(filtradasH_R2, verbose = TRUE)

names(derepM_R1) <- filtradasMS_R1
names(derepM_R2) <- filtradasMS_R2

names(derepH_R1) <- filtradasH_R1
names(derepH_R2) <- filtradasH_R2

# CERRAR ARCHIVOS

# ALMACENAMIENTO DE LOS DATOS EN UNA BASE DE DATOS MONGODB
# Abrir consola y escribir mongodb
datos_ENA <- mongo(collection = "datos_ENA", url = "mongodb://localhost:27017")

muestras <- list.files("INPUT/DATA", pattern = "\\.fastq\\.gz$", full.names = TRUE)

for (i in muestras) {
  contenido <- list(nombre = basename(i), contenido = i)
  datos_ENA$insert(contenido) 
}

datos_ENA$find(sort = '{"nombre": 1}')
datos_ENA$disconnect()
# Para usar el algoritmo dada() también tenemos que calcular la tasa de error