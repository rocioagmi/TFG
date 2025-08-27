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

# QuasR
BiocManager::install("QuasR")

# svDialogs
install.packages("svDialogs")

# TAF
install.packages("TAF")

# ShortRead
BiocManager::install("ShortRead", force = TRUE)

# Rqc
BiocManager::install("Rqc")

# DT
install.packages("DT")

# BiocParallel
BiocManager::install("BiocParallel")

library(readr)
library(dplyr)
library(Biostrings)
library(ShortRead)
library(mongolite)
library(svDialogs)
library(TAF)
library(dada2)
library(QuasR)
library(DT)
library(Rqc)
library(BiocParallel)

# PRIMERO BÚSQUEDA PROGRAMÁTICA EN ENA
  # El estudio que escojamos en la búsqueda será el que se introduzca a continuación para 
  # la descarga automática de las muestras.

# DESCARGA DE LOS DATOS
source("FUNC/Descargas_ENA.R")
nAcceso <- dlgInput(message = "Introduzca el número de acceso al proyecto ENA: ")$res
descargas_ENA(nAcceso)

# PREPROCESAMIENTO DE LOS DATOS

  # Informe de calidad
dir.create("OUTPUT/REPORT")

source("FUNC/InformeCalidad.R")
directorioMuestras <- dir("INPUT/DATA", "\\.fastq\\.gz$", full = TRUE)
informeCalidad(directorioMuestras)

  # Filtrado

listadoMuestras <- sort(list.files("INPUT/DATA", pattern = "\\.fastq\\.gz$", full.names = TRUE))

muestrasMS <- sort(list.files("INPUT/DATA", pattern = "MS", full.names = TRUE))
muestrasHealty <- sort(list.files("INPUT/DATA", pattern = "Healthy", full.names = TRUE))

MS_R1 <- muestrasMS[grepl("R1", muestrasMS)]
MS_R2 <- muestrasMS[grepl("R2", muestrasMS)]

Healthy_R1 <- muestrasHealthy[grepl("R1", muestrasHealthy)]
Healthy_R2 <- muestrasHealthy[grepl("R2", muestrasHealthy)]

  # --- PACKAGE DADA2 ---
dir.create("OUTPUT/FILTRADO")

source("FUNC/FiltrarMuestras.R")
filtrarMuestras(MS_R1, MS_R2)
filtrarMuestras(Healthy_R1, Healthy_R2)

filtradasMS <- sort(list.files("OUTPUT/FILTRADO", pattern = "MS", full.names = TRUE))
filtradasHealthy <- sort(list.files("OUTPUT/FILTRADO", pattern = "Healthy", full.names = TRUE))

filtradasMS_R1 <- filtradasMS[grepl("R1", filtradasMS)]
filtradasMS_R2 <- filtradasMS[grepl("R2", filtradasMS)]

filtradasH_R1 <- filtradasHealthy[grepl("R1", filtradasHealthy)]
filtradasH_R2 <- filtradasHealthy[grepl("R2", filtradasHealthy)]

  # Confirmar la calidad
plotQualityProfile(filtradasMS_R1[1:10])
plotQualityProfile(filtradasMS_R2[1:10])

plotQualityProfile(filtradasH_R1[1:10])
plotQualityProfile(filtradasH_R2[1:10])

  # --- PACKAGE SHORTREAD ---
source("FUNC/FiltrarMuestras.R")
filtradoSR(listadoMuestras)

  # Evalúa la calidad de las muestras ya filtradas
source("FUNC/InformeCalidad.R")
directorioFiltradas <- dir("OUTPUT/FILTRADO", "\\.fastq\\.gz$", full = TRUE)
informeCalidad(directorioFiltradas)


  # REVISAR READ COUNTS, DEDUPLICACIÓN, CUTADAPT

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
