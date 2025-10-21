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

# httr
install.packages("httr")

# readr
install.packages("readr")

# Jsonlite
install.packages("jsonlite")


library(readr)
library(httr)
library(jsonlite)
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
library(ggplot2)

# PRIMERO BÚSQUEDA PROGRAMÁTICA EN ENA - IMPORTANTE CAMBIAR ESTO
source("FUNC/Busqueda_ENA.R")
construirConsultaENA()

# DESCARGA DE LOS DATOS
source("FUNC/Descargas_ENA.R")
nAcceso <- dlgInput(message = "Introduzca el número de acceso al proyecto ENA: ")$res
descargas_ENA(nAcceso)


# PREPROCESAMIENTO DE LOS DATOS

listadoMuestras <- sort(list.files("INPUT/DATA", pattern = "\\.fastq\\.gz$", full.names = TRUE))

muestrasMS <- sort(list.files("INPUT/DATA", pattern = "MS", full.names = TRUE))
muestrasHealthy <- sort(list.files("INPUT/DATA", pattern = "Healthy", full.names = TRUE))

muestrasR1 <- sort(list.files("INPUT/DATA", pattern = "R1", full.names = TRUE))
muestrasR2 <- sort(list.files("INPUT/DATA", pattern = "R2", full.names = TRUE))


  # Informe de calidad
dir.create("OUTPUT/REPORT")
dir.create("OUTPUT/FIGURES")
dir.create("OUTPUT/RDS")
dir.create("OUTPUT/FILTRADO")

source("FUNC/InformeCalidad.R")
informeCalidad(listadoMuestras)

source("FUNC/GraficosCalidad.R")
graficosCalidad(muestrasR1, muestrasR2)


  # Filtrado DADA2
source("FUNC/FiltrarMuestras.R")
filtrarMuestras(MS_R1, MS_R2)
filtrarMuestras(Healthy_R1, Healthy_R2)

listadoFiltrado <- sort(list.files("OUTPUT/FILTRADO", pattern = "\\.fastq\\.gz$", full.names = TRUE))

filtradasMS <- sort(list.files("OUTPUT/FILTRADO", pattern = "MS", full.names = TRUE))
filtradasHealthy <-  sort(list.files("OUTPUT/FILTRADO", pattern = "Healthy", full.names = TRUE))

filtradoR1 <- listadoFiltrado[grepl("R1", listadoFiltrado)]
filtradoR2 <- listadoFiltrado[grepl("R2", listadoFiltrado)]

# Evalúa la calidad de las muestras ya filtradas
source("FUNC/InformeCalidad.R")
informeCalidad(listadoFiltrado)

source("FUNC/GraficosCalidad.R")
graficosCalidad(filtradoR1, filtradoR2)


# PROCESAMIENTO DE LOS DATOS 

  # Aplica la función dada 
source("FUNC/FlujoTrabajoDada.R")
union <- flujoTrabajoDada(filtradoR1, filtradoR2)


  # Construye la tabla de secuencias
source("FUNC/ConstruirTablaSecuencias.R")
construirTablaSecuencias(union)
tablaSecuencias <- "OUTPUT/RDS/seqtab.Rds"


  # Elimina las quimeras
source("FUNC/EliminarQuimeras.R")
eliminarQuimeras(tablaSecuencias)
tablaSinQuim <- "OUTPUT/RDS/tabSinQuim.Rds"


  # Asignación taxonómica
# C:/ANTIGUA_D/TodoTFG/BB_DD/hitdb_v1.00.fa.gz
# C:/ANTIGUA_D/TodoTFG/BB_DD/rdp_19_toGenus_trainset.fa.gz /rdp_19_toSpecies_trainset.fa.gz
# C:/ANTIGUA_D/TodoTFG/BB_DD/silva_nr99_v138.2_toGenus_trainset.fa.gz

source("FUNC/AsignarTaxonomia.R")
asignarTaxonomia(tablaSinQuim)







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
