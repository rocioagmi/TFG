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
library(ggplot2)


# PRIMERO BÚSQUEDA PROGRAMÁTICA EN ENA
  # El estudio que escojamos en la búsqueda será el que se introduzca a continuación para 
  # la descarga automática de las muestras.

# DESCARGA DE LOS DATOS
source("FUNC/Descargas_ENA.R")
nAcceso <- dlgInput(message = "Introduzca el número de acceso al proyecto ENA: ")$res
descargas_ENA(nAcceso)


# PREPROCESAMIENTO DE LOS DATOS

listadoMuestras <- sort(list.files("INPUT/DATA", pattern = "\\.fastq\\.gz$", full.names = TRUE))

muestrasMS <- sort(list.files("INPUT/DATA", pattern = "MS", full.names = TRUE))
muestrasHealthy <- sort(list.files("INPUT/DATA", pattern = "Healthy", full.names = TRUE))

MS_R1 <- muestrasMS[grepl("R1", muestrasMS)]
MS_R2 <- muestrasMS[grepl("R2", muestrasMS)]

Healthy_R1 <- muestrasHealthy[grepl("R1", muestrasHealthy)]
Healthy_R2 <- muestrasHealthy[grepl("R2", muestrasHealthy)]

  # Informe de calidad
dir.create("OUTPUT/REPORT")
dir.create("OUTPUT/FIGURES")
dir.create("OUTPUT/RDS")

source("FUNC/InformeCalidad.R")
informeCalidad(listadoMuestras)

source("FUNC/GraficosCalidad.R")
graficosCalidad(MS_R1,MS_R2)
graficosCalidad(Healthy_R1,Healthy_R2)


  # Filtrado DADA2
source("FUNC/FiltrarMuestras.R")
filtrarMuestras(MS_R1, MS_R2)
filtrarMuestras(Healthy_R1, Healthy_R2)

listadoFiltrado <- sort(list.files("OUTPUT/FILTRADO", pattern = "\\.fastq\\.gz$", full.names = TRUE))

filtradasMS <- sort(list.files("OUTPUT/FILTRADO", pattern = "MS", full.names = TRUE))
filtradasHealthy <-  sort(list.files("OUTPUT/FILTRADO", pattern = "Healthy", full.names = TRUE))

filtradoR1 <- listadoFiltrado[grepl("R1", listadoFiltrado)]
filtradoR2 <- listadoFiltrado[grepl("R2", listadoFiltrado)]

filtradasMS_R1 <- filtradasMS[grepl("R1", filtradasMS)]
filtradasMS_R2 <- filtradasMS[grepl("R2", filtradasMS)]

filtradasH_R1 <- filtradasHealthy[grepl("R1", filtradasHealthy)]
filtradasH_R2 <- filtradasHealthy[grepl("R2", filtradasHealthy)]

# Evalúa la calidad de las muestras ya filtradas
source("FUNC/InformeCalidad.R")
informeCalidad(listadoFiltrado)

source("FUNC/GraficosCalidad.R")
graficosCalidad(filtradasMS_R1, filtradasMS_R2)
graficosCalidad(filtradasHealthy_R1, filtradasHealthy_R2)


  # Filtrado SHORTREAD
dir.create("OUTPUT/FILTRADOSR")
source("FUNC/FiltrarMuestras.R")
filtradoSR(listadoMuestras)

filtradoSR <- sort(list.files("OUTPUT/FILTRADOSR", pattern = "\\.fastq\\.gz$", full.names = TRUE))

source("FUNC/InformeCalidad.R")
informeCalidad(filtradoSR)

filtradosrMS <- sort(list.files("OUTPUT/FILTRADOSR", pattern = "MS", full.names = TRUE))
filtradosrH <-  sort(list.files("OUTPUT/FILTRADOSR", pattern = "Healthy", full.names = TRUE))

filtradosrR1 <- filtradoSR[grepl("R1", filtradoSR)]
filtradosrR2 <- filtradoSR[grepl("R2", filtradoSR)]

filtradosrMS_R1 <- filtradosrMS[grepl("R1", filtradosrMS)]
filtradosrMS_R2 <- filtradosrMS[grepl("R2", filtradosrMS)]

filtradosrH_R1 <- filtradosrH[grepl("R1", filtradosrH)]
filtradosrH_R2 <- filtradosrH[grepl("R2", filtradosrH)]

source("FUNC/GraficosCalidad.R")
graficosCalidad(filtradosrMS_R1, filtradosrMS_R2)
graficosCalidad(filtradosrH_R1, filtradosrH_R2)


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
dir.create("INPUT/BB_DD")

source("FUNC/AsignarTaxonomia.R")
asignarTaxonomia(tablaSinQuim)

taxHit <- assignTaxonomy(tabSinQuim, "INPUT/BB_DD/hitdb_v1.00.fa.gz", multithread = TRUE)
saveRDS(taxHit, paste0("OUTPUT/RDS", "/taxHit.Rds"))

taxRDP_S <- assignTaxonomy(tabSinQuim, "INPUT/BB_DD/rdp_19_toSpecies_trainset.fa.gz", multithread = TRUE)
saveRDS(taxRDP_S, paste0("OUTPUT/RDS", "/taxRDP_S.Rds"))

taxRDP_G <- assignTaxonomy(tabSinQuim, "INPUT/BB_DD/rdp_19_toGenus_trainset.fa.gz", multithread = TRUE)
saveRDS(taxRDP_G, paste0("OUTPUT/RDS", "/taxRDP_G.Rds"))

taxSilva <- assignTaxonomy(tabSinQuim, "INPUT/BB_DD/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread = TRUE)
saveRDS(taxSilva, paste0("OUTPUT/RDS", "/taxSilva.Rds"))



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
