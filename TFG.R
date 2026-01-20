# BIOCONDUCTOR
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

# Biostrings
#BiocManager::install("Biostrings", force = TRUE)

# DADA2
#BiocManager::install("dada2")

# phyloseq
#BiocManager::install("phyloseq")

# DESeq2
#BiocManager::install("DESeq2")

# microbiome
#BiocManager::install("microbiome")

# qiime2R
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#  devtools::install_github("jbisanz/qiime2R")

# ANCOMBC
#BiocManager::install("ANCOMBC") 

# microbiotaProcess
#BiocManager::install("MicrobiotaProcess")

# mongolite
#install.packages("mongolite")

# QuasR
#BiocManager::install("QuasR")

# svDialogs
#install.packages("svDialogs")

# TAF
#install.packages("TAF")

# ShortRead
#BiocManager::install("ShortRead", force = TRUE)

# Rqc
#BiocManager::install("Rqc")

# DT
#install.packages("DT")

# BiocParallel
#BiocManager::install("BiocParallel")

# httr
#install.packages("httr")

# httr2
#install.packages("httr2")

# readr
#install.packages("readr")

# Jsonlite
#install.packages("jsonlite")

# stringr
#install.packages("stringr")

# ===============================================
# LIBRERÍAS
# ===============================================

library(httr)
library(jsonlite)
library(svDialogs)
library(DT)
library(readr)
library(ShortRead)
library(dada2)
#library(dplyr)
#library(Biostrings)
#library(mongolite)
#library(TAF)
#library(QuasR)
#library(Rqc)
#library(BiocParallel)
#library(ggplot2)



# ================================================
# BÚSQUEDA PROGRAMÁTICA (EBI/ENA)
# ================================================
source("FUNC/Busqueda_ENA.R")
muestrasEBI <- construirConsulta()


# ================================================
# DESCARGA DE LOS DATOS
# ================================================
source("FUNC/Descargas_ENA.R")
nAcceso <- dlgInput(message = "Introduzca el número de acceso al proyecto ENA: ")$res
descargas_ENA(nAcceso)


# ================================================
# PREPROCESAMIENTO DE LOS DATOS
# ================================================

# Obtiene la lista ordenada de todos los archivos .fastq.gz en la carpeta de entrada
listadoMuestras <- sort(list.files("INPUT/DATA", pattern = "\\.fastq\\.gz$", full.names = TRUE))

# Separa las muestras por tipo: MS (enfermos) y Healthy (sanos)
muestrasMS <- sort(list.files("INPUT/DATA", pattern = "MS", full.names = TRUE))
muestrasHealthy <- sort(list.files("INPUT/DATA", pattern = "Healthy", full.names = TRUE))

# Separa los archivos por lectura: R1 (forward) y R2 (reverse)
muestrasR1 <- sort(list.files("INPUT/DATA", pattern = "R1", full.names = TRUE))
muestrasR2 <- sort(list.files("INPUT/DATA", pattern = "R2", full.names = TRUE))

# -------------------------------------------------
# CREACIÓN DE CARPETAS DE SALIDA
# -------------------------------------------------
# Crea directorios para organizar los resultados: informes, gráficos, datos filtrados, ...
dir.create("OUTPUT/REPORT")
dir.create("OUTPUT/FIGURES")
dir.create("OUTPUT/RDS")
dir.create("OUTPUT/FILTRADO")

# -------------------------------------------------
# INFORME DE CALIDAD DE SECUENCIACIÓN
# -------------------------------------------------
# Carga la función para generar informes FastQC
source("FUNC/InformeCalidad.R")
# Genera informe de calidad para todas las muestras
informeCalidad(listadoMuestras)

# Carga función para gráficos de calidad 
source("FUNC/GraficosCalidad.R")
# Genera gráficos de calidad para lecturas R1 (forward) y R2 (reverse)
graficosCalidad(muestrasR1, muestrasR2)

# ----------------------------------------------
# FILTRADO DE CALIDAD CON DADA2
# ----------------------------------------------
# Carga función de filtrado 
source("FUNC/FiltrarMuestras.R")

# Aplica filtrado a muestras MS y Healthy por separado
filtrarMuestras(MS_R1, MS_R2)
filtrarMuestras(Healthy_R1, Healthy_R2)

# Lista todos los archivos filtrados generados
listadoFiltrado <- sort(list.files("OUTPUT/FILTRADO", pattern = "\\.fastq\\.gz$", full.names = TRUE))

# Clasifica archivos filtrados por tipo de muestra
filtradasMS <- sort(list.files("OUTPUT/FILTRADO", pattern = "MS", full.names = TRUE))
filtradasHealthy <-  sort(list.files("OUTPUT/FILTRADO", pattern = "Healthy", full.names = TRUE))

# Separa lecturas forward y reverse de los archivos filtrados
filtradoR1 <- listadoFiltrado[grepl("R1", listadoFiltrado)]
filtradoR2 <- listadoFiltrado[grepl("R2", listadoFiltrado)]

# -------------------------------------------
# INFORME DE CALIDAD POST FILTRADO
# -------------------------------------------
# Reutiliza funciones para evaluar mejora tras filtrado
source("FUNC/InformeCalidad.R")
informeCalidad(listadoFiltrado)

source("FUNC/GraficosCalidad.R")
graficosCalidad(filtradoR1, filtradoR2)


# ===========================================
# PROCESAMIENTO DE LOS DATOS 
# ===========================================

# -------------------------------------------
# FLUJO PRINCIPAL DADA2
# -------------------------------------------

# Carga la función que ejecuta: learnErrors, derep, dada, mergePairs
source("FUNC/FlujoTrabajoDada.R")
# Procesa lecturas filtradas y devuelve la tabla de conteos de secuencias
union <- flujoTrabajoDada(filtradoR1, filtradoR2)

# ------------------------------------------
# CONSTRUCCIÓN DE LA TABLA DE SECUENCIAS -- ESTE PASO IGUAL LO JUNTO AL ANTERIOR
# ------------------------------------------
# Genera y guarda la matriz de conteos
source("FUNC/ConstruirTablaSecuencias.R")
construirTablaSecuencias(union)
# Ruta al objeto Rds con la tabla final de secuencias
tablaSecuencias <- "OUTPUT/RDS/seqtab.Rds"

# -----------------------------------------
# ELIMINACIÓN DE QUIMERAS
# -----------------------------------------
# Detecta y remueve secuencias quiméricas
source("FUNC/EliminarQuimeras.R")
eliminarQuimeras(tablaSecuencias)
# Ruta a la tabla limpia
tablaSinQuim <- "OUTPUT/RDS/tabSinQuim.Rds"

# ------------------------------------------
# ASIGNACIÓN TAXONÓMICA
# ------------------------------------------
# Asigna taxonomóia usando bases de datos locales (HTTdb, RDP, SILVA)
# C:/ANTIGUA_D/TodoTFG/BB_DD/hitdb_v1.00.fa.gz
# C:/ANTIGUA_D/TodoTFG/BB_DD/rdp_19_toGenus_trainset.fa.gz /rdp_19_toSpecies_trainset.fa.gz
# C:/ANTIGUA_D/TodoTFG/BB_DD/silva_nr99_v138.2_toGenus_trainset.fa.gz

source("FUNC/AsignarTaxonomia.R")
# Aplica asignación taxonómica a la tabla sin quimeras
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
