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

#library(dplyr)
#library(Biostrings)
#library(mongolite)
#library(TAF)
#library(QuasR)
#library(Rqc)
#library(BiocParallel)
#library(ggplot2)
library(httr)
library(jsonlite)
library(svDialogs)
library(DT)
library(readr)
library(ShortRead)
library(dada2)



# ================================================
# BÚSQUEDA PROGRAMÁTICA (EBI/ENA)                 
# ================================================
source("FUNC/BusquedaHibrida.R")
source("FUNC/FiltradoBusqueda.R")

query <- "(multiple sclerosis) OR MS OR RRMS OR SPMS OR PPMS OR (relapsing remitting) OR (primary progressive) OR (secondary progressive)"
muestras_totales <- busquedaHibrida(query, limit = 20000)

if (is.null(muestras_totales) || nrow(muestras_totales) == 0) {
  message("No se encontraron resultados para esta búsqueda.")
  return(NULL)
}

respuesta_filtro <- dlg_message( 
  message = sprintf("Se encontraron %d muestras. ¿Deseas filtrarlas?", nrow(muestras_totales)),
  type = "yesno")$res

if (respuesta_filtro == "yes") {
  busqueda_filtrada <- filtrarBusqueda(muestras_totales)
} else {
  busqueda_filtrada <- muestras_totales
  cat("Usando todas las muestras sin filtrar.\n")
}

if (is.null(busqueda_filtrada) || nrow(busqueda_filtrada) == 0) {
  message("No quedan muestras después del filtrado.")
  return(NULL)
}

suppressWarnings(
  print(datatable(busqueda_filtrada, 
                  caption = sprintf("Muestras seleccionadas: %d", nrow(busqueda_filtrada)),
                  options = list(pageLength = 10, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv", "excel"), deferRender = TRUE, scrollY = 400, scroller = TRUE),
                  filter = 'top',
                  selection = 'none',
                  extensions = c('Buttons', 'Scroller')))
)

guardar <- dlg_message(message = "¿Deseas guardar estos resultados en un archivo CSV?",
                       type = "yesno")$res

if (guardar == "yes") {
  nombre_archivo <- sprintf("busqueda_ENA_%s_%s.csv", 
                            gsub(" ", "_", query),
                            format(Sys.Date(), "%Y%m%d"))
  write_csv(busqueda_filtrada, nombre_archivo)
  cat(sprintf("\n Resultados guardados en: %s\n", nombre_archivo))
}



# ================================================
# DESCARGA DE LOS DATOS - RETOCAR ESTA FUNCIÓN    <-----
# ================================================
source("FUNC/Descargas_ENA.R")
nAcceso <- dlgInput(message = "Introduzca el número de acceso al proyecto ENA: ")$res
descargas_ENA(nAcceso)

# Obtiene una lista ordenada con todos los archivos .fastq.gz de la carpeta de entrada
listadoMuestras <- sort(list.files("INPUT/DATA", pattern = "\\.fastq\\.gz$", full.names = TRUE))

if(any(grepl("_R1", listadoMuestras))){
  R1 <- listadoMuestras[grepl("_R1", listadoMuestras)]
  R2 <- listadoMuestras[grepl("_R2", listadoMuestras)]
} else if (any(grepl("_1\\.fastq", listadoMuestras))){
  R1 <- listadoMuestras[grepl("_1\\.fastq", listadoMuestras)]
  R2 <- listadoMuestras[grepl("_2\\.fastq", listadoMuestras)]
} else {
  stop("No se reconoce el patrón de nomenclatura de los archivos.")
}

# Separa las muestras por tipo: MS (enfermos) y Healthy (sanos) ---- TODAVIA NO LO ESTOY USANDO
#muestrasMS <- sort(list.files("INPUT/DATA", pattern = "MS", full.names = TRUE))
#muestrasHealthy <- sort(list.files("INPUT/DATA", pattern = "Healthy", full.names = TRUE))

# Separa los archivos por lectura: R1 (forward) y R2 (reverse)
#muestrasR1 <- sort(list.files("INPUT/DATA", pattern = "R1", full.names = TRUE))
#muestrasR2 <- sort(list.files("INPUT/DATA", pattern = "R2", full.names = TRUE))



# ================================================
# PREPROCESAMIENTO DE LOS DATOS
# ================================================

# -------------------------------------------------
# CREACIÓN DE CARPETAS DE SALIDA
# -------------------------------------------------
# Crea directorios para organizar los resultados: informes, gráficos, datos filtrados, ...
dir.create("OUTPUT/REPORT")
dir.create("OUTPUT/FIGURES")
dir.create("OUTPUT/RDS")
dir.create("OUTPUT/FILTRADO")

# -------------------------------------------------
# INFORME DE CALIDAD DE SECUENCIACIÓN              <-----
# -------------------------------------------------
# Carga la función para generar informes FastQC
source("FUNC/InformeCalidad.R")
# Genera informe de calidad para todas las muestras
informeCalidad(listadoMuestras, umbral_calidad = 20)

# ------------------------------------------------
# FILTRADO DE CALIDAD CON DADA2                   <-----
# ------------------------------------------------
# Carga función de filtrado 
source("FUNC/FiltrarMuestras.R")
filtrarMuestras(R1, R2)

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
informeCalidad(listadoFiltrado, umbral_calidad = 20)


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
