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

# shiny
#install.packages("shiny")

# miniUI
#install.packages("miniUI")
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



# ================================================
# BÚSQUEDA PROGRAMÁTICA (EBI/ENA)                 
# ================================================
source("FUNC/BusquedaHibrida.R")
source("FUNC/FiltradoBusqueda.R")

query <- "(multiple sclerosis) OR MS OR RRMS OR SPMS OR PPMS OR (relapsing remitting) OR (primary progressive) OR (secondary progressive)"
muestras_totales <- busquedaHibrida(query, limit = 200000)

if (is.null(muestras_totales) || nrow(muestras_totales) == 0) {
  message("No se encontraron resultados para esta búsqueda.")
  stop()
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
  stop()
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
# DESCARGA DE LOS DATOS    
# ================================================

# ------------------------------------------------
# CREAR CARPETA PARA DATOS DE ENTRADA
# ------------------------------------------------
dir.create("INPUT/DATA", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# FUNCIÓN DE DESCARGA DE MUESTRAS
# ------------------------------------------------
source("FUNC/DescargarMuestras.R")

descarga <- dlgInput(message = "Introduce los IDs de los estudios que deseas descargar separados por comas:",
                     default = "")$res

if (length(descarga) == 0 || nchar(trimws(descarga)) == 0) {
  cat("No se ha descargado nada.\n")
  return(busqueda_filtrada)
}

descargarMuestras(busqueda_filtrada, descarga)


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



# ================================================
# PREPROCESAMIENTO DE LOS DATOS
# ================================================

# -------------------------------------------------
# CREACIÓN DE CARPETAS DE SALIDA
# -------------------------------------------------
# Crea directorios para organizar los resultados: informes, gráficos, datos filtrados, ...
dir.create("OUTPUT/REPORT", recursive = TRUE, showWarnings = FALSE)
dir.create("OUTPUT/RDS", recursive = TRUE, showWarnings = FALSE)
dir.create("OUTPUT/FILTRADO", recursive = TRUE, showWarnings = FALSE)

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

if(any(grepl("_R1", listadoFiltrado))){
  fR1 <- listadoFiltrado[grepl("_R1", listadoFiltrado)]
  fR2 <- listadoFiltrado[grepl("_R2", listadoFiltrado)]
} else if (any(grepl("_1\\.fastq", listadoFiltrado))){
  fR1 <- listadoFiltrado[grepl("_1\\.fastq", listadoFiltrado)]
  fR2 <- listadoFiltrado[grepl("_2\\.fastq", listadoFiltrado)]
} else {
  stop("No se reconoce el patrón de nomenclatura de los archivos.")
}


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
union <- flujoTrabajoDada(fR1, fR2)


# --------------------------------------------------------------------
# CONSTRUCCIÓN DE LA TABLA DE SECUENCIAS Y ELIMINACIÓN DE QUIMERAS
# --------------------------------------------------------------------
# Construye la matriz de conteos de ASVs y elimina secuencias quiméricas
source("FUNC/procesarSecuencias.R")
tabSinQuim <- procesarSecuencias(union)


# ------------------------------------------
# ASIGNACIÓN TAXONÓMICA
# ------------------------------------------
source("FUNC/AsignarTaxonomia.R")
# Aplica asignación taxonómica a la tabla sin quimeras
asignarTaxonomia(tabSinQuim)


# -----------------------------------------
# INFORME FINAL
# -----------------------------------------
source("FUNC/GenerarInformeFinal.R")
generarInformeFinal()