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
library(rmarkdown)


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

# ------------------------------------------------
# FILTRADO DE RESULTADOS DE BÚSQUEDA
# ------------------------------------------------
respuesta_filtro <- dlg_message( 
  message = sprintf("Se encontraron %d estudios. ¿Deseas filtrarlos?", nrow(muestras_totales)),
  type = "yesno")$res

if (respuesta_filtro == "yes") {
  busqueda_filtrada <- filtrarBusqueda(muestras_totales)
} else {
  busqueda_filtrada <- muestras_totales
  cat("Usando todas los estudios sin filtrar.\n")
}

if (is.null(busqueda_filtrada) || nrow(busqueda_filtrada) == 0) {
  message("No quedan estudios después del filtrado.")
  stop()
}

# ---------------------------------------------------------------------------
# MUESTRA TABLA CON LOS RESULTADOS DE BÚSQUEDA FILTRADOS PARA SU SELECCIÓN
# ---------------------------------------------------------------------------
cat("Abriendo aplicación Shiny para explorar y descargar resultados.\n")
source("FUNC/mostrarTabla.R")

df_final <- mostrarTabla(busqueda_filtrada)


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

muestras <- obtenerMuestras(df_final)

estudios_ids <- unique(muestras$study_accession)

for (estudio in estudios_ids) {
  df_estudio <- muestras %>% filter(study_accession == estudio)
  descargarMuestras(df_estudio, estudio_id = estudio)
}


# ================================================
# PREPROCESAMIENTO DE LOS DATOS
# ================================================

# Obtiene una lista ordenada con todos los estudios y sus muestras .fastq.gz

estudio_ids <- basename(list.dirs("INPUT/DATA", recursive = FALSE))
estudios <- list()

for (estudio in estudio_ids) {
  carpeta_estudio <- file.path("INPUT/DATA", estudio)
  muestras_estudio <- sort(list.files(carpeta_estudio, 
                                      pattern = "\\.fastq\\.gz$",
                                      full.names = TRUE))
  if(any(grepl("_R1", muestras_estudio))){
    muestra_R1 <- muestras_estudio[grepl("_R1", muestras_estudio)]
    muestra_R2 <- muestras_estudio[grepl("_R2", muestras_estudio)]
  } else if (any(grepl("_1\\.fastq", muestras_estudio))){
    muestra_R1 <- muestras_estudio[grepl("_1\\.fastq", muestras_estudio)]
    muestra_R2 <- muestras_estudio[grepl("_2\\.fastq", muestras_estudio)]
  } else {
    stop("No se reconoce el patrón de nomenclatura de los archivos.")
  }
  
  estudios[[estudio]] <- list(R1 = muestra_R1, R2 = muestra_R2)
}

# -------------------------------------------------
# INFORME DE CALIDAD DE SECUENCIACIÓN              
# -------------------------------------------------
# Carga la función para generar informes FastQC
source("FUNC/InformeCalidad.R")

# Genera informe de calidad para cada estudio
estudio <- names(estudios)[1]
muestras <- estudios[[estudio]]

datos_qa <- informeCalidad(estudio_id = estudio, 
               muestras = muestras,
               etapa = "prefiltrado",
               umbral_calidad = 20,
               umbral_ee = 2)

# ------------------------------------------------
# FILTRADO DE CALIDAD CON DADA2                   
# ------------------------------------------------
# Carga función de filtrado 
source("FUNC/FiltrarMuestras.R")

propuesta <- proponerParametros(datos_qa, umbral_calidad = 20)
parametros <- confirmarParametros(propuesta, interactivo = TRUE)
filtrarMuestras(muestras, parametros, estudio)


# Lista todos los archivos filtrados generados para el estudio actual

dir_filtrado <- file.path("OUTPUT", estudio, "FILTRADO")
listadoFiltrado <- sort(list.files(dir_filtrado, pattern = "\\.fastq\\.gz$",
                                   full.names = TRUE))

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
muestras_filtradas <- list(R1 = fR1, R2 = fR2)

datos_qa_post <- informeCalidad(estudio_id = estudio,
                                muestras = muestras_filtradas,
                                etapa = "postfiltrado",
                                umbral_calidad = 20,
                                umbral_ee = 2)


# ===========================================
# PROCESAMIENTO DE LOS DATOS 
# ===========================================

# -------------------------------------------
# FLUJO PRINCIPAL DADA2
# -------------------------------------------
# Carga la función que ejecuta: learnErrors, derep, dada, mergePairs
source("FUNC/FlujoTrabajoDada.R")
# Procesa lecturas filtradas y devuelve la tabla de conteos de secuencias
union <- flujoTrabajoDada(muestras_filtradas$R1, muestras_filtradas$R2, estudio)


# --------------------------------------------------------------------
# CONSTRUCCIÓN DE LA TABLA DE SECUENCIAS Y ELIMINACIÓN DE QUIMERAS
# --------------------------------------------------------------------
# Construye la matriz de conteos de ASVs y elimina secuencias quiméricas
source("FUNC/procesarSecuencias.R")
tabSinQuim <- procesarSecuencias(union, estudio)


# ------------------------------------------
# ASIGNACIÓN TAXONÓMICA
# ------------------------------------------
source("FUNC/AsignarTaxonomia.R")
# Aplica asignación taxonómica a la tabla sin quimeras
asignarTaxonomia(tabSinQuim)


# -----------------------------------------
# INFORME FINAL
# -----------------------------------------
source("FUNC/InformeFinal.R")
informeFinal()