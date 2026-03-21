library(rmarkdown)
library(svDialogs)

generarInformeFinal <- function() {
  
  cat("\nGenerando Informe Final del Pipeline...\n")
  
  # Definir rutas
  dir_salida <- "OUTPUT/REPORT"
  dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)
  archivo_rmd <- file.path(dir_salida, "InformeFinal.Rmd")
  
  # Buscar el archivo de taxonomía (el que no es seqtab ni tabSinQuim)
  archivos_rds <- list.files("OUTPUT/RDS", pattern = "\\.Rds$", full.names = TRUE)
  archivo_taxa <- archivos_rds[!grepl("seqtab\\.Rds|tabSinQuim\\.Rds", archivos_rds)]
  
  # Si hay varios, cogemos el primero o dejamos vacío si no hay
  ruta_taxa <- if(length(archivo_taxa) > 0) archivo_taxa[1] else "NA"
  
  # Escribir el contenido del RMarkdown
  contenido_rmd <- c(
    '---',
    'title: "Resumen del Procesamiento Metagenómico (TFG)"',
    'author: "Rocío"',
    paste0('date: "', format(Sys.time(), "%Y-%m-%d"), '"'),
    'output:',
    '  html_document:',
    '    theme: flatly',
    '    toc: true',
    '    toc_float: true',
    '---',
    '',
    '```{r setup, include=FALSE}',
    'knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)',
    'library(DT)',
    'library(dplyr)',
    'library(ggplot2)',
    '```',
    '',
    '## 1. Introducción',
    'Este documento resume los resultados finales del pipeline de procesamiento de datos de secuenciación (metadatos recuperados de ENA, procesados con DADA2).',
    '',
    '## 2. Construcción de la Tabla de Secuencias (ASVs)',
    '',
    '```{r carga_datos}',
    '# Cargar objetos guardados en RDS',
    'ruta_seqtab <- "../../RDS/seqtab.Rds"',
    'ruta_quim <- "../../RDS/tabSinQuim.Rds"',
    '',
    'seqtab <- if(file.exists(ruta_seqtab)) readRDS(ruta_seqtab) else NULL',
    'tabSinQuim <- if(file.exists(ruta_quim)) readRDS(ruta_quim) else NULL',
    '```',
    '',
    '```{r resumen_secuencias, results="asis"}',
    'if(!is.null(seqtab) & !is.null(tabSinQuim)){',
    '  cat(sprintf("- **Muestras procesadas:** %d\\n", nrow(seqtab)))',
    '  cat(sprintf("- **ASVs identificadas inicialmente:** %d\\n", ncol(seqtab)))',
    '  cat(sprintf("- **ASVs tras eliminar quimeras:** %d\\n", ncol(tabSinQuim)))',
    '  ',
    '  retenidas <- sum(tabSinQuim)/sum(seqtab) * 100',
    '  cat(sprintf("- **Porcentaje de lecturas no quiméricas retenidas:** %.2f%%\\n", retenidas))',
    '} else {',
    '  cat("No se encontraron los archivos RDS de secuencias.")',
    '}',
    '```',
    '',
    '## 3. Distribución de Longitud de Lecturas',
    '',
    '```{r grafico_longitud, fig.height=4, fig.width=8}',
    'if(!is.null(seqtab)){',
    '  longitudes <- nchar(getSequences(seqtab))',
    '  df_long <- data.frame(Longitud = longitudes)',
    '  ',
    '  ggplot(df_long, aes(x = Longitud)) +',
    '    geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +',
    '    labs(title = "Distribución de longitudes de las ASVs", x = "Longitud (pb)", y = "Frecuencia") +',
    '    theme_minimal()',
    '}',
    '```',
    '',
    '## 4. Asignación Taxonómica',
    '',
    '```{r taxonomia}',
    paste0('ruta_taxa <- "', ruta_taxa, '"'),
    'if(ruta_taxa != "NA" && file.exists(ruta_taxa)){',
    '  taxa <- readRDS(ruta_taxa)',
    '  # Mostramos solo los primeros 15 registros para no saturar el HTML',
    '  datatable(head(taxa, 15), options = list(scrollX = TRUE, pageLength = 5), caption = "Vista previa de las primeras 15 asignaciones taxonómicas")',
    '} else {',
    '  cat("No se encontró archivo de taxonomía generado o la ruta es inválida.")',
    '}',
    '```'
  )
  
  writeLines(contenido_rmd, archivo_rmd)
  
  # Renderizar a HTML
  cat("Compilando HTML...\n")
  rmarkdown::render(archivo_rmd, quiet = TRUE, output_file = "InformeFinal.html")
  
  ruta_final <- file.path(dir_salida, "InformeFinal.html")
  browseURL(ruta_final)
  
  cat(sprintf("\n¡Informe generado con éxito en: %s!\n", ruta_final))
}