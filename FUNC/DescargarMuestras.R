library(readr)
library(httr)
library(dplyr)
library(svDialogs)

descargarMuestras <- function(df) {
  
  if (is.null(df) || nrow(df) == 0) {
    stop("No hay muestras para descargar.")
  }
  
  modo <- dlg_list(choices = c("Descargar todas las muestras filtradas",
                               "Seleccionar estudios concretos",
                               "Seleccionar muestras concretas (run_accession)"),
                   title = "¿Qué deseas descargar?")$res
  
  if (length(modo) == 0) {
    message("Descarga cancelada.")
    return(invisible(NULL))
  }
  
  df_descarga <- switch(modo, 
                        "Descargar todas las muestras filtradas" = {
                          confirmar <- dlg_message(
                            sprintf("Vas a descargar %d muestras. ¿Continuar?", nrow(df)),
                            type = "yesno")$res
                          if (confirmar != "yes") return(invisible(NULL))
                          df_filtrado
                        },
                        
                        "Seleccionar estudios concretos" = {
                          estudios_disponibles <- unique(df$study_accession)
                          estudios_disponibles <- estudios_disponibles[!is.na(estudios_disponibles)]
                          seleccion <- dlg_list(
                            choices = estudios_disponibles,
                            title = "Selecciona uno o varios estudios (Ctrl+clic para varios)",
                            multiple = TRUE)$res
                          if (length(seleccion) == 0) return(invisible(NULL))
                          df %>% filter(study_accession %in% seleccion)
                        },
                        
                        "Seleccionar muestras concretas (run_accession)" = {
                          runs_disponibles <- unique(df$run_accession)
                          runs_disponibles <- runs_disponibles[!is.na(runs_disponibles)]
                          seleccion <- dlg_list(
                            choices = runs_disponibles,
                            title = "Selecciona las muestras (Ctrl+clic para varias)",
                            multiple = TRUE)$res
                          if (length(seleccion) == 0) return(invisible(NULL))
                          df %>% filter(run_accession %in% seleccion)
                        }
  )
  
  if (is.null(df_descarga) || nrow(df_descarga) == 0) {
    message("No hay muestras seleccionadas.")
    return(invisible(NULL))
  }
  
  descargar_muestra <- function(run_accession, enlaces_ftp) {
    if (is.na(enlaces_ftp) || enlaces_ftp == "") {
      warning(sprintf("Sin enlace FTP para %s — omitida.", run_accession))
      return(FALSE)
    }
    
    enlaces <- strsplit(enlaces_ftp, ";")[[1]]
    exito <- TRUE
    
    for (enlace in enlaces) {
      enlace <- trimws(enlace)
      if (!startsWith(enlace, "ftp://")) {
        enlace <- paste0("ftp://", enlace)
      }
      nombre_archivo <- basename(enlace)
      destino <- file.path("INPUT/DATA", nombre_archivo)
      
      if (file.exists(destino)) {
        cat(sprintf("  Ya existe: %s — omitido.\n", nombre_archivo))
        next
      }
      
      resultado <- tryCatch(
        download.file(enlace, destino, mode = "wb", quiet = TRUE),
        error = function(e) {
          warning(sprintf("Error descargando %s: %s", nombre_archivo, e$message))
          return(1)
        }
      )
      if (resultado != 0) exito <- FALSE
    }
    return(exito)
  }
  
  cat(sprintf("\nIniciando descarga de %d muestras...\n", nrow(df_descarga)))
  resultados <- vector("logical", nrow(df_descarga))
  
  for (i in seq_len(nrow(df_descarga))) {
    cat(sprintf("[%d/%d] %s\n", i, nrow(df_descarga), df_descarga$run_accession[i]))
    resultados[i] <- descargar_muestra(
      df_descarga$run_accession[i],
      df_descarga$fastq_ftp[i]      
    )
  }
  
  n_ok  <- sum(resultados)
  n_err <- sum(!resultados)
  cat(sprintf("\nDescarga completada: %d correctas, %d con errores.\n", n_ok, n_err))
  
  if (n_err > 0) {
    fallidas <- df_descarga$run_accession[!resultados]
    cat("Muestras con error:\n")
    cat(paste(" -", fallidas, collapse = "\n"), "\n")
  }
  
  return(invisible(df_descarga[resultados, ]))
}