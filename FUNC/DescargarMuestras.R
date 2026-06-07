
obtenerMuestras <- function(estudios, batch_size = 50){
  require(httr)
  require(readr)
  require(dplyr)
  
  ids <- unique(estudios$study_accession)
  cat(sprintf("Recuperando muestras de %d estudios seleccionados.\n", length(ids)))
  
  url <- "https://www.ebi.ac.uk/ena/portal/api/search"
  datos_completos <- data.frame()
  
  batches_total <- ceiling(length(ids) / batch_size)
  pb <- txtProgressBar(min = 0, max = batches_total, style = 3)
  
  for (i in seq(1, length(ids), by = batch_size)) {
    
    final     <- min(i + batch_size - 1, length(ids))
    batch_ids <- ids[i:final]
    
    parametros <- list(
      result = "read_run",
      query  = paste0("study_accession=", batch_ids, collapse = " OR "),
      fields = paste("study_accession","sample_accession","run_accession",
                "scientific_name","experiment_title","study_title","sample_title",
                "run_alias","sample_alias","submitted_ftp","fastq_ftp",sep = ","),
      format = "tsv",
      limit  = 0
    )
    
    exito <- FALSE
    
    for (intento in 1:3) {
      intento_actual <- intento
      tryCatch({
        respuesta <- GET(url, query = parametros, timeout(60))
        
        if (status_code(respuesta) == 200) {
          contenido <- httr::content(respuesta, "text", encoding = "UTF-8")
          
          if (nchar(contenido) > 100) {
            df_batch <- read_tsv(contenido, show_col_types = FALSE)
            if (nrow(df_batch) > 0) {
              datos_completos <- bind_rows(datos_completos, df_batch)
            }
          }
          exito <- TRUE
          
        } else {
          if (intento < 3) Sys.sleep(2)
        }
        
      }, error = function(e) {
        if (intento_actual < 3) {
          warning(sprintf("Error en intento %d para batch %d: %s",
                          intento_actual, ceiling(i / batch_size), e$message))
          Sys.sleep(2)
        }
      })
      if (exito) break
    }
    
    if (!exito) {
      warning(sprintf("Batch %d fallido después de 3 intentos - continuando",
                      ceiling(i / batch_size)))
    }
    
    setTxtProgressBar(pb, ceiling(i / batch_size))
    Sys.sleep(0.3)
  }
  
  close(pb)
  
  # ELIMINAR DUPLICADOS run_accession
  if (nrow(datos_completos) > 0) {
    n_antes <- nrow(datos_completos)
    datos_completos <- datos_completos %>%
      distinct(run_accession, .keep_all = TRUE)
    n_despues <- nrow(datos_completos)
    
    if (n_antes > n_despues) {
      cat(sprintf("Duplicados eliminados: %d\n", n_antes - n_despues))
    }
  }
  
  cat(sprintf("Total muestras recuperadas: %d\n", nrow(datos_completos)))
  return(datos_completos)
}

clasificarEtiqueta <- function(texto) {
  if(length(texto) == 0) return(character(0))

  palabras_ms <- c("MS", "RRMS", "SPMS", "PPMS", "Patient", "PT", "pat")
  palabras_hc <- c("HC", "HD", "Control", "Healthy", "ctrl")
  
  sep <- "[ ._0-9-]"
  patron <- function(palabras) {
    paste0("(^|", sep, ")(", paste(palabras, collapse = "|"), ")(", sep, "|$)") 
  }
  
  ms <- grepl(patron(palabras_ms), texto, ignore.case = TRUE)
  hc <- grepl(patron(palabras_hc), texto, ignore.case = TRUE)
  
  ifelse(ms & !hc, "MS",
         ifelse(hc & !ms, "HC", NA_character_))
}

detectarEtiqueta <- function(df) {
  
  if (nrow(df) <= 4) {
    idx <- seq_len(nrow(df))
  } else {
    centro <- ceiling(nrow(df) / 2)
    medio  <- floor(4 / 2)
    ini    <- max(1, centro - medio)
    fin    <- min(nrow(df), ini + 4 - 1)
    idx    <- ini:fin
  }
  
  insp <- df[idx, , drop = FALSE]
  
  contieneEtiqueta <- function(x) any(!is.na(clasificarEtiqueta(x)))
  
  nombre_archivo <- basename(sub(";.*$", "", insp$enlaces))
  
  if (contieneEtiqueta(nombre_archivo))    return("nombre")
  if (contieneEtiqueta(insp$sample_title)) return("sample_title")
  if (contieneEtiqueta(insp$sample_alias)) return("sample_alias")
  "ninguno"
}

descargarMuestras <- function(df, estudio_id){
  
  require(dplyr)
  
  df_procesado <- df %>%
    mutate(
      submitted_ftp = na_if(submitted_ftp, ""),
      fastq_ftp = na_if(fastq_ftp, ""),
      enlaces =  coalesce(submitted_ftp, fastq_ftp)
    ) %>%
    filter(!is.na(enlaces)) %>%
    arrange(enlaces)
  
  if (nrow(df_procesado) == 0) {
    cat("No hay enlaces válidos para descargar.\n")
    return(invisible(NULL))
  }
  
  carpeta_estudio <- file.path("INPUT", "DATA", estudio_id)
  if (!dir.exists(carpeta_estudio)) dir.create(carpeta_estudio, recursive = TRUE)
  
  campo <- detectarEtiqueta(df_procesado)
  renombrar <- campo %in% c("sample_title", "sample_alias")
  
  fallidos <- c()
  descargados <- 0
  contador <- 0
  
  for (idx in seq_along(df_procesado)){
    
    fila <- df_procesado[idx, ]
    
    grupo <- if (renombrar) clasificarGrupo(fila[[campo]]) else NA_character_
    
    nombre_archivo <- trimws(strsplit(fila$enlaces, ";")[[1]])
    nombre_archivo <- nombre_archivo[nchar(nombre_archivo) > 0]
    
    for (i in nombre_archivo) {
      contador <- contador + 1
      nombre_completo <- ifelse(grepl("^http|^ftp", nombre_archivo), 
                                nombre_archivo, paste0("http://", nombre_archivo))
      original <- basename(nombre_archivo) 
      
      nuevo <- if (!is.na(grupo)) paste0(grupo, "_", original) else original
      
      destino <- file.path(carpeta_estudio, nuevo)
      cat(sprintf("[%d / %d] %s\n", contador, 
                  sum(lengths(strsplit(df_procesado$enlaces, ";"))), nuevo))
      
      exito <- FALSE
      for (intento in 1:3) {
        intento_actual <- intento
        tryCatch({
          download.file(nombre_completo, destino, mode = "wb", timeout = 300, method = "libcurl")
          exito <- TRUE
        }, error = function(e) {
          if (intento_actual < 3) {
            Sys.sleep(5)
          } else {
            cat(sprintf("\nError en %s: %s\n", original, e$message))
            fallidos <<- c(fallidos, nombre_completo)
          }
        })
        if (exito) break
      }
      if (exito) descargados <- descargados + 1
    }
  }
  
  cat("--------------------------------\n")
  cat(sprintf("\n DESCARGA COMPLETADA %s\n", estudio_id))
  cat("--------------------------------\n")
  cat(sprintf("  Archivos descargados : %d\n", descargados))
  cat(sprintf("  Archivos fallidos    : %d\n", length(fallidos)))
  
  if (length(fallidos) > 0) {
    cat(sprintf("\nEnlaces fallidos:\n"))
    for (f in fallidos) cat(sprintf("  - %s\n", f))
  }
  
  return(invisible(df_procesado)) 
}
