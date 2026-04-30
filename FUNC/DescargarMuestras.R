
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



descargarMuestras <- function(df, estudio_id){
  
  require(dplyr)
  
  # nombrar archivos descargados por ¿alias o nombre archivo?
  
  df_procesado <- df %>%
    mutate(
      submitted_ftp = na_if(submitted_ftp, ""),
      fastq_ftp = na_if(fastq_ftp, ""),
      enlaces =  coalesce(submitted_ftp, fastq_ftp)
    ) %>%
    arrange(enlaces)
  
  enlaces_validos <- df_procesado$enlaces
  vacios <- sum(is.na(enlaces_validos))
  
  if (vacios > 0) {
    cat(sprintf("Advertencia: %d muestras sin enlace de descarga disponible.\n", vacios))
    enlaces_validos <- enlaces_validos[!is.na(enlaces_validos)]
  }
  
  if (length(enlaces_validos) == 0) {
    cat("No hay enlaces válidos para descargar.\n")
    return(invisible(NULL))
  }
  
  enlaces_expandidos <- unlist(strsplit(enlaces_validos, ";"))
  enlaces_expandidos <- trimws(enlaces_expandidos)
  enlaces_expandidos <- enlaces_expandidos[nchar(enlaces_expandidos) > 0]
  
  cat(sprintf("Iniciando descarga de %d archivos\n", length(enlaces_expandidos)))
  
  fallidos <- c()
  descargados <- 0
  pb <- txtProgressBar(min = 0, max = length(enlaces_expandidos), style = 3)
  
  for (idx in seq_along(enlaces_expandidos)){
    i <- enlaces_expandidos[idx]
    nombre_archivo <- tail(strsplit(i, "/")[[1]], 1)
    
    carpeta_estudio <- file.path("INPUT", "DATA", estudio_id)
    if (!dir.exists(carpeta_estudio)) dir.create(carpeta_estudio, recursive = TRUE)
    destino <- file.path(carpeta_estudio, nombre_archivo)   
    
    exito <- FALSE
    for (intento in 1:3) {
      intento_actual <- intento
      tryCatch({
        download.file(i, destino, mode = "wb", timeout = 120)
        exito <- TRUE
      }, error = function(e) {
        if (intento_actual < 3) {
          Sys.sleep(3)
        } else {
          cat(sprintf("\nError en %s: %s\n", nombre_archivo, e$message))
          fallidos <<- c(fallidos, i)
        }
      })
      if (exito) break
    }
    if (exito) descargados <- descargados + 1
    setTxtProgressBar(pb, idx)
  }
  
  close(pb)
  
  cat(sprintf("\n--- DESCARGA COMPLETADA ---\n"))
  cat(sprintf("Archivos descargados : %d\n", descargados))
  cat(sprintf("Archivos fallidos    : %d\n", length(fallidos)))
  
  if (length(fallidos) > 0) {
    cat(sprintf("\nEnlaces fallidos:\n"))
    for (f in fallidos) cat(sprintf("  - %s\n", f))
  }
  
  return(invisible(df_procesado)) 
}
