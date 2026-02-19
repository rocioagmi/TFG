library(httr)     
library(jsonlite)
library(dplyr)
library(readr)
#library(svDialogs)  
#library(DT)         
#library(utils)
#library(stringr)

busquedaEBI <- function(query, limit = 20000) {
  
  url <- "https://www.ebi.ac.uk/ebisearch/ws/rest/sra-run"
  
  print(paste("Consultando EBI..."))
  print(paste("Query:", query))
  
  todas_muestras <- data.frame()
  start <- 0
  
  repeat{
    parametros <- list(
      query = query,
      fields = "run_accession",
      format = "json",
      size = ifelse(limit > 100, 100, limit), 
      start = start)
    
    respuesta <- GET(url, query = parametros)
    
    if (status_code(respuesta) != 200) {
      cat(sprintf("\nError de conexión con EBI. Código: %d\n", status_code(respuesta)))
      if (nrow(todas_muestras) > 0) {
        cat(sprintf("Se recuperaron %d registros antes del error.\n", nrow(todas_muestras)))
      }
      break
    }
    
    contenido <- httr::content(respuesta, "text", encoding = "UTF-8")
    json <- fromJSON(contenido, flatten = TRUE)
    
    if (json$hitCount == 0 || length(json$entries) == 0){
      print(paste("No se encontraron resultados para esta consulta:", respuesta$url))
      break
    } 
    
    df_temp <- data.frame( 
      run_accession = json$entries$acc,
      stringsAsFactors = FALSE)
  
    todas_muestras <- bind_rows(todas_muestras, df_temp)
    
    cat(sprintf("\rIDs recuperados: %d de %d...", nrow(todas_muestras), json$hitCount))
    
    if(nrow(todas_muestras) >= limit || nrow(todas_muestras) >= json$hitCount) break
    start <- start + 100
  }
  
  cat("\nBúsqueda en EBI finalizada.\n")
  return(todas_muestras)
}


continuarConENA <- function(run_accessions, batch_size = 100) {
  cat("Obteniendo metadatos completos desde ENA\n")
  
  url <- "https://www.ebi.ac.uk/ena/portal/api/search"
  
  datos_completos <- data.frame()
  total_batches <- ceiling(length(run_accessions) / batch_size)
  
  pb <- txtProgressBar(min = 0, max = total_batches, style = 3)
  
  for (i in seq(1, length(run_accessions), by = batch_size)) {
    
    end_idx <- min(i + batch_size - 1, length(run_accessions))
    batch <- run_accessions[i:end_idx]
    
    query_batch <- paste0("run_accession=", batch, collapse = " OR ")
    
    parametros <- list(
      result = "read_run",
      query = query_batch,
      fields = "study_accession,sample_accession,run_accession,scientific_name,experiment_title,study_title,sample_title,run_alias,sample_alias,fastq_ftp",
      format = "tsv",
      limit = 0)
    
    exito <- FALSE
    max_intentos <- 3
    
    for (intento in 1:max_intentos) {
      tryCatch({
        respuesta <- GET(url, query = parametros, timeout(60))
        
        if (status_code(respuesta) == 200) {
          contenido <- httr::content(respuesta, "text", encoding = "UTF-8")
          
          if (nchar(contenido) > 100) {
            df_batch <- read_tsv(contenido, show_col_types = FALSE)
            
            if (nrow(df_batch) > 0) {
              if("tax_id" %in% colnames(df_batch)){
                df_batch$tax_id <- as.character(df_batch$tax_id)
              }
              
              datos_completos <- bind_rows(datos_completos, df_batch)
            }
          }
          exito <- TRUE
          break  
        } else {
          
          if (intento < max_intentos) {
            Sys.sleep(2)  
          }
        }
        
      }, error = function(e) {
        
        if (intento < max_intentos) {
          cat(sprintf("\nError en intento %d para batch %d: %s\n", 
                      intento, ceiling(i/batch_size), e$message))
          Sys.sleep(2)
        }
      })
      
      if (exito) break  
    }
    
    if (!exito) {
      warning(sprintf("Batch %d fallido después de %d intentos - continuando con el siguiente", 
                      ceiling(i/batch_size), max_intentos))
    }
    
    setTxtProgressBar(pb, ceiling(i/batch_size))
    Sys.sleep(0.3)
  }
  
  close(pb)
  
  cat(sprintf("\nMuestras con metadatos: %d\n", nrow(datos_completos)))
  return(datos_completos)
}


busquedaHibrida <- function(query, limit = 20000) {
  ids_ebi <- busquedaEBI(query, limit)
  
  if (is.null(ids_ebi) || nrow(ids_ebi) == 0) {
    return(NULL)
  }
  
  datos_completos <- continuarConENA(ids_ebi$run_accession)
  
  if (is.null(datos_completos) || nrow(datos_completos) == 0) {
    warning("No se pudieron obtener metadatos de ENA. Usando solo IDs de EBI.")
    return(ids_ebi)
  }
  
  campos_con_datos <- c()
  campos_vacios <- c()
  
  for (col in names(datos_completos)) {
    valores_validos <- datos_completos[[col]][!is.na(datos_completos[[col]]) & 
                                                datos_completos[[col]] != ""]
    if (length(valores_validos) > 0) {
      porcentaje <- (length(valores_validos) / nrow(datos_completos)) * 100
      campos_con_datos <- c(campos_con_datos, 
                            sprintf("%s (%.1f%%)", col, porcentaje))
    } else {
      campos_vacios <- c(campos_vacios, col)
    }
  }
  
  cat(sprintf("Campos CON DATOS (%d):\n", length(campos_con_datos)))
  for (campo in campos_con_datos) {
    cat(sprintf("  - %s\n", campo))
  }
  
  if (length(campos_vacios) > 0) {
    cat(sprintf("\n Campos VACÍOS (%d): %s\n", 
                length(campos_vacios), 
                paste(campos_vacios, collapse = ", ")))
  }
  cat("-----------------------------\n\n")
  
  return(datos_completos)
}