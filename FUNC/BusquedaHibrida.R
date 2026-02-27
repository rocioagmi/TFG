library(httr)     
library(jsonlite)
library(dplyr)
library(readr)


busquedaEBI <- function(query, limit = 20000) {
  
  base_url <- "https://www.ebi.ac.uk/ebisearch/ws/rest"
  
  print(paste("Consultando EBI..."))
  print(paste("Query:", query))
  
  # PRIMERO BUSCAR DOMINIO: sra-study
  
  url_study <- paste0(base_url, "/sra-study")
  
  parametros_study <- list(
    query = query,
    fields = "acc",
    format = "json",
    size = 100, 
    start = 0)
  
  respuesta_study <- GET(url_study, query = parametros_study)
  
  resultados_study <- c()
  
  if (status_code(respuesta_study) == 200) {
    cont_study <- content(respuesta_study, "text", encoding = "UTF-8")
    json_study <- fromJSON(cont_study, flatten = TRUE)
    
    if (json_study$hitCount > 0 && length(json_study$entries) > 0) {
      resultados_study <- json_study$entries$acc
      cat(sprintf("Estudios encontrados: %d\n", length(resultados_study)))
    } else {
      cat("No se encontraron estudios en este dominio.\n")
    }
  } else {
    cat(sprintf("Error al consultar sra-study: código %d\n", status_code(respuesta_study)))
  }
  
  # SEGUNDO BUSCAR MUESTRA: sra-run
  
  url_run <- paste0(base_url, "sra-run")
  
  parametros_run <- list(
    query = query,
    fields = "run_accession",
    format = "json",
    size = 100,
    start = 0)
  
  respuesta_run <- GET(url_run, query = parametros_run)
  
  resultados_run <- c()
  
  if (status_code(respuesta_run) == 200) {
    cont_sample <- content(respuesta_run, "text", encoding = "UTF-8")
    json_sample <- fromJSON(cont_sample, flatten = TRUE)
    
    if (json_sample$hitCount > 0 && length(json_sample$entries) > 0) {
      resultados_run <- json_sample$entries$acc
      cat(sprintf("Muestras encontradas: %d\n", length(resultados_run)))
    } else {
      cat("No se encontraron muestras en este dominio.\n")
    }
  } else {
    cat(sprintf("Error al consultar sra-sample: código %d\n", status_code(respuesta_run)))
  }
  
  return(list(study_accessions = resultados_study,run_accessions = resultados_run))
}


continuarConENA <- function(lista, batch_size = 100) {
  cat("Obteniendo metadatos completos desde ENA\n")
  
  url <- "https://www.ebi.ac.uk/ena/portal/api/search"
  
  datos_completos <- data.frame()
  
  # PROCESAR STUDY-ACCESSIONS
  
  if (!is.null(lista$study_accessions) && length(lista$study_accessions) > 0) {
    
    cat(sprintf("\nProcesando %d estudios...\n", length(lista$study_accessions)))
    pb_studies <- txtProgressBar(min = 0, max = length(lista$study_accessions), style = 3)
    
    for (i in seq_along(lista$study_accessions)) {
      study_acc <- lista$study_accessions[i]
      
      parametros <- list(
        result = "read_run",
        query = paste0("study_accession=", study_acc),
        fields = "study_accession,sample_accession,run_accession,scientific_name,experiment_title,study_title,sample_title,run_alias,sample_alias,fastq_ftp",
        format = "tsv",
        limit = 0)
      
      tryCatch({
          respuesta <- GET(url, query = parametros, timeout(60))
          
          if (status_code(respuesta) == 200) {
            contenido <- content(respuesta, "text", encoding = "UTF-8")
            
            if (nchar(contenido) > 100) {
              df_batch <- read_tsv(contenido, show_col_types = FALSE)
              
              if (nrow(df_batch) > 0) {
                datos_completos <- bind_rows(datos_completos, df_batch)
              }
            }
          } 
        }, error = function(e) {
          warning(sprintf("Estudio al procesar el estudio %s: %s", study_acc, e$message))
      })
      
      setTxtProgressBar(pb_studies, i)
      Sys.sleep(0.3)
    }
    
    close(pb_studies)
  }
  
  # PROCESAR RUN-ACCESSIONS
  
  if (!is.null(lista$run_accessions) && length(lista$run_accessions) > 0) {
    
    cat(sprintf("\nProcesando %d muestras...\n", length(lista$run_accessions)))
    
    pb_run <- txtProgressBar(min = 0, max = length(lista$run_accessions), style = 3)
    
    for (i in seq_along(lista$run_accessions)) {
      run_acc <- lista$run_accessions[i]
      
      parametros <- list(
        result = "read_run",
        query = paste0("run_accession=", run_acc),
        fields = "study_accession,sample_accession,run_accession,scientific_name,experiment_title,study_title,sample_title,run_alias,sample_alias,fastq_ftp",
        format = "tsv",
        limit = 0)
      
      tryCatch({
        respuesta <- GET(url, query = parametros, timeout(60))
        
        if (status_code(respuesta) == 200) {
          contenido <- content(respuesta, "text", encoding = "UTF-8")
          
          if (nchar(contenido) > 100) {
            df_batch <- read_tsv(contenido, show_col_types = FALSE)
            
            if (nrow(df_batch) > 0) {
              datos_completos <- bind_rows(datos_completos, df_batch)
            }
          }
        }
      }, error = function(e) {
        warning(sprintf("Error al procesar muestra %s: %s", run_acc, e$message))
      })
      
      setTxtProgressBar(pb_samples, i)
      Sys.sleep(0.3)
    }
    
    close(pb_samples)
  }
  
  # ELIMINAR DUPLICADOS
  
  if (nrow(datos_completos) > 0) {
    n_antes <- nrow(datos_completos)
    
    datos_completos <- datos_completos %>%
      distinct(run_accession, .keep_all = TRUE)
    
    n_despues <- nrow(datos_completos)
    
    if (n_antes > n_despues) {
      cat(sprintf("\nDuplicados eliminados: %d\n", n_antes - n_despues))
    }
  }
  
  cat(sprintf("\nMuestras con metadatos (deduplicadas): %d\n", nrow(datos_completos)))
  
  return(datos_completos)
}
  

busquedaHibrida <- function(query, limit = 20000) {
  lista_ids <- busquedaEBI(query, limit)
  
  if (is.null(lista_ids) || (length(lista_ids$study_accessions) == 0 && 
                             length(lista_ids$run_accessions) == 0)) {
    cat("No se encontraron resultados.\n")
    return(NULL)
  }
  
  datos_completos <- continuarConENA(lista_ids)
  
  if (is.null(datos_completos) || nrow(datos_completos) == 0) {
    warning("No se pudieron obtener metadatos de ENA.")
    return(NULL)
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