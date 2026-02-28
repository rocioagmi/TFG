library(httr)     
library(jsonlite)
library(dplyr)
library(readr)


busquedaEBI <- function(query, limit = 200000) {
  
  base_url <- "https://www.ebi.ac.uk/ebisearch/ws/rest"
  
  print(paste("Consultando EBI..."))
  print(paste("Query:", query))
  
  # PRIMERO BUSCAR DOMINIO: sra-study
  
  url_study <- paste0(base_url, "/sra-study")
  resultados_study <- c()
  start <- 0
  
  repeat {
    parametros_study <- list(
      query = query,
      fields = "acc",
      format = "json",
      size = ifelse(limit > 100, 100, limit), 
      start = start)
    
    respuesta_study <- GET(url_study, query = parametros_study)
  
    if (status_code(respuesta_study) != 200) {
      cat(sprintf("\nError de conexión con EBI al consultar sra-study. Código: %d\n", status_code(respuesta_study)))
      break
    }
    
    cont_study <- httr::content(respuesta_study, "text", encoding = "UTF-8")
    json_study <- fromJSON(cont_study, flatten = TRUE)
    
    if (json_study$hitCount == 0 || length(json_study$entries) == 0){
      print(paste("No se encontraron resultados para esta consulta:", respuesta_study$url_study))
      break
    } 
    
    resultados_study <- c(json_study$entries$acc)
    
    cat(sprintf("\rEstudios recuperados: %d de %d", length(resultados_study), json_study$hitCount))
    
    if(length(resultados_study) >= limit || length(resultados_study) >= json_study$hitCount) break
    start <- start + 100
    Sys.sleep(0.2)
  }

  
  # SEGUNDO BUSCAR MUESTRA: sra-run
  
  url_run <- paste0(base_url, "/sra-run")
  resultados_run <- c()
  start <- 0
  
  repeat {
    parametros_run <- list(
      query = query,
      fields = "run_accession",
      format = "json",
      size = ifelse(limit > 100, 100, limit),
      start = start)
  
    respuesta_run <- GET(url_run, query = parametros_run)
    
    if (status_code(respuesta_run) != 200) {
      cat(sprintf("\nError de conexión con EBI al consultar sra-run. Código: %d\n", status_code(respuesta_run)))
      break
    }
    
    cont_run <- httr::content(respuesta_run, "text", encoding = "UTF-8")
    json_run <- fromJSON(cont_run, flatten = TRUE)
    
    if (json_run$hitCount == 0 || length(json_run$entries) == 0){
      print(paste("No se encontraron resultados para esta consulta:", respuesta_run$url_run))
      break
    }
    
    resultados_run <- c(resultados_run, json_run$entries$acc)
    
    cat(sprintf("\rRuns recuperados: %d de %d", length(resultados_run), json_run$hitCount))
    
    if (length(resultados_run) >= limit || length(resultados_run) >= json_run$hitCount) break
    start <- start + 100
    
    Sys.sleep(0.2)
  }
  
  df_study <- if (length(resultados_study) > 0) {
    data.frame(accession = resultados_study, tipo = "study", stringsAsFactors = FALSE)
  } else {
    data.frame(accession = character(0), tipo = character(0))
  }
  
  df_run <- if (length(resultados_run) > 0) {
    data.frame(accession = resultados_run, tipo = "run", stringsAsFactors = FALSE)
  } else {
    data.frame(accession = character(0), tipo = character(0))
  }
  
  resultado_final <- bind_rows(df_study, df_run)
  
  cat(sprintf("Total IDs recuperados: %d\n", nrow(resultado_final)))
  return(resultado_final)
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
            contenido <- httr::content(respuesta, "text", encoding = "UTF-8")
            
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
          contenido <- httr::content(respuesta, "text", encoding = "UTF-8")
          
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
      
      setTxtProgressBar(pb_run, i)
      Sys.sleep(0.3)
    }
    
    close(pb_run)
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
  

busquedaHibrida <- function(query, limit = 200000) {
  ids_ebi <- busquedaEBI(query, limit)
  
  if (is.null(ids_ebi) || nrow(ids_ebi) == 0) {
    cat("No se encontraron resultados.\n")
    return(NULL)
  }
  
  study_ids <- ids_ebi$accession[ids_ebi$tipo == "study"]
  run_ids <- ids_ebi$accession[ids_ebi$tipo == "run"]
  
  datos_completos <- continuarConENA(list(study_accessions = study_ids, run_accessions = run_ids))
  
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