library(httr)     
library(jsonlite)
library(dplyr)
library(readr)


busquedaEBI <- function(query, limit = 200000) {
  
  base_url <- "https://www.ebi.ac.uk/ebisearch/ws/rest"
  
  print(paste("Consultando EBI..."))
  print(paste("Query:", query))
  
  # PRIMERO BUSCAR ESTUDIOS
  
  url_study <- paste0(base_url, "/project")
  todos_estudios <- data.frame()
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
      cat(sprintf("\nError de conexión con EBI al consultar project. Código: %d\n", 
                  status_code(respuesta_study)))
      break
    }
    
    cont_study <- httr::content(respuesta_study, "text", encoding = "UTF-8")
    json_study <- fromJSON(cont_study, flatten = TRUE)
    
    if (json_study$hitCount == 0 || length(json_study$entries) == 0){
      print(paste("No se encontraron resultados para esta consulta:", 
                  respuesta_study$url_study))
      break
    } 
    
    df_study <- data.frame(accession = json_study$entries$acc, tipo = "study",
                           stringsAsFactors = FALSE)
    todos_estudios <- bind_rows(todos_estudios, df_study)
    
    if(nrow(todos_estudios) >= limit || nrow(todos_estudios) >= json_study$hitCount) break
    start <- start + 100
  }

  
  # SEGUNDO BUSCAR MUESTRAS
  
  url_run <- paste0(base_url, "/sra-run")
  todas_muestras <- data.frame()
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
      cat(sprintf("\nError de conexión con EBI al consultar sra-run. Código: %d\n",
                  status_code(respuesta_run)))
      break
    }
    
    cont_run <- httr::content(respuesta_run, "text", encoding = "UTF-8")
    json_run <- fromJSON(cont_run, flatten = TRUE)
    
    if (json_run$hitCount == 0 || length(json_run$entries) == 0){
      print(paste("No se encontraron resultados para esta consulta:", 
                  respuesta_run$url_run))
      break
    }
    
    df_run <- data.frame(accession = json_run$entries$acc, tipo = "run",
                         stringsAsFactors = FALSE)
    todas_muestras <- bind_rows(todas_muestras, df_run)
    
    if (nrow(todas_muestras) >= limit || nrow(todas_muestras) >= json_run$hitCount) break
    start <- start + 100
  }
  
  resultado_final <- bind_rows(todos_estudios, todas_muestras)
  
  cat(sprintf("\rTotal IDs recuperados: %d estudios y %d muestras.\n", nrow(todos_estudios), nrow(todas_muestras)))
  return(resultado_final)
}



continuarConENA <- function(lista, batch_size = 100) {
  cat("Obteniendo metadatos completos desde ENA\n")
  
  url <- "https://www.ebi.ac.uk/ena/portal/api/search"
  
  datos_completos <- data.frame()
  
  # PROCESAR STUDY-ACCESSIONS
  
  if (!is.null(lista$study_accessions) && length(lista$study_accessions) > 0) {
    
    cat(sprintf("\nProcesando %d estudios...\n", length(lista$study_accessions))) 
    batches_study <- ceiling(length(lista$study_accessions)/ batch_size)
    pb_study <- txtProgressBar(min = 0, max = batches_study, style = 3)
    
    for (i in seq(1, length(lista$study_accessions), by = batch_size)) {
      
      finals <- min(i + batch_size - 1, length(lista$study_accessions))
      study_acc <- lista$study_accessions[i:finals]
      
      parametros <- list(
        result = "read_run",
        query = paste0("study_accession=", study_acc, collapse = " OR "),
        fields = "study_accession,sample_accession,run_accession,scientific_name,
                  experiment_title,study_title,sample_title,run_alias,sample_alias,
                  submitted_ftp,fastq_ftp",
        format = "tsv",
        limit = 0)
      
      exito <- FALSE
      for (intento in 1:3) {
        tryCatch({
          respuesta_study <- GET(url, query = parametros, timeout(60))
          
          if (status_code(respuesta_study) == 200) {
            contenido_study <- httr::content(respuesta_study, "text", encoding = "UTF-8")
            
            if (nchar(contenido_study) > 100) {
              df_batch_study <- read_tsv(contenido_study, show_col_types = FALSE)
              
              if (nrow(df_batch_study) > 0) {
                datos_completos <- bind_rows(datos_completos, df_batch_study)
              }
            }
            exito <- TRUE
          } else {
            if (intento < 3) {
              Sys.sleep(2)
            }
          }
        }, error = function(e) {
          if (intento < 3) {
            warning(sprintf("\nError en intento %d para batch %d: %s\n",
                            intento, ceiling(i/batch_size), e$message))
            Sys.sleep(2)
          }
        })
        if (exito) break
      }
      
      if (!exito) {
        warning(sprintf("Batch %d fallido después de %d intentos - continuando con el siguiente", 
                        ceiling(i/batch_size), 3))
      }
      
      setTxtProgressBar(pb_study, ceiling(i/batch_size))
      Sys.sleep(0.3)
    }
    
    close(pb_study)
  }
  
  # PROCESAR RUN-ACCESSIONS
  
  if (!is.null(lista$run_accessions) && length(lista$run_accessions) > 0) {
    
    cat(sprintf("\nProcesando %d muestras...\n", length(lista$run_accessions)))
    batches_run <- ceiling(length(lista$run_accessions)/ batch_size)
    pb_run <- txtProgressBar(min = 0, max = batches_run, style = 3)
    
    for (i in seq(1, length(lista$run_accessions), by = batch_size)) {
      
      finalr <- min(i + batch_size -1, length(lista$run_accessions))
      run_acc <- lista$run_accessions[i:finalr]
      
      parametros <- list(
        result = "read_run",
        query = paste0("run_accession=", run_acc, collapse = " OR "),
        fields = "study_accession,sample_accession,run_accession,scientific_name,
                  experiment_title,study_title,sample_title,run_alias,sample_alias,
                  submitted_ftp,fastq_ftp",
        format = "tsv",
        limit = 0)
      
      exito <- FALSE
      for (intento in 1:3) {
        tryCatch({
          respuesta_run <- GET(url, query = parametros, timeout(60))
          
          if (status_code(respuesta_run) == 200) {
            contenido_run <- httr::content(respuesta_run, "text", encoding = "UTF-8")
            
            if (nchar(contenido_run) > 100) {
              df_batch_run <- read_tsv(contenido_run, show_col_types = FALSE)
              
              if (nrow(df_batch_run) > 0) {
                datos_completos <- bind_rows(datos_completos, df_batch_run)
              }
            }
            exito <- TRUE
          } else {
            if (intento < 3) {
              Sys.sleep(2)
            }
          }
        }, error = function(e) {
            if (intento < 3) {
              warning(sprintf("\nError en intento %d para batch %d: %s\n", 
                              intento, ceiling(i/batch_size), e$message))
              Sys.sleep(2)
            }
          })
          if (exito) break
        }
      
        if (!exito) {
          warning(sprintf("Batch %d fallido después de %d intentos - continuando con el siguiente", 
                          ceiling(i/batch_size), 3))
        }
      
        setTxtProgressBar(pb_run, ceiling(i/batch_size))
        Sys.sleep(0.3)
      }
      close(pb_run)
  }
  
  cat(sprintf("\nEn total se han encontrado %d muestras.\n", nrow(datos_completos)))
  
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
  
  cat(sprintf("\nMuestras deduplicadas: %d\n", nrow(datos_completos)))
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
  
  # ESTO LO HE USADO PARA HACER PRUEBAS, ¿LO ELIMINO?
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
  cat("-----------------------------\n\n")
  cat(sprintf("Campos CON DATOS (%d):\n", length(campos_con_datos)))
  cat("-----------------------------\n\n")
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