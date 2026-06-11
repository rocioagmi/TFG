
busquedaEBI <- function(query, limit = 200000) {
  
  require(httr)     
  require(jsonlite)
  require(dplyr)
  
  base_url <- "https://www.ebi.ac.uk/ebisearch/ws/rest"
  
  print(paste("Consultando EBI..."))
  print(paste("Query:", query))
  
  dominios = c("project", "sra-study")
  
  resultado_final <- data.frame()
  conteo_dom <- list()
  
  for (dom in dominios) {
    
    url <- paste0(base_url, "/", dom)
    start <- 0
    ids_dom <- data.frame()
    
    repeat {
      pag <- min(100, limit - nrow(ids_dom))
      if (pag <= 0) break
      
      parametros <- list(
        query = query,
        fields = "acc",
        format = "json",
        size = pag, 
        start = start)
      
      respuesta <- GET(url, query = parametros)
      
      if (status_code(respuesta) != 200) {
        cuerpo <- httr::content(respuesta, "text", encoding = "UTF-8")
        cat(sprintf("\nError de conexión con EBI al consultar %s. Código: %d | %s\n", 
                    dom, status_code(respuesta), substr(cuerpo, 1, 200)))
        break
      }
      
      cont <- httr::content(respuesta, "text", encoding = "UTF-8")
      json <- fromJSON(cont, flatten = TRUE)
      
      if (is.null(json$hitCount) || json$hitCount == 0 || 
          is.null(json$entries) || length(json$entries) == 0) {
        cat(sprintf("No se encontraron resultados para %s:", dom))
        break
      } 
      
      nuevos <- data.frame(accession = json$entries$acc, 
                           stringsAsFactors = FALSE)
      ids_dom <- bind_rows(ids_dom, nuevos)
      
      if(nrow(ids_dom) >= json$hitCount || nrow(ids_dom) >= limit) break
      
      start <- start + pag
    }
    
    conteo_dom[[dom]] <- nrow(ids_dom)
    
    resultado_final <- bind_rows(resultado_final, ids_dom)
  }
  
  cat(sprintf(
    "\nTotal IDs recuperados: %d (project: %d | sra-study: %d)\n",
    nrow(resultado_final),
    conteo_dom[["project"]],
    conteo_dom[["sra-study"]]
  ))
  
  return(resultado_final)
}  
  


continuarConENA <- function(ids_ebi, batch_size = 100) {
  
  require(httr)
  require(readr)
  require(dplyr)
  
  cat("Obteniendo metadatos completos desde ENA\n")
  
  ids <- unique(ids_ebi$accession)
  
  cat(sprintf("\nProcesando %d estudios.\n", length(ids)))
  
  url <- "https://www.ebi.ac.uk/ena/portal/api/search"
  datos_completos <- data.frame()
  
  batches_total <- ceiling(length(ids)/ batch_size)
  pb <- txtProgressBar(min = 0, max = batches_total, style = 3)
    
  for (i in seq(1, length(ids), by = batch_size)) {
      
      final <- min(i + batch_size - 1, length(ids))
      batch_ids <- ids[i:final]
      
      parametros <- list(
        result = "study",
        query = paste0('study_accession="', batch_ids, '"', collapse = " OR "),
        fields = paste("study_accession","secondary_study_accession",
                       "study_title", "description","scientific_name","tax_id",
                       "center_name", sep = ","),
        format = "tsv",
        limit = 0)
      
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
            codigo <- status_code(respuesta)
            cuerpo <- httr::content(respuesta, "text", encoding = "UTF-8")
            message(sprintf("  [ENA] Batch %d intento %d -> HTTP %d | %s",
                            ceiling(i/batch_size), intento, codigo, 
                            substr(cuerpo, 1, 200)))
            if (codigo >= 400 && codigo < 500) break
            if (intento < 3) Sys.sleep(2 * intento)
          }
        }, error = function(e) {
          if (intento_actual < 3) {
            warning(sprintf("\nError en intento %d para batch %d: %s\n",
                            intento_actual, ceiling(i/batch_size), e$message))
            Sys.sleep(2)
          }
        })
        if (exito) break
      }
      
      if (!exito) {
        warning(sprintf("Batch %d fallido después de 3 intentos - continuando con el siguiente", 
                        ceiling(i/batch_size)))
      }
      
      setTxtProgressBar(pb, ceiling(i/batch_size))
      Sys.sleep(0.3)
    }
    
    close(pb)
    cat(sprintf("\nEstudios encontrados: %d\n", nrow(datos_completos)))
  
  # ELIMINAR DUPLICADOS
  if (nrow(datos_completos) > 0) {
    n_antes <- nrow(datos_completos)
    datos_completos <- datos_completos %>% 
      distinct(study_accession, .keep_all = TRUE)
    n_despues <- nrow(datos_completos)
    
    if (n_antes > n_despues) {
      cat(sprintf("\nDuplicados eliminados: %d\n", n_antes - n_despues))
    }
  }
  
  cat(sprintf("\nEstudios únicos obtenidos: %d\n", nrow(datos_completos)))
  return(datos_completos)
}
  


busquedaHibrida <- function(query, limit = 200000) {
  
  ids_ebi <- busquedaEBI(query, limit)
  
  if (is.null(ids_ebi) || nrow(ids_ebi) == 0) {
    cat("No se encontraron resultados.\n")
    return(NULL)
  }
  
  datos_completos <- continuarConENA(ids_ebi)
  
  if (is.null(datos_completos) || nrow(datos_completos) == 0) {
    warning("No se pudieron obtener metadatos de ENA.")
    return(NULL)
  }
  
  return(datos_completos)
}