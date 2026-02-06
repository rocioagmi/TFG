library(httr)     
library(jsonlite)   
library(svDialogs)  
library(DT)         
library(utils)
library(stringr)
library(dplyr)
library(readr)

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


filtrarDatos <- function(df) {
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df_filtrado <- df
  
  cat("INICIANDO FILTRADO\n")
  cat(sprintf("Total de muestras inicial: %d\n\n", nrow(df_filtrado)))
  
  campos_texto <- c("scientific_name", "experiment_title", "study_title",
                    "sample_title", "run_alias", "sample_alias")
  campos_disponibles <- campos_texto[campos_texto %in% colnames(df_filtrado)]
  
  # Solo los que tengan al menos un valor no nulo y no vacío
  campos_con_datos <- campos_disponibles[sapply(campos_disponibles, function(col) {
    any(!is.na(df_filtrado[[col]]) & df_filtrado[[col]] != "")
  })]
  
  if (length(campos_con_datos) == 0) {
    cat("No hay campos de texto con datos para filtrar.\n\n")
    return(df_filtrado)
  }
  
  cat(sprintf("Campos disponibles para buscar: %s\n\n",
              paste(campos_con_datos, collapse = ", ")))
  
  metodo <- dlg_list(
    choices = c(
      "Usar términos sugeridos del diccionario",
      "Introducir mis propios términos",
      "Buscar en campos específicos"),
    title = "¿Cómo deseas filtrar?",
    multiple = FALSE)$res
  
  if (length(metodo) == 0) {
    cat("No se seleccionó método. Se devuelven las muestras sin filtrar.\n")
    return(df_filtrado)
  }
  
  # TÉRMINOS SUGERIDOS DEL DICCIONARIO - TERMINAR DE MIRAR ESTO
  # filtro automático que garantiza muestras humanas, intestinales o fecales de microbioma y/o 16S

  if (metodo == "Usar términos sugeridos del diccionario") {
    
    # Cada elemento de la lista es un grupo de patrones que se combinen con OR entre sí.
    # Los GRUPOS se combinan con AND (la muestra debe pasar todos los grupos).
    filtros_obligatorios <- list(
      "Humano" = "\\b(human|homo\\s*sapiens)\\b",
      "Intestinal/Fecal" = "\\b(fecal|stool|feces|faecal|gut|intestin(al|e)|colon|rectal|bowel)\\b")
    #"Microbioma/Metagenómica" = "\\b(microbiom(e|a)|metagenomic(s|a)?|16S)\\b"
    
    # Función auxiliar: aplica un patrón regex contra todos los campos de texto
    # y devuelve un vector lógico de largo nrow(df)
    aplica_patron <- function(df_local, patron, campos) {
      matches <- rep(FALSE, nrow(df_local))
      for (col in campos) {
        valores <- df_local[[col]]
        valores[is.na(valores)] <- ""
        matches <- matches | grepl(patron, valores, ignore.case = TRUE, perl = TRUE)
      }
      matches
    }
    
    cat("\n--- Filtro automático: microbioma intestinal humano ---\n")
    
    for (nombre_grupo in names(filtros_obligatorios)) {
      patron <- filtros_obligatorios[[nombre_grupo]]
      mascara <- aplica_patron(df_filtrado, patron, campos_con_datos)
      
      n_antes  <- nrow(df_filtrado)
      df_filtrado <- df_filtrado[mascara, , drop = FALSE]
      n_despues <- nrow(df_filtrado)
      
      cat(sprintf("  [%s] %d -> %d muestras (%.1f%% superaron)\n",
                  nombre_grupo, n_antes, n_despues,
                  ifelse(n_antes > 0, (n_despues / n_antes) * 100, 0)))
      
      if (nrow(df_filtrado) == 0) {
        cat("\n  ⚠ No quedan muestras tras aplicar el filtro obligatorio.\n")
        break
      }
    }
    
    # Si tras los tres filtros obligatorios queda nada, se relaja el filtro
    # "Microbioma/Metagenómica" (el más restrictivo normalmente) y se avisa.
    if (nrow(df_filtrado) == 0) {
      cat("\n  → Relajando criterio 'Microbioma/Metagenómica'...\n")
      df_filtrado <- df  # volver a principio
      
      filtros_relajados <- filtros_obligatorios[c("Humano", "Intestinal/Fecal")]
      for (nombre_grupo in names(filtros_relajados)) {
        patron <- filtros_relajados[[nombre_grupo]]
        mascara <- aplica_patron(df_filtrado, patron, campos_con_datos)
        df_filtrado <- df_filtrado[mascara, , drop = FALSE]
        
        cat(sprintf("  [%s (relajado)] -> %d muestras\n",
                    nombre_grupo, nrow(df_filtrado)))
        
        if (nrow(df_filtrado) == 0) break
      }
    }
    
    # INTRODUCIR MIS PROPIOS TÉRMINOS
    # El usuario escribe uno o varios términos separados por comas y se buscan en TODOS los campos con datos.
  
  } else if (metodo == "Introducir mis propios términos") {
    
    termino_custom <- dlgInput(
      message = "Introduce términos separados por comas\n(se combinan con OR: la muestra debe contener al menos uno):",
      default = ""
    )$res
    
    if (length(termino_custom) == 0 || nchar(trimws(termino_custom)) == 0) {
      cat("No se introdujeron términos. Se devuelven las muestras sin filtrar.\n")
      return(df_filtrado)
    }
    
    terminos_usuario <- trimws(unlist(strsplit(termino_custom, ",")))
    terminos_usuario <- terminos_usuario[nchar(terminos_usuario) > 0]  # eliminar vacíos
    
    # Cada término se encierra en \b para buscar palabras completas
    patron <- paste0("\\b(", paste(terminos_usuario, collapse = "|"), ")\\b")
    
    cat(sprintf("\nBuscando términos: %s\n", paste(terminos_usuario, collapse = ", ")))
    cat(sprintf("Campos donde se busca: %s\n", paste(campos_con_datos, collapse = ", ")))
    
    df_filtrado <- df_filtrado %>%
      filter(if_any(all_of(campos_con_datos),
                    ~grepl(patron, ., ignore.case = TRUE, perl = TRUE)))
    
    cat(sprintf("Muestras tras filtro: %d\n", nrow(df_filtrado)))
    
    # BUSCAR EN CAMPOS ESPECÍFICOS
    # El usuario selecciona uno o varios campos y después introduce uno o varios términos (separados por comas).

  } else if (metodo == "Buscar en campos específicos") {
    
    campos_selec <- dlg_list(
      choices = campos_con_datos,
      title = "Selecciona los campos donde buscar (puedes elegir más de uno)",
      multiple = TRUE
    )$res
    
    if (length(campos_selec) == 0) {
      cat("No se seleccionaron campos. Se devuelven las muestras sin filtrar.\n")
      return(df_filtrado)
    }
    
    termino_campo <- dlgInput(
      message = sprintf("Introduce términos separados por comas\n(se buscarán en: %s):",
                        paste(campos_selec, collapse = ", ")),
      default = ""
    )$res
    
    if (length(termino_campo) == 0 || nchar(trimws(termino_campo)) == 0) {
      cat("No se introdujeron términos. Se devuelven las muestras sin filtrar.\n")
      return(df_filtrado)
    }
    
    terminos <- trimws(unlist(strsplit(termino_campo, ",")))
    terminos <- terminos[nchar(terminos) > 0]
    
    patron <- paste0("\\b(", paste(terminos, collapse = "|"), ")\\b")
    
    cat(sprintf("\nBuscando términos: %s\n", paste(terminos, collapse = ", ")))
    cat(sprintf("En campos: %s\n", paste(campos_selec, collapse = ", ")))
    
    df_filtrado <- df_filtrado %>%
      filter(if_any(all_of(campos_selec),
                    ~grepl(patron, ., ignore.case = TRUE, perl = TRUE)))
    
    cat(sprintf("Muestras tras filtro: %d\n", nrow(df_filtrado)))
  }
  
  # ===========================================================
  # RESUMEN FINAL
  # ===========================================================
  cat(sprintf("\n=== FILTRADO COMPLETADO ===\n"))
  cat(sprintf("Muestras iniciales : %d\n", nrow(df)))
  cat(sprintf("Muestras finales   : %d\n", nrow(df_filtrado)))
  cat(sprintf("Muestras eliminadas: %d (%.1f%%)\n",
              nrow(df) - nrow(df_filtrado),
              ifelse(nrow(df) > 0, ((nrow(df) - nrow(df_filtrado)) / nrow(df)) * 100, 0)))
  
  return(df_filtrado)
}


construirConsulta <- function() {
  
  query <- "(multiple sclerosis) OR MS OR RRMS OR SPMS OR PPMS OR (relapsing remitting) OR (primary progressive) OR (secondary progressive)"
  
  datos_raw <- busquedaHibrida(query, limit = 20000)
  
  if (is.null(datos_raw) || nrow(datos_raw) == 0) {
    message("No se encontraron resultados para esta búsqueda.")
    return(NULL)
  }
  
  respuesta_filtro <- dlg_message( message = sprintf("Se encontraron %d muestras. ¿Deseas filtrarlas?", nrow(datos_raw)),
                                   type = "yesno")$res
  
  if (respuesta_filtro == "yes") {
    datos_filtrados <- filtrarDatos(datos_raw)
  } else {
    datos_filtrados <- datos_raw
    cat("Usando todas las muestras sin filtrar.\n")
  }
  
  if (is.null(datos_filtrados) || nrow(datos_filtrados) == 0) {
    message("No quedan muestras después del filtrado.")
    return(NULL)
  }
  
  suppressWarnings(
    print(datatable(datos_filtrados, 
                    caption = sprintf("Muestras seleccionadas: %d", nrow(datos_filtrados)),
                    options = list(pageLength = 10, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv", "excel"), deferRender = TRUE, scrollY = 400, scroller = TRUE),
                    filter = 'top',
                    selection = 'none',
                    extensions = c('Buttons', 'Scroller')))
    )
    
  guardar <- dlg_message(message = "¿Deseas guardar estos resultados en un archivo CSV?",
                         type = "yesno")$res
 
  if (guardar == "yes") {
    nombre_archivo <- sprintf("muestras_ENA_%s_%s.csv", 
                              gsub(" ", "_", query),
                              format(Sys.Date(), "%Y%m%d"))
    write_csv(datos_filtrados, nombre_archivo)
    cat(sprintf("\n Resultados guardados en: %s\n", nombre_archivo))
  }
  
  return(datos_filtrados)
}
