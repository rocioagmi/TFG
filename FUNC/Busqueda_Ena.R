library(httr)     
library(jsonlite)   
library(svDialogs)  
library(DT)         
library(utils)
library(stringr)
library(dplyr)


busquedaENA <- function(query, limit = 2000) {
  
  url <- "https://www.ebi.ac.uk/ebisearch/ws/rest/sra-run"
  fields <- "run_accession,study_accession,sample_accession,library_strategy,scientific_name,sample_alias,description,experiment_title,fastq_ftp"
  
  print(paste("Consultando EBI..."))
  print(paste("Query:", query))
  
  todas_muestras <- data.frame()
  start <- 0
  
  repeat{
    parametros <- list(
      query = query,
      fields = fields,
      format = "json",
      size = ifelse(limit > 100, 100, limit), 
      start = start
    )
    
    respuesta <- GET(url, query = parametros)
    
    if (status_code(respuesta) != 200) {
      warning("Error de conexión con EBI. Código:", status_code(respuesta))
      break
    }
    
    contenido <- httr::content(respuesta, "text", encoding = "UTF-8")
    json <- fromJSON(contenido, flatten = TRUE)
    
    if (json$hitCount == 0 || length(json$entries) == 0){
      print(paste("No se encontraron resultados para esta consulta:", respuesta$url))
      break
    } 
    
    df_temp <- as.data.frame(json$entries)
    
    colnames(df_temp) <- gsub("fields\\.", "", colnames(df_temp))
  
    todas_muestras <- bind_rows(todas_muestras, df_temp)
    
    cat(sprintf("\rMetadatos recuperados: %d de %d...", nrow(todas_muestras), json$hitCount))
    
    if(nrow(todas_muestras) >= limit || nrow(todas_muestras) >= json$hitCount) break
    start <- start + 100
  }
  
  cat("\nConsulta finalizada.\n")
  return(todas_muestras)
}


filtrarDatos <- function(df) {
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df_filtrado <- df
  
  cat("--- INICIANDO EL FILTRADO ---\n")
  
  # ===========================================================
  # FILTRAR POR ORGANISMO
  # ===========================================================
  if ("scientific_name" %in% colnames(df_filtrado)) {
    organismo <- unique(df_filtrado$scientific_name)
    
    if (length(organismo) > 1) {
      cat(sprintf("Organismos encontrados: %d\n", length(organismo)))
      
      lista_org <- dlg_list(choices = organismo, 
                            preselect = organismo,
                            multiple = TRUE,
                            title = "Selecciona organismos a MANTENER")$res
      
      if (length(lista_org) > 0) {
        df_filtrado <- df_filtrado %>% filter(scientific_name %in% lista_org)
        cat(sprintf("Organismos seleccionados: ", paste(lista_org, collapse = ","), "\n"))
      } else {
        cat("Se mantienen todos los organismos.\n")
      }
    } else {
      cat(sprintf("Solo hay un organismo: %s\n", organismo[1]))
    }
  }
  
  # ===============================================
  # FILTRAR POR ESTRATEGIA
  # ===============================================
  if ("library_strategy" %in% colnames(df_filtrado)) {
    estrategia <- unique(df_filtrado$library_strategy)
  
    if (length(estrategia) > 1) {
      
      lista_est <- dlg_list(choices = estrategia,
                            preselect = estrategia,
                            multiple = TRUE,
                            title = "Selecciona estrategias a MANTENER")$res
      
      if (length(lista_est) > 0) {
        df_filtrado <- df_filtrado %>% filter(library_strategy %in% lista_est)
        cat(sprintf("Estrategias seleccionadas: ", paste(lista_est, collapse = ","), "\n"))
      } else {
        cat("Se mantienen todos las estrategias.\n")
      }
    } else {
      cat(sprintf("Solo hay una estrategia: %s\n", estrategia[1]))
    }
  }
  
  # =====================================================
  # FILTRAR POR PALABRA CLAVE
  # =====================================================
  terminos_sugeridos <- list(
    "\\b(fecal|stool|feces)\\b" = "Muestras fecales",
    "\\b(gut|intestinal|colon|rectal)\\b" = "Muestras intestinales",
    
    "multiple\\s*sclerosis|\\bMS\\b|rrms|spms|ppms" = "Esclerosis Múltiple",
    
    "\\bcontrol\\b|\\bhealthy\\b|\\bhc\\|normobiota" = "Controles sanos",
    
    "16S\\s*(rRNA|RNA|ribosomal)?\\s*(gene|sequencing)?" = "Secuenciación 16S", 
    "metagenom(ic|e|ics)" = "Metagenómica", 
    "\\b(WGS|shotgun)\\b" = "Shotgun \ WGS"
  )
  
  metodo <- dlg_list(choices = c(
    "Usar términos sugeridos del diccionario",
    "Introducir mis propios términos",
    "Buscar en campos específicos",
    "No filtrar por palabras clave"
  ), title = "¿Cómo deseas filtrar por contenido?",
  multiple = FALSE)$res
  
  if (length(metodo) > 0) {
    
    if (metodo == "Usar términos sugeridos del diccionario") {
      
      terminos_encontrados <- list()
      
      for (patron in names(terminos_sugeridos)) {
        campos_texto <- c("library_strategy, scientific_name, sample_alias, description, experiment_title")
        campos_disponibles <- campos_texto[campos_texto %in% colnames(df_filtrado)]
        
        encontrado <- FALSE
        for (campo in campos_disponibles) {
          if (any(grepl(patron, df_filtrado[[campo]], ignore.case = TRUE, perl = TRUE))) {
            encontrado <- TRUE
            break
          }
        }
        
        if (encontrado) {
          terminos_encontrados[[terminos_sugeridos[[patron]]]] <- patron
        }
      }
      
      if (length(terminos_encontrados) > 0) {
        cat(sprintf("Términos del diccionario encontrados en los datos: %d\n", 
                    lenght(terminos_encontrados)))
        
        terminos_selec <- dlg_list(choices = names(terminos_encontrados),
                                   preselect = NULL,
                                   multiple = TRUE,
                                   title = "Selecciona términos para FILTRAR (mantener muestras que los contengan)")$res
        
        if (length(terminos_selec) > 0) {
          
          patrones_regrex <- unlist(terminos_encontrados[terminos_selec])
          
          campos_texto <- c("library_strategy, scientific_name, sample_alias, description, experiment_title")
          campos_disponibles <- campos_texto[campos_texto %in% colnames(df_filtrado)]
          
          patron_final <- paste(patrones_regrex, collapse = "|")
          
          df_filtrado <- df_filtrado %>%
            filter(if_any(all_of(campos_disponibles),
                          ~grepl(patron_final, ., ignore.case = TRUE, perl = TRUE)))
          
          cat(sprintf("Filtrado por términos personalizados: %s\n", paste(terminos_selec, collapse = ",")))
          cat(sprintf("Muestras restantes: %d\n", nrow(df_filtrado)))
        }
      } else {
        cat("No se encontraron términos del diccionario en los datos.\n")
      }
      
    } else if (metodo == "Introducir mis porpios términos") {
      
      termino_custom <- dlgInput(message = "Introduce términos separados por comas:",
                                 default = "")$res
      
      if (length(termino_custom) > 0 && nchar(termino_custom) > 0) {
        terminos_usuario <- trimws(unlist(strsplit(termino_custom, ",")))
        
        terminos_b <- paste0("\\b", terminos_usuario, "\\b")
        
        campos_texto <- c("library_strategy, scientific_name, sample_alias, description, experiment_title")
        campos_disponibles <- campos_texto[campos_texto %in% colnames(df_filtrado)]
        
        patron <- paste(terminos_b, collapse = "|")
        
        df_filtrado <- df_filtrado %>%
          filter(if_any(all_of(campos_disponibles), 
                   ~grepl(patron, ., ignore.case = TRUE, perl = TRUE)))
        
        cat(sprintf("Filtrado por términos personalizados: %s\n", 
                    paste(terminos_usuario, collapse = ", ")))
        cat(sprintf("Muestras restantes: %d\n\n", nrow(df_filtrado)))
      }
      
    } else if (metodo == "Buscar en campos específicos") {
      
      campos_disponibles <- c("run_accesion", "sample_accession", "sample_alias", "description", "experiment_title", 
                              "study_accession", "scientific_name", "library_strategy")
      campos_existentes <- campos_disponibles[campos_disponibles %in% colnames(df_filtrado)]
      
      campo_selec <- dlg_list(
        choices = campos_existentes,
        title = "Selecciona el campo donde buscar",
        multiple = FALSE)$res
      
      if (length(campo_selec) > 0) {
        termino_campo <- dlgInput(
          message = sprintf("Introduce el término a buscar en '%s':", campo_selec),
          default = "")$res
        
        if (length(termino_campo) > 0 && nchar(termino_campo) > 0) {
          df_filtrado <- df_filtrado %>%
            filter(grepl(termino_campo, .data[[campo_selec]], ignore.case = TRUE))
          
          cat(sprintf("Filtrado en campo '%s' por: %s\n", campo_selec, termino_campo))
          cat(sprintf("Muestras restantes: %d\n\n", nrow(df_filtrado))) 
        }
      }
    }  
  }
  
  # ============================================
  # RESUMEN FINAL
  # ============================================
  cat(sprintf("FILTRADO COMPLETADO\n"))
  cat(sprintf("Muestras iniciales: %d\n", nrow(df)))
  cat(sprintf("Muestras finales: %d\n", nrow(df_filtrado)))
  
  return(df_filtrado)
}


construirConsulta <- function() {
  
  query <- dlgInput(message = "Introduce el término de búsqueda:", 
                    default = "multiple sclerosis")$res
  
  if (!length(query)) return(NULL)
  
  datos_raw <- busquedaENA(query, limit = 2000)
  
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
  
  print(datatable(datos_filtrados, 
                  caption = sprintf("Muestras seleccionadas: %d", nrow(datos_filtrados)),
                  options = list(pageLength = 10, scrollX = TRUE, dom = "Bftrip", buttons = c("copy", "csv", "excel")),
                  filter = 'top',
                  selection = 'none',
                  extensions = 'Buttons')) 
  
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


explorarResultado <- function(df) {
  return(invisible(df))
}

