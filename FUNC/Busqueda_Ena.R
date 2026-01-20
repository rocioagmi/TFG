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
  
  # FILTRAR POR ORGANISMO
  organismo <- unique(df_filtrado$scientific_name)
  
  if(length(organismo) > 0){
    lista_org <- dlg_list(choices = organismo, 
                      preselect = organismo,
                      multiple = TRUE,
                      title = "Selecciona organismos a MANTENER")$res
    if(length(lista_org) > 0) {
      df_filtrado <- df_filtrado %>% filter(scientific_name %in% lista_org)
      cat(paste("Organismos seleccionados: ", paste(lista_org, collapse = ","), "\n"))
    } else {
      cat("Se mantienen todos los organismos.\n")
    }
  }
  
  # FILTRAR POR ESTRATEGIA
  estrategia <- unique(df_filtrado$library_strategy)
  
  if(length(estrategia) > 0){
    lista_est <- dlg_list(choices = estrategia,
                          preselect = estrategia,
                          multiple = TRUE,
                          title = "Selecciona estrategias a MANTENER")$res
    if(length(lista_est) > 0){
      df_filtrado <- df_filtrado %>% filter(library_strategy %in% lista_est)
      cat(paste("Estrategias seleccionadas: ", paste(lista_est, collapse = ","), "\n"))
    } else {
      cat("Se mantienen todos las estrategias.\n")
    }
  }
    
  # FILTRAR POR PALABRA CLAVE
}


construirConsulta <- function() {
  
  query <- dlgInput(message = "Término general (ej: multiple sclerosis):", default = "multiple sclerosis")$res
  if (!length(query)) return(NULL)
  
  datos_raw <- busquedaENA(query, limit = 2000)
  
  datos_filtrados <- filtrarDatos(datos_raw)
  
  if (nrow(datos_filtrados) == 0) {
    message("No se encontraron muestras que cumplan todos los criterios.")
    return(NULL)
  }
  
  print(table(datos_filtrados$Grupo))
  
  print(datatable(datos_filtrados, 
                  caption = "Muestras seleccionadas. Revisa y filtra si es necesario.",
                  options = list(pageLength = 10, scrollX = TRUE),
                  filter = 'top',
                  selection = 'none')) 
  
  return(datos_filtrados)
}


explorarResultado <- function(df) {
  return(invisible(df))
}

