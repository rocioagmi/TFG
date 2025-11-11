busquedaENA <- function(dominio, query, fields, filter = filter, limit = 1000) {
  url <- paste0("https://www.ebi.ac.uk/ebisearch/ws/rest/", dominio)
  fields <- gsub(",\\s*", ",", fields)
  
  print(paste("Realizando la consulta EBI Search (Dominio:", dominio, "):"))
  print(paste("Query:", query))
  print(paste("Filtros:", filter))
  
  todas_las_muestras <- data.frame()
  start <- 0
  
  repeat{
    parametros <- list(
      query = query,
      fields = fields,
      filter = filter,
      format = "json",
      size = limit,
      start = start
    )
    
    print(paste("Descargando página desde start = ", start))
    
    respuesta <- GET(url = url, query = parametros)
    
    if(http_error(respuesta)) {
      print(paste("Error en página:", respuesta$url))
      break
    }
    
    stop_for_status(respuesta)
    
    contenidoRespuesta <- httr::content(respuesta, "text", encoding = "UTF-8")
    dataJson <- fromJSON(contenidoRespuesta, flatten = TRUE, simplifyDataFrame = TRUE)
    
    if(dataJson$hitCount == 0 || length(dataJson$entries) == 0){
      print(paste("No se encontraron resultados para esta consulta:", respuesta$url))
      break
    }
    
    temporalDF <- as.data.frame(dataJson$entries)
    
    todas_las_muestras <- rbind(todas_las_muestras, temporalDF)
    
    if(nrow(temporalDF) < limit) break
    start <- start + limit
    Sys.sleep(0.3)
  }
  
  cat("ÉXITO:", dataJson$hitCount, "muestras encontradas\n")
  return(todas_las_muestras)
}



construirConsulta <- function(limit = 1000) {
  dominios_validos <- c("nucleotideSequences", "project", "sra", "biosamples", 
                        "sra-study", "sra-sample", "sra-run", "sra-experiment")
  
  dominio <- dlgInput(message = "Introduzca el dominio de búsqueda:")$res
  if (!is.character(dominio) || length(dominio) == 0) {
    print("Operación (dominio) cancelada por el usuario.")
    return(NULL)
  } else if (!dominio %in% dominios_validos){
    print("Dominio no válido.")
    return(NULL)
  }
  
  query <- dlgInput(message = "Ingrese un término de búsqueda:")$res
  if (!is.character(query) || length(query) == 0) {
    print("Operación (query) cancelada por el usuario.")
    return(NULL)
  }
  
  filter <- dlgInput(message = "Ingrese un término para filtrar la búsqueda:")$res
  if (!is.character(filter) || length(filter) == 0) {
    print("Operación (filter) cancelada por el usuario.")
    return(NULL)
  }
  
  dominio <- trimws(dominio)
  query <- trimws(query)
  filter <- trimws(filter)
  
  fields <- "description"
  
  muestrasDF <- busquedaENA(dominio = dominio, query = query, fields = fields, filter = filter, limit = limit)
  return(muestrasDF)
}


explorarResultado <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    cat("El dataframe está vacío.\n")
    return(NULL)
  }
  
  print(datatable(df, options = list(pageLength = 50, scrollX = TRUE)))
  return(invisible(df))
}