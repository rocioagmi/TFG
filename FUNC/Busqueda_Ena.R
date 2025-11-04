busquedaENA <- function(dominio, query, fields, limit = 100) {
  url <- paste0("https://www.ebi.ac.uk/ebisearch/ws/rest/", dominio)
  
  query <- gsub('"', '', query)
  fields <- gsub(",\\s*", ",", fields)
  
  parametros <- list(
    query = query,
    fields = fields,
    format = "json",
    size = limit
  )
  
  print(paste("Realizando la consulta EBI Search (Dominio:", dominio, "):"))
  print(paste("Query:", query))
  
  respuesta <- GET(url = url, query = parametros)
  
  if(http_error(respuesta)) {
    print(paste("URL generada y fallida:", respuesta$url))
  }
  
  stop_for_status(respuesta, task = "consultar API de EBI")
  
  contenidoRespuesta <- httr::content(respuesta, "text", encoding = "UTF-8")
  
  dataJson <- fromJSON(contenidoRespuesta, flatten = TRUE)
  
  if (dataJson$hitCount == 0) {
    print("No se encontraron resultados para esta consulta.")
  }
  
  resultadosDF <- dataJson$entries
  return(resultadosDF)
}

construirConsulta <- function(limit = 1000) {
  dominio <- dlgInput(message = "Introduzca el dominio de búsqueda (biosamples, embl):")$res
  if (!is.character(dominio) || length(dominio) == 0) {
    print("Operación cancelada por el usuario (dominio).")
    return(NULL)
  }
  
  query <- dlgInput(message = "Ingrese un término de búsqueda:")$res
  if (!is.character(query) || length(query) == 0) {
    print("Operación cancelada por el usuario (query).")
    return(NULL)
  }
  
  dominio <- trimws(dominio)
  query <- trimws(query)
  
  fields <- "id,acc,description,sample_title"
  
  muestrasDF <- busquedaENA(dominio = dominio, query = query, fields = fields, limit = limit)
  
  return(muestrasDF)
}