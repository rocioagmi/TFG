busquedaENA <- function(dominio, query, fields, limit = 1000) {
  url <- paste0("https://www.ebi.ac.uk/ebisearch/ws/rest/", dominio)
  
  # query <- gsub('"', '', query)
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
    print(paste("No se encontraron resultados para esta consulta:", respuesta$url))
  }
  
  resultadosDF <- as.data.frame(dataJson$entries)
  return(resultadosDF)
}


construirConsulta <- function(limit = 1000) {
  dominio <- dlgInput(message = "Introduzca el dominio de búsqueda (biosamples, sra-sample):")$res
  if (!is.character(dominio) || length(dominio) == 0) {
    print("Operación (dominio) cancelada por el usuario.")
    return(NULL)
  }
  
  query <- dlgInput(message = "Ingrese un término de búsqueda:")$res
  if (!is.character(query) || length(query) == 0) {
    print("Operación (query) cancelada por el usuario.")
    return(NULL)
  }
  
  dominio <- trimws(dominio)
  query <- trimws(query)
  
  # "id,acc,description,sample_title"
  fields <- "sample_accession,study_accession,description,disease,subject"
  
  muestrasDF <- busquedaENA(dominio = dominio, query = query, fields = fields, limit = limit)
  
  return(muestrasDF)
}