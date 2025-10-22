busquedaENA <- function(dominio, query, fields, limit = 100) {
  url <- paste0("https://www.ebi.ac.uk/ebisearch/ws/rest/", dominio)
  
  parametros <- list(
    query = query,
    fields = fields,
    format = "json",
    size = limit
  )
  
  print(paste("Realizando la consulta EBI Search (Dominio:", dominio, "):"))
  print(paste("Query:", query))
  
  respuesta <- GET(url = url, query = parametros)
  
  stop_for_status(respuesta, task = "consultar API de EBI")
  
  contenidoRespuesta <- content(respuesta, "text", encoding = "UTF-8")
  
  dataJson <- fromJSON(contenidoRespuesta, flatten = TRUE)
  
  if (dataJson$hitCount == 0) {
    print("No se encontraron resultados para esta consulta.")
  }
  
  resultadosDF <- dataJson$entries
  
  return(resultadosDF)
}

construirConsulta <- function(limit = 1000) {
  dominio <- dlgInput(message = "Introduzca el dominio de búsqueda (ena_sample, ena_run, ena_study):", default = "ena_sample")
  if (!is.character(dominio$res) || length(dominio$res) == 0) {
    print("Operación cancelada por el usuario (dominio).")
    return(NULL)
  }
  
  query <- dlgInput(message = "Ingrese un término de búsqueda:")
  if (!is.character(query$res) || length(query$res) == 0) {
    print("Operación cancelada por el usuario (query).")
    return(NULL)
  }
  
  dominio <- trimws(dominio$res)
  query <- trimws(query$res)
  
  fields <- "id, acc, description, ena_study, ena_run, sample_title"
  
  muestrasDF <- busquedaENA(dominio = dominio, query = query, fields = fields, limit = limit)
  
  return(muestrasDF)
}