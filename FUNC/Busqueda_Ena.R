busquedaENA <- function(dominio, query, fields, limit = 100) {
  url <- paste0("https://www.ebi.ac.uk/ebisearch/ws/rest/", dominio)
  
  fields <- gsub(",\\s*", ",", fields)
  
  parametros <- list(
    query = query,
    fields = fields,
    format = "json",
    size = 100
  )
  
  print(paste("Realizando la consulta EBI Search (Dominio:", dominio, "):"))
  print(paste("Query:", parametros$query))
  
  respuesta <- GET(url = url, query = parametros)
  stop_for_status(respuesta)
  
  if(http_error(respuesta)) {
    print(paste("URL generada y fallida:", respuesta$url))
  }
  
  contenidoRespuesta <- httr::content(respuesta, "text", encoding = "UTF-8")
  dataJson <- fromJSON(contenidoRespuesta, flatten = TRUE)
  
  if (dataJson$hitCount == 0) {
    print(paste("No se encontraron resultados para esta consulta:", respuesta$url))
  }
  
  cat("ÉXITO:", dataJson$hitCount, "muestras encontradas\n")
  #return(as.data.frame(dataJson$entries))
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
  
  fields <- "sample_accession,study_accession,description,disease"
  
  busquedaENA(dominio = dominio, query = query, fields = fields, limit = limit)
}


explorarResultado <- function(df) {
  if (nrow(df) == 0) {
    cat("El dataframe esta vacío. No hay nada que explorar.\n")
    return(NULL)
  }
  
  datatable(df, options = list(papeLength = 50, scrollX = TRUE))
  
  return(df)
}