busquedaENA <- function(result_type, query, fields, limit, format){
  url <- "https://www.ebi.ac.uk/ena/portal/api/search"
  parametros <- list(
    result = result_type,
    query = query,
    fields = fields,
    limit = limit,
    format = format
  )
  
  response <- POST(url, body = parametros, encode = "form")
  
  if(status_code(response) == 200){
    data_text <- httr::content(response, "text", encoding = "UTF-8")
    df <- fromJSON(data_text)
    return(df)
  } else if (status_code(response) == 429) {
    stop("Error 429: Demasiadas solicitudes. Espera y reintenta.")
  } else {
    stop("Error en API: ", status_code(response))
  }
}

construirConsultaENA <- function(){
  
  result_type <- dlgInput("Tipo de resultado (sample, study):")$res
  
  input <- dlgInput("Introduce un término de búsqueda:")$res
  query <- paste0('description="', input,'"')
  
  if(result_type == "sample"){
    fields <- "sample_accession, study_accession, sample_title, sample_alias, description, sample_description"
  } else if(result_type == "study"){
    fields <- "study_accession, study_name, study_title, description, study_description"
  } 
  
  cat("\nEjecutando búsqueda...\n")
  resultados <- busquedaENA(result_type, query, fields, limit = 1000, format = "json")
  
  if (nrow(resultados) == 0) {
    cat("No se encontraron resultados.\n")
    return(NULL)
  } else {
    head(resultados)
  }
}