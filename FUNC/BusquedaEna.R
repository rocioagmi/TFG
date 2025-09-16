busquedaEna <- function(result_type, query, fields, limit, format) {
  base_url <- "https://www.ebi.ac.uk/ena/portal/api/search"
  params <- list(
    result = result_type,
    query = query,
    fields = fields,
    limit = limit,
    format = format
  )
  
  response <- POST(base_url, body = params, encode = "form")
  data_text <- content(response, "text", encoding = "UTF-8")
  df <- read_tsv(data_text, show_col_types = FALSE)
  
  return(df)
}


busquedaInteractivaENA <- function() {
  cat("---Búsqueda interactiva en ENA API ---\n")
  
  cat("\n1. Consulta (ej:description = '16S rRNA' AND description = 'multiple sclerosis'): ")
  query <- readline()
  
  cat("\nEjecutando búsqueda...\n")
  resultados <- busquedaEna(result_type = study, query, fields = "study_accesion, study_title, description", limit = 1000, format = "tsv")
  
  if (nrow(resultados) == 0) {
    cat("No se encontraron resultados.\n")
    return(NULL)
  }
  
  filename <- paste0("INPUT/ena_resultados_", Sys.time(), ".tsv")
  write_tsv(resultados, filename)
  cat("Resultados guardados en: ", filename, "\n")
  cat("Número de resultados: ", nrow(resultados), "\n")
  print(head(resultados))
  
  return(resultados)
}