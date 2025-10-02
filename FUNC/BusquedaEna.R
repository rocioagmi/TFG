# CAMBIAR ESTO

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
  
  if (status_code(response) == 200) {
    data_text <- httr::content(response,"text", encoding = "UTF-8")
    df <- read_tsv(data_text, show_col_types = FALSE)
    return(df)
  } else if (status_code(response) == 429) {
    stop("Error 429: Demasiadas solicitudes. Espera y reintenta.")
  } else {
    stop("Error en API: ", status_code(response), "\n", httr::content(response, "text", encoding = "UTF-8"))
  }
}


busquedaInteractivaENA <- function() {
  cat("---Búsqueda interactiva en ENA API ---\n")
  
  result_type <- "raw_reads"
  # 'description="16S rRNA" AND (description="multiple slerosis" OR description="MS")'
  query <- 'description="16S rRNA"'
  fields <- "sample_accession, sample_description, description, study_accession"
  limit <- "1000"
  
  cat("\nEjecutando búsqueda...\n")
  resultados <- busquedaEna(result_type, query, fields, limit, format = "tsv")
  
  if (nrow(resultados) == 0) {
    cat("No se encontraron resultados.\n")
    return(NULL)
  }
  
  filename <- paste0("INPUT/ena_resultados_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".tsv")
  write_tsv(resultados, filename)
  cat("Resultados guardados en: ", filename, "\n")
  cat("Número de resultados: ", nrow(resultados), "\n")
  print(head(resultados))
  
  return(resultados)
}