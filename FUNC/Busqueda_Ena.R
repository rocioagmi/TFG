busquedaENA <- function(dominio, query, size = 10, fields, start = 0) {
  url <- paste0("https://www.ebi.ac.uk/ebisearch/rest/", dominio,"/search")
  
  parametros <- list(
    query = query,
    format = "json",
    size = size,
    fields = fields,
    start = start
  )
  
}

construirBusqueda <- function() {
  dominio <- dlgInput(message = "Introduzca el dominio de búsqueda:", default = "ena_sample, ena_run, ena_study")
  query <- dlgInput(message = "Ingresa un término de búsqueda:")
  fields <- "ID, Name, Description"
}