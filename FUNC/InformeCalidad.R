informeCalidad <- function(directorioMuestras){
  qaSummary <- qa(directorioMuestras, type = "fastq")
  
  timestamp <- format(Sys.time(), "%d-%m-%Y_%H:%M")
  directorio_salida <- file.path("OUTPUT/REPORT", paste0("InformeCalidad_", timestamp))
  
  report(qaSummary, dest = directorio_salida, type = "html")
  browseURL(file.path(directorio_salida, "index.html")) 
}

