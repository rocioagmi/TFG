informeCalidadSR <- function(directorioMuestras){
  qaSummary <- qa(directorioMuestras, type = "fastq")
  browseURL(report(qaSummary)) 
}

informeCalidadR <- function(directorioMuestras){
  qaSummary <- rqc(path = directorioMuestras, pattern = "fastq.gz", openBrowser = TRUE)
}