informeCalidadSR <- function(directorioMuestras){
  qaSummary <- qa(directorioMuestras, type = "fastq")
  browseURL(report(qaSummary)) 
}

