informeCalidad <- function(directorioMuestras){
  qaSummary <- qa(directorioMuestras, type = "fastq")
  report(qaSummary, dest = "OUTPUT/REPORT", type = "html")
  browseURL(file.path("OUTPUT/REPORT", "index.html")) 
}

