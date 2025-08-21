informeCalidad <- function(directorioMuestras){
  qaSummary <- qa(directorioMuestras, type = "fastq")
  direcrtorio <- mkdir("OUTPUT/REPORT")
  report(qaSummary, dest = "OUTPUT/REPORT", type = "html")
  browseURL(file.path("OUTPUT/REPORT", "index.html")) 
}

