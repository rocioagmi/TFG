informeCalidad <- function(directorioMuestras){
  qaSummary <- qa(directorioMuestras, type = "fastq")
  report(qaSummary, dest = "INPUT/DATA/REPORT", type = "html")
  browseURL(file.path("INPUT/DATA/REPORT", "index.html")) 
}

