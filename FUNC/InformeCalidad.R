informeCalidadSR <- function(directorioMuestras){
  qaSummary <- qa(directorioMuestras, type = "fastq")
  browseURL(report(qaSummary)) 
}

informeCalidadD <- function(R1, R2){
  plotQualityProfile(R1)
  plotQualityProfile(R2)
}