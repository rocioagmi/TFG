# FILTRADO DE LAS MUESTRAS POR CALIDAD Q = 20

library(Biostrings)
library(ShortRead)
library(dada2)

filtrarMuestras <- function(R1, R2){
  nombres_R1 <- paste("INPUT/DATA", strsplit(basename(R1), "/M"),sep = "/FILTRADAS/")
  nombres_R2 <- paste("INPUT/DATA", strsplit(basename(R2), "/M"),sep = "/FILTRADAS/")
  
  filterAndTrim(R1, nombres_R1, R2, nombres_R2, truncLen = c(240,200),
                 trimLeft = c(15, 15), maxN = 0, maxEE = c(2,2), truncQ = 10, 
                 rm.phix = TRUE, compress = TRUE, multithread = FALSE)
  
  return(list(nombres_R1 = nombres_R1, nombres_R2 = nombres_R2))
}
