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


filtradoSR <- function(listadoMuestras){
  for (muestra in listadoMuestras){
    destino <- file.path("INPUT/DATA/FILTRADAS", paste0(basename(muestra), "_filtered.fastq.gz"))
    
    stream <- open(FastqStreamer(muestra))
    on.exit(close(stream))
    
    repeat {
      fq <- yield(stream)
      if (length(fq) == 0)
        break
      
      fq <- fq[nFilter()(fq)] 
      fq <- trimTailw(fq, 2, "4", 2)
      fq <- fq[width(fq) >= 36]
      
      writeFastq(fq, destino, "a")
    }
  }
}