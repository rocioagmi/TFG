# FILTRADO DE LAS MUESTRAS POR CALIDAD Q = 20

library(Biostrings)
library(ShortRead)
library(dada2)

filtrarMuestras <- function(R1, R2){
  
  nombres_R1 <- file.path("OUTPUT/FILTRADO", basename(R1))
  nombres_R2 <- file.path("OUTPUT/FILTRADO", basename(R2))
  
  filterAndTrim(R1, nombres_R1, R2, nombres_R2, truncLen = c(230,160),
                trimLeft = c(10, 10), maxN = 0, maxEE = c(2,5), truncQ = 2, 
                rm.phix = TRUE, compress = TRUE, multithread = FALSE)
  
  return(list(nombres_R1 = nombres_R1, nombres_R2 = nombres_R2))
}



filtradoSR <- function(listadoMuestras){
  for (muestra in listadoMuestras){
    
    stream <- open(FastqStreamer(muestra))
    on.exit(close(stream))
    
    repeat {
      fq <- yield(stream)
      if (length(fq) == 0)
        break
      
      fq <- fq[nFilter()(fq)] 
      fq <- trimTailw(fq, 2, "4", 2)
      fq <- fq[width(fq) >= 36]
      
      destino <- file.path("OUTPUT/FILTRADOSR", basename(muestra))
      
      if (!file.exists(destino)) {
        writeFastq(fq, destino, "w")
      } else {
        writeFastq(fq, destino, "a")
      }
    }
  }
}