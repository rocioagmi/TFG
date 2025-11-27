library(dada2)

filtrarMuestras <- function(R1, R2){
  
  nombres_R1 <- file.path("OUTPUT/FILTRADO", basename(R1))
  nombres_R2 <- file.path("OUTPUT/FILTRADO", basename(R2))
  
  filterAndTrim(R1, nombres_R1, R2, nombres_R2, truncLen = c(230,160),
                trimLeft = c(10, 10), maxN = 0, maxEE = c(2,5), truncQ = 2, 
                rm.phix = TRUE, compress = TRUE, multithread = FALSE)
  
  return(list(nombres_R1 = nombres_R1, nombres_R2 = nombres_R2))
}
