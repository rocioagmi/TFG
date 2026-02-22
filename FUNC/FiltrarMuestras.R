library(dada2)

filtrarMuestras <- function(R1, R2){
  
  if(length(R1) == 0 || length(R2) == 0) stop("No se encontraron archivos de entrada.")
  
  nombres_R1 <- file.path("OUTPUT/FILTRADO", basename(R1))
  nombres_R2 <- file.path("OUTPUT/FILTRADO", basename(R2))
  
  filtrado <- filterAndTrim(R1, nombres_R1, R2, nombres_R2, truncLen = c(230,160),
                trimLeft = c(10, 10), maxN = 0, maxEE = c(2,2), truncQ = 2, 
                rm.phix = TRUE, compress = TRUE, multithread = FALSE)
  
  cat(sprintf("\nFiltrado completado: %d/%d muestras conservaron reads.", 
              sum(filtrado[,2] > 0), nrow(filtrado)))
  
  muestras_vacias <- rownames(filtrado)[filtrado[,2] == 0]
  if (length(muestras_vacias) > 0) {
    warning("Las siguientes muestras quedaron vacías tras el filtrado:\n",
            paste(muestras_vacias, collapse = "\n"))
  }
  
  print(filtrado)
  return(list(R1 = nombres_R1[filtrado[,2] > 0], 
              R2 = nombres_R2[filtrado[,2] > 0],
              stats = filtrado))
  
}
