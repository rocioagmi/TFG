
flujoTrabajoDada <- function(R1, R2, estudio){
  
  require(dada2)
  
  # Error Rates
  err_R1 <- learnErrors(R1, multithread =  FALSE)
  err_R2 <- learnErrors(R2, multithread =  FALSE)
  
  # Derreplicar
  derep_R1 <- derepFastq(R1, verbose = TRUE)
  derep_R2 <- derepFastq(R2, verbose = TRUE)
  
  nombres <- sub("^((HC|MS)_)?([^_]+).*$", "\\1\\3", basename(R1))
  
  names(derep_R1) <- nombres
  names(derep_R2) <- nombres
  
  # Aplica algoritmo DADA2
  dadaR1 <- dada(derep_R1, err = err_R1, pool = "pseudo", multithread = FALSE)
  dadaR2 <- dada(derep_R2, err = err_R2, pool = "pseudo", multithread = FALSE)
  
  # Junta las secuencias R1 y R2
  union <- mergePairs(dadaR1, derep_R1, dadaR2, derep_R2, verbose = TRUE)
  
  tasa_union <- sapply(seq_along(union), function(i){
    sum(union[[i]]$abundance) / sum(getUniques(dadaR1[[i]]))
  })
  
  directorio_rds <- file.path("OUTPUT", estudio, "RDS")
  dir.create(directorio_rds, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(list(err_R1 = err_R1, err_R2 = err_R2), file.path(directorio_rds, "errores.Rds"))
  saveRDS(union, file.path(directorio_rds, "union.Rds"))
  
  cat(sprintf("Tasa de unión R1/R2 - %s: min %.1f%%, media %.1f%%, max %.1f&&\n",
              estudio, min(tasa_union)*100, mean(tasa_union)*100, max(tasa_union)*100))
  
  return(union)
} 
