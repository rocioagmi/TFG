flujoTrabajoDada <- function(R1, R2){
  # Error Rates
  err_R1 <- learnErrors(R1)
  err_R2 <- learnErrors(R2)
  
  # Derreplicar
  derep_R1 <- derepFastq(R1, verbose = TRUE)
  derep_R2 <- derepFastq(R2, verbose = TRUE)
  
  nombres <- sapply(strsplit(basename(R1), "_"), `[`, 1)
  
  names(derep_R1) <- nombres
  names(derep_R2) <- nombres
  
  # Aplica algoritmo DADA2
  dadaR1 <- dada(derep_R1, err = err_R1, multithread = FALSE)
  dadaR2 <- dada(derep_R2, err = err_R2, multithread = FALSE)
  
  # Junta las secuencias R1 y R2
  union <- mergePairs(dadaR1, derep_R1, dadaR2, derep_R2, verbose = TRUE, justConcatenate = TRUE)

  return(union)
} 
