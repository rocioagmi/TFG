library(dada2)

procesarSecuencias <- function(union) {
  
  # CONSTRUIR LA TABLA DE SECUENCIAS
  message("Construyendo tabla de secuencias")
  seqtab <- makeSequenceTable(union)
  
  message(sprintf("Dimensiones: %d muestras x %d ASV's", nrow(seqtab), ncol(seqtab)))
  message("Distribución de longitudes de amplicón:")
  print(table(nchar(getSequences(seqtab))))
  
  saveRDS(seqtab, "OUTPUT/RDS/seqtab.Rds")

  # ELIMINAR QUIMERAS
  message("Eliminando quimeras")
  tabSinQuim <- removeBimeraDenovo(seqtab, method = "consensus", 
                                   multithread = TRUE, verbose = TRUE)
  pct_limpias <- round(sum(tabSinQuim) / sum(seqtab) * 100, 2)
  message(sprintf("Lecturas conservadas tras eliminar quimeras: %s%%", pct_limpias))
  message(sprintf("ASV's: %d  - %d", ncol(seqtab), ncol(tabSinQuim)))
  
  saveRDS(tabSinQuim, "OUTPUT/RDS/tabSinQuim.Rds")
  message("Tabla de secuencias sin quimeras guardada en OUTPUT/RDS/tabSinQuim.Rds")
  
  return(tabSinQuim)
}